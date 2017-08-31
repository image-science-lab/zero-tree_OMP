function [D,Gamma,err,gerr] = box_ksvd_scale(params,varargin)
% Locate all the necessary information in box_ksvd.m . This section will 
% contain further information for this script.
% Extra information needed in params:
%	'l1_dict' -- Low resolution dictionary.
%	'l1_sparse' -- Low resolution sparse coded signal.
%	'l0_resolve' -- Number of high resolution atoms per each low resoluttion atom.

global CODE_SPARSITY CODE_ERROR codemode
global MEM_LOW MEM_NORMAL MEM_HIGH memusage
global ompfunc ompparams exactsvd

CODE_SPARSITY = 1;
CODE_ERROR = 2;

MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;


%%%%% parse input parameters %%%%%
data = params.data;
D1 = params.l1_dict;
X1 = params.l1_sparse;
l0_resolve = params.l0_resolve;
W = params.W;
U = params.U;

ompparams = {'checkdict','off'};

% coding mode %

if (isfield(params,'codemode'))
  switch lower(params.codemode)
    case 'sparsity'
      codemode = CODE_SPARSITY;
      thresh = params.Tdata;
    case 'error'
      codemode = CODE_ERROR;
      thresh = params.Edata;
    otherwise
      error('Invalid coding mode specified');
  end
elseif (isfield(params,'Tdata'))
  codemode = CODE_SPARSITY;
  thresh = params.Tdata;
elseif (isfield(params,'Edata'))
  codemode = CODE_ERROR;
  thresh = params.Edata;

else
  error('Data sparse-coding target not specified');
end


% max number of atoms %

if (codemode==CODE_ERROR && isfield(params,'maxatoms'))
  ompparams{end+1} = 'maxatoms';
  ompparams{end+1} = params.maxatoms;
end


% memory usage %

if (isfield(params,'memusage'))
  switch lower(params.memusage)
    case 'low'
      memusage = MEM_LOW;
    case 'normal'
      memusage = MEM_NORMAL;
    case 'high'
      memusage = MEM_HIGH;
    otherwise
      error('Invalid memory usage mode');
  end
else
  memusage = MEM_NORMAL;
end


% iteration count %

if (isfield(params,'iternum'))
  iternum = params.iternum;
else
  iternum = 10;
end


% omp function %

if (codemode == CODE_SPARSITY)
  ompfunc = @omp;
else
  ompfunc = @omp2;
end


% status messages %

printiter = 0;
printreplaced = 0;
printerr = 0;
printgerr = 0;

verbose = 't';
msgdelta = -1;

for i = 1:length(varargin)
  if (ischar(varargin{i}))
    verbose = varargin{i};
  elseif (isnumeric(varargin{i}))
    msgdelta = varargin{i};
  else
    error('Invalid call syntax');
  end
end

for i = 1:length(verbose)
  switch lower(verbose(i))
    case 'i'
      printiter = 1;
    case 'r'
      printiter = 1;
      printreplaced = 1;
    case 't'
      printiter = 1;
      printerr = 1;
      if (isfield(params,'testdata'))
        printgerr = 1;
      end
  end
end

if (msgdelta<=0 || isempty(verbose))
  msgdelta = -1; 
end

ompparams{end+1} = 'messages';
ompparams{end+1} = msgdelta;



% compute error flag %

comperr = (nargout>=3 || printerr);


% validation flag %

testgen = 0;
if (isfield(params,'testdata'))
  testdata = params.testdata;
  if (nargout>=4 || printgerr)
    testgen = 1;
  end
end


% data norms %

XtX = []; XtXg = [];
if (codemode==CODE_ERROR && memusage==MEM_HIGH)
  XtX = colnorms_squared(data);
  if (testgen)
    XtXg = colnorms_squared(testdata);
  end
end


% mutual incoherence limit %

if (isfield(params,'muthresh'))
  muthresh = params.muthresh;
else
  muthresh = 0.99;
end
if (muthresh < 0)
  error('invalid muthresh value, must be non-negative');
end


% exact svd computation %

exactsvd = 0;
if (isfield(params,'exact') && params.exact~=0)
  exactsvd = 1;
end


% determine dictionary size %

if (isfield(params,'initdict'))
  if (any(size(params.initdict)==1) && all(iswhole(params.initdict(:))))
    dictsize = length(params.initdict);
  else
    dictsize = size(params.initdict,2);
  end
end
if (isfield(params,'dictsize'))    % this superceedes the size determined by initdict
  dictsize = params.dictsize;
end

if (size(data,2) < dictsize)
  error('Number of training signals is smaller than number of atoms to train');
end


% initialize the dictionary %

if (isfield(params,'initdict'))
  if (any(size(params.initdict)==1) && all(iswhole(params.initdict(:))))
    D = data(:,params.initdict(1:dictsize));
  else
    if (size(params.initdict,1)~=size(data,1) || size(params.initdict,2)<dictsize)
      error('Invalid initial dictionary');
    end
    D = params.initdict(:,1:dictsize);
  end
else
  data_ids = find(colnorms_squared(data) > 1e-6);   % ensure no zero data elements are chosen
  perm = randperm(length(data_ids));
  D = data(:,data_ids(perm(1:dictsize)));
end


% normalize the dictionary %

D = normcols(D);

err = zeros(1,iternum);
gerr = zeros(1,iternum);

if (codemode == CODE_SPARSITY)
  errstr = 'RMSE';
else
  errstr = 'mean atomnum';
end



%%%%%%%%%%%%%%%%%  main loop  %%%%%%%%%%%%%%%%%


for iter = 1:iternum
  
  G = [];
  if (memusage >= MEM_NORMAL)
    G = D'*D;
  end
  
  
  %%%%%  sparse coding  %%%%%
  
  %Gamma = sparsecode(data,D,XtX,G,thresh);
  Gamma = constrained_sparsecode(data, D, X1, l0_resolve, thresh);
  
  
  %%%%%  dictionary update  %%%%%
  
  replaced_atoms = zeros(1,dictsize);  % mark each atom replaced by optimize_atom
  
  unused_sigs = 1:size(data,2);  % tracks the signals that were used to replace "dead" atoms.
                                 % makes sure the same signal is not selected twice
  
  p = randperm(dictsize);
  tid = timerinit('updating atoms', dictsize);
  for j = 1:dictsize
	%[D(:,p(j)),gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom(data,D,p(j),Gamma,unused_sigs,replaced_atoms, l0_resolve, X1, D1, R);
	[D(:, p(j)), gamma_j, data_indices] = optimize_atom(data, D, p(j), Gamma, l0_resolve, X1, D1, W, U);
    Gamma(p(j),data_indices) = gamma_j;
    if (msgdelta>0)
      timereta(tid, j, msgdelta);
    end
  end
  if (msgdelta>0)
    printf('updating atoms: iteration %d/%d', dictsize, dictsize);
  end
  
  
  %%%%%  compute error  %%%%%
  
  if (comperr)
    err(iter) = compute_err(D,Gamma,data);
  end
  if (testgen)
    if (memusage >= MEM_NORMAL)
      G = D'*D;
    end
    GammaG = sparsecode(testdata,D,XtXg,G,thresh);
    gerr(iter) = compute_err(D,GammaG,testdata);
  end
   
  
  %%%%%  clear dictionary  %%%%%
  
  [D,cleared_atoms] = cleardict(D,Gamma,data,muthresh,unused_sigs,replaced_atoms);
  
  
  %%%%%  print info  %%%%%
  
  info = sprintf('Iteration %d / %d complete', iter, iternum);
  if (printerr)
    info = sprintf('%s, %s = %.4g', info, errstr, err(iter));
  end
  if (printgerr)
    info = sprintf('%s, test %s = %.4g', info, errstr, gerr(iter));
  end
  if (printreplaced)
    info = sprintf('%s, replaced %d atoms', info, sum(replaced_atoms) + cleared_atoms);
  end
  
  if (printiter)
    disp(info);
    if (msgdelta>0), disp(' '); end
  end
  
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            optimize_atom             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%function [atom,gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom(X,D,j,Gamma,unused_sigs,replaced_atoms, l0_resolve, X1, D1, R)
function [atom, gamma_j, data_indices] = optimize_atom(X, D, j, Gamma, l0_resolve, X1, D1, W, U)
global exactsvd

% data samples which use the atom, and the corresponding nonzero
% coefficients in Gamma
[gamma_j, data_indices] = sprow(Gamma, j);

% Modify the dictionary replacement algorithm a bit to accommodate for the
% block nature of the high resolution dictionary.
if (length(data_indices) < 1)
	fprintf('|');
	if true
		% Get the index of the low resolution dictionary atom.
		atom_idx = ceil(j/l0_resolve);
		% Get all the low resolution sparse vectors where this dictionary atom
		% has contribution.
		member_set = find(X1(atom_idx, :) ~= 0);
		X1_members = X1(:, member_set);
		% We want the contribution from only this particular atom, so subtract
		% all other components.
		X1_members(atom_idx, :) = 0;
		Y_non_members = D1*X1_members;
		Y_members = X(:, member_set) - U*Y_non_members;
    	perm = randperm(size(Y_members, 2));
        if ~isempty(perm)
            atom = Y_members(:, perm(1));
        else
            error = sum((X - D*Gamma).^2, 1);
            [~, max_idx] = max(error);
            atom = X(:, max_idx);
        end
  		atom = atom./norm(atom);
  		gamma_j = zeros(size(gamma_j));
	end
	% Tailored for HF dictionary training.
	if false
		% Find the atoms where this block makes non-zero contribution.
		atom_idx = floor(j/l0_resolve);
		member_set = find(X1(atom_idx, :) ~= 0);
		Y_members = X(:, member_set);
		% Get the residue for this block.
		Y_residue = Y_members - D*Gamma(:, member_set);
		res_norm = diag(Y_residue'*Y_residue);
		[~, max_idx] = max(res_norm);
		atom = Y_members(:, max_idx);
		atom = atom./norm(atom);
		gamma_j = zeros(size(gamma_j));
	end
  return;
end

smallGamma = Gamma(:,data_indices);
Dj = D(:,j);

if (exactsvd)

  [atom,s,gamma_j] = svds(X(:,data_indices) - D*smallGamma + Dj*gamma_j, 1);
  gamma_j = s*gamma_j;
  
else
  
  atom = collincomb(X,data_indices,gamma_j') - D*(smallGamma*gamma_j') + Dj*(gamma_j*gamma_j');
  atom = atom/norm(atom);
  gamma_j = rowlincomb(atom,X,1:size(X,1),data_indices) - (atom'*D)*smallGamma + (atom'*Dj)*gamma_j;

end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             sparsecode               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Gamma = sparsecode(data,D,XtX,G,thresh)

global CODE_SPARSITY codemode
global MEM_HIGH memusage
global ompfunc ompparams

if (memusage < MEM_HIGH)
  Gamma = ompfunc(D,data,G,thresh,ompparams{:});
  
else  % memusage is high
  
  if (codemode == CODE_SPARSITY)
    Gamma = ompfunc(D'*data,G,thresh,ompparams{:});
    
  else
    Gamma = ompfunc(D'*data,XtX,G,thresh,ompparams{:});
  end
  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			constrained_sparsecode	   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Gamma = constrained_sparsecode(data, D, X1, l0_resolve, thresh)
	% Global constants
	global CODE_SPARSITY codemode;
	global ompfunc ompparams;
	
	% since the constrain for each data vector is separate, we 
	% have to iterate through each data signal.
    % thresh = k0;
    
	%Gamma = zeros(size(D,2), size(data, 2));
    
    % Store data in sparse format.
    row_idx = zeros(thresh*size(data, 2), 1);
    col_idx = zeros(thresh*size(data, 2), 1);
    values = zeros(thresh*size(data, 2), 1);
    iidx = 1;
	for idx = 1:1:size(data, 2)
		y = data(:, idx);
		x1 = X1(:, idx);
		cset = find(x1 ~= 0);
        if isempty(cset)
            cset = (1:length(x1))';
        end
		single_idx = 0:l0_resolve-1;
		repeated_idx = repmat(single_idx, 1, length(cset(:)));
		temp = repmat((cset-1)*l0_resolve+1, 1, l0_resolve)';
		repeated_cset = temp(:)' + repeated_idx;
		D_sub = D(:, repeated_cset);
		G = D_sub'*D_sub;
		X0_sub = omp(D_sub'*y, G, thresh, ompparams{:});
		%Gamma(repeated_idx + repeated_cset, idx) = X0_sub;
		%keyboard();
		xsup = repeated_cset(X0_sub ~= 0);
		xlen = length(xsup);

        row_idx(iidx:iidx+xlen-1) = xsup;
        col_idx(iidx:iidx+xlen-1) = idx*ones(xlen, 1);
        values(iidx:iidx+xlen-1) = X0_sub(X0_sub ~= 0);
		iidx = iidx + xlen;
	end
	% Save the matrix as sparse matrix.
	%keyboard();
    Gamma = sparse(row_idx(1:iidx-1),...
   			       col_idx(1:iidx-1),...
				   values(1:iidx-1),...
				   size(D, 2),...
				   size(data, 2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             compute_err              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err = compute_err(D,Gamma,data)
  
global CODE_SPARSITY codemode

if (codemode == CODE_SPARSITY)
  %err = sqrt(sum(reperror2(data,D,Gamma))/numel(data));
  err = -20*log10(norm(data - D*Gamma, 'fro')/norm(data, 'fro'));
else
  err = nnz(Gamma)/size(data,2);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           cleardict                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [D,cleared_atoms] = cleardict(D,Gamma,X,muthresh,unused_sigs,replaced_atoms)

use_thresh = 4;  % at least this number of samples must use the atom to be kept

dictsize = size(D,2);

% compute error in blocks to conserve memory
err = zeros(1,size(X,2));
blocks = [1:3000:size(X,2) size(X,2)+1];
for i = 1:length(blocks)-1
  err(blocks(i):blocks(i+1)-1) = sum((X(:,blocks(i):blocks(i+1)-1)-D*Gamma(:,blocks(i):blocks(i+1)-1)).^2);
end

cleared_atoms = 0;
usecount = sum(abs(Gamma)>1e-7, 2);

for j = 1:dictsize
  
  % compute G(:,j)
  Gj = D'*D(:,j);
  Gj(j) = 0;
  
  % replace atom
  if ( (max(Gj.^2)>muthresh^2 || usecount(j)<use_thresh) && ~replaced_atoms(j) )
    [y,i] = max(err(unused_sigs));
    D(:,j) = X(:,unused_sigs(i)) / norm(X(:,unused_sigs(i)));
    unused_sigs = unused_sigs([1:i-1,i+1:end]);
    cleared_atoms = cleared_atoms+1;
  end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            misc functions            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err2 = reperror2(X,D,Gamma)

% compute in blocks to conserve memory
err2 = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  err2(blockids) = sum((X(:,blockids) - D*Gamma(:,blockids)).^2);
end

end


function Y = colnorms_squared(X)

% compute in blocks to conserve memory
Y = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  Y(blockids) = sum(X(:,blockids).^2);
end

end
