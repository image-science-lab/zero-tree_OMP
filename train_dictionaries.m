%% Master script for training dictionaries for any signal class.

% Parameters and data is loaded from the parent script.

% Create the directory to save output.
mkdir(save_dir);

%% Train a small dictionary.
params.data = data;
params.Tdata = k0;
params.dictsize = N1;
params.iternum = niters;
params.exact = 0;

fprintf('Training small dictionary\n');
[D, X, ~, ~] = box_ksvd(params);

% Save data.
save([save_dir '/small_dict.mat'], 'D', 'k0');

%% Train a large dictionary.
clear params;
params.data = data;
params.Tdata = k1;
params.dictsize = N1*l0_resolve;
params.iternum = niters;
params.exact = 0;

fprintf('Training large` dictionary\n');
[D, X, ~, ~] = box_ksvd(params);

% Save data.
save([save_dir 'data/large_dict.mat'], 'D', 'k1');

%% Train a low resolution dictionary.
clear params;
params.data = W*data;
params.Tdata = k1;
params.dictsize = N1;
params.iternum = niters;
params.exact = 0;

fprintf('Training low resolution dictionary\n');
[D1, X1, ~, ~] = box_ksvd(params);

% Train a high resolution dictionary.
clear params;
params.data = data;
params.Tdata = k0;
params.dictsize = N1*l0_resolve;
params.iternum = niters;
params.l1_dict = D1;
params.l1_sparse = X1;
params.initdict = get_init_dict(params.data, D1, X1, N1, l0_resolve, U, 'col');
params.l0_resolve = l0_resolve;
params.W = W;
params.U = U;

fprintf('Training high resolution dictionary\n');
[D0, ~, ~, ~] = box_ksvd_scale(params);

% Save
save([save_dir 'data/two_scale_dict.mat'], 'D0', 'D1', 'k0', 'k1', 'l0_resolve', ...
	'W', 'U');
