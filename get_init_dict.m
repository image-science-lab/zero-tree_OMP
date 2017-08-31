function init_dict = get_init_dict(Y_train, D1, X1, N, l0_resolve, U, method)
	% Function to get an initial dictionary for training a multiscale
	% hig resolution dictionary.
	%
	% Inputs:
	%	Y_train: Training data
	%	D1: Low resolution dictionary.
	%	X1: Low resolution sparse representation.
	%	N:  Number of low resolution atoms.
	%	l0_resolve: Number of high resolution atoms per each low resolution
	%		atoms.
	%	U: Upscale operator.
    %   method: 'kmeans' or 'randcol'
	%
	% Output:
	%	init_dict: Initialized dictionary.

	D0 = zeros(size(Y_train,1), N*l0_resolve);
	for atom_idx = 1:1:N
    	% Get all those vectors where the atom has non zero contribution.
    	member_set = find(X1(atom_idx, :) ~= 0);
    	X1_members = X1(:, member_set);
    	% We want the contribution only from this particular atom, so subtract
    	% all other components.
    	X1_members(atom_idx, :) = 0;
    	Y_non_members = D1*X1_members;
    	Y_members = Y_train(:, member_set) - U*Y_non_members;
        if size(Y_members, 2) == 0
            residue = Y_train - U*D1*X1;
            residue_norms = sum(residue.^2);
            [~, sorted_idx] = sort(residue_norms, 'descend');
            % Extremely tricky. Never knew that such a thing can happen.
            D0(:, (atom_idx-1)*l0_resolve+1:(atom_idx)*l0_resolve) = ...
                Y_train(:, sorted_idx(1:l0_resolve));
        elseif size(Y_members, 2) < l0_resolve
            fprintf('%d -- %d\n', atom_idx, size(Y_members, 2));
        	rand_idx = randi([1, size(Y_members, 2)], 1, l0_resolve);
        	D0(:, (atom_idx-1)*l0_resolve+1:(atom_idx)*l0_resolve) = Y_members(:, rand_idx);
        else
            if strcmp(method, 'kmeans')
                % Now get l0_resolve random atoms from this
                [~, init_dict] = kmeans(Y_members', l0_resolve, 'EmptyAction', 'singleton');
                D0(:, (atom_idx-1)*l0_resolve+1:(atom_idx)*l0_resolve) = ...
                init_dict';
            else
                rand_idx = randperm(size(Y_members, 2));
                D0(:, (atom_idx-1)*l0_resolve+1:(atom_idx)*l0_resolve) = ...
                Y_members(:, rand_idx(1:l0_resolve));
            end
    	end
	end
	norms = sqrt(diag(D0'*D0));
	init_dict = bsxfun(@rdivide, D0, norms');
end
