function X0 = solveXConstrained(Y, D0, X1, k0, l0_resolve)
    % Function to calculate the sparse approximation of a set of signals
	% The further constrain is given by the lower resolution sparse
	% signal X1.
	% Input:
	%		Y: Set of signal vectors which need to be approximated.
	%		D0: High resolution dictionary. We essentially are solving
	%			Y = D0*X0 but with some further constraints on the
	%			support set of X0.
	%		X1: Low resolution sparse approximation. The non zero
	%			entries will decide the sub dictionary to use.
	%		k0: Sparsity of X0.
	%		l0_resolve: Number of high resolution atoms per low
	%			resolution atom.

    % We want to store X as a sparse matrix.
    row_idx = zeros(k0, size(Y,2));
    col_idx = zeros(k0, size(Y,2));
    values = zeros(k0, size(Y, 2));

    parfor idx = 1:size(Y,2)
        if norm(Y(:, idx)) ~= 0
            [x, support_set] = sparse_comp(D0, Y(:,idx), X1(:, idx), k0, l0_resolve);
        else
            x = zeros(k0, 1);
            support_set = -ones(k0, 1);
        end
		row_idx(:, idx) = support_set';
        col_idx(:, idx) = idx*ones(k0, 1);
        values(:, idx) = x;
    end
    t = find(row_idx ~= -1);
    row_idx = row_idx(t);
    col_idx = col_idx(t);
    values = values(t);
    X0 = sparse(row_idx(:), col_idx(:), values(:), size(D0,2), size(Y, 2));
end
