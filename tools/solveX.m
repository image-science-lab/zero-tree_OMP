function X = solveX(Y, D, k)
    % Data parallel implementation of sparse coding using OMP.
    %
    % Inputs:
    %   Y: Data
    %   D: Dictionary
    %   k: Sparsity
    %
    % Outputs:
    %   X: Sparse coefficients.

    T = size(Y, 2);
    N = size(D, 2);
    row_idx = zeros(k, T);
    col_idx = zeros(k, T);
    vals = zeros(k, T);

    for idx=1:T
        [x, sup] = sparse_omp(D, Y(:, idx), k, 1e-5);
        row_idx(:, idx) = sup;
        col_idx(:, idx) = ones(k,1)*idx;
        vals(:, idx) = x;
    end

    row_idx = row_idx(:);
    col_idx = col_idx(:);
    vals = vals(:);
    nnz_idx = find(row_idx ~= -1);

    X = sparse(row_idx(nnz_idx), col_idx(nnz_idx), vals(nnz_idx), N, T);
end
