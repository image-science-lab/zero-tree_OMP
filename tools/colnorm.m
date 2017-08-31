function [normalized_mat, norms] = colnorm(mat)
    % Function to normalize columns of a matrix.
    norms = sqrt(sum(mat.*mat, 1));
    norms(norms == 0) = 1;
    N = length(norms);
    normalized_mat = mat*sparse(1:N, 1:N, 1./norms, N, N);
end