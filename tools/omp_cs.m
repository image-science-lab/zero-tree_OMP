function X = omp_cs(D, Y, A, k, thres)
	% Function to extract sparse coefficients with compressed
	% measurements.
	%
	% 	Inputs
	% 		D: Dictionary.
	% 		Y: Data.
	% 		A: Operator. Should be sparse.
	% 		k: Sparsity.
	% 		thres: OMP threshold.
	%
	% 	Output:
	% 		X: Sparse coefficients for data.

	% Modify the dictionary
	Dmod = A*D;
    [Dmod, norms] = colnorm(Dmod);

	% Solve for X now.
	X = solveX(Y, Dmod, k);

	% Renormalize X
    parfor idx=1:size(X, 2)
        X(:, idx) = X(:, idx)./norms';
    end
end
