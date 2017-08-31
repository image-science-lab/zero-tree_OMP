function [X1, X0] = comp_cs(D1, D0, Y, U, A, k1, k0, thres, hf)
	% Function to find sparse coefficients with compressed measurements
	% using multi-scale dictionary.
	%
	% Inputs:
	% 	D1: Coarse level dictionary.
	% 	D0: Fine level dictionary.
	% 	Y: Compressed data.
	% 	U: Upscale operator. Should be sparse.
	% 	A: Measurement operator. Should be sparse.
	% 	k1: Coarse sparsity.
	% 	k0: Fine sparsity.
	% 	thres: OMP threshold.
	% 	hf: If 1, use high frequency model.
	%
	% Output:
	% 	X1: Coarse level sparse coefficients.
	% 	X0: Fine level sparse coeffcients.

	% Modify dictionaries.
	D1mod = A*U*D1;
    [D1mod, norms1] = colnorm(D1mod);

	D0mod = A*D0;
    [D0mod, norms0] = colnorm(D0mod);

	% Solve coarse approximation problem.
	X1 = solveX(Y, D1mod, k1);
	
	% Solve the fine approximation problem.
	if hf
		Y = Y - D1mod*X1;
	end
	X0 = solveXConstrained(Y, D0mod, X1, k0, size(D0, 2)/size(D1, 2));

	% Renormalize the sparse coefficients.
	X1 = X1./(norms1'*ones(1, size(Y, 2)));
	X0 = X0./(norms0'*ones(1, size(Y, 2)));
end
