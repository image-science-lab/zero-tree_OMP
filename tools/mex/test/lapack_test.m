% Testing least squares routines in LAPACK.

% Compile the code.
if cc == true
	mex -v test_lapack.c -llapack -lblas -lgfortran -L../../lib
end

A = [1 2; 3 1; 9 4; 7 1];
X = [1 2 9; 7 3 5;];

y = A*X;

x_matlab = A\y;

x_mex = test_lapack(A, y);
