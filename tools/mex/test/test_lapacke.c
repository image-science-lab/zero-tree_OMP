/*
 * C routine to check the usage of BLAS functions for C.
 */

#include<mex.h>
#include<stdio.h>
#include<lapacke.h>

void test_lapack(double A[], double *B, double *X, int rsize, int csize)
{
	/*Function to find the linear least square solutions min.||B - AX||^2*/
	char trans = 'N'; 	/* Don't take transpose */
	int nrhs = 1;
	int lda = rsize, ldb = rsize;
	int lwork = -1;
	double wkopt;
	double *work;
	int info; int i;

	/* Call LAPACK_dgels from lapacke library */
	info = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', rsize, csize, nrhs, A, lda, B, ldb);
	printf("info = %d\n", info);
	for(i = 0; i < csize; i++){
		X[i] = B[i];
		printf("%f %f\n", B[i], B[i*csize]);
	}
	/*memcpy(X, B, csize*sizeof(double));*/
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Need to make sure that the number of arguments are correct.*/
	if ((nlhs != 1) || (nrhs != 2))
		return;

	/* Matrix and array holder */
	double *A;
	double *B;
	double *X;
	
	mwSize nrows, ncols;

	nrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);

	/* Copy the data */
	A = mxGetPr(prhs[0]);
	B = mxGetPr(prhs[1]);

	/* Create output data pointer */
	plhs[0] = mxCreateDoubleMatrix(ncols, 1, mxREAL);
	X = mxGetPr(plhs[0]);

	/* Now invoke the function */
	test_lapack(A, B, X, nrows, ncols);

	return;
}
