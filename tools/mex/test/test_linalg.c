/*
 * C routine to check the usage of BLAS functions for C.
 */

#include<mex.h>
#include<utils.h>
#include<linalg.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Need to make sure that the number of arguments are correct.*/
	if ((nlhs != 1) || (nrhs != 2))
		return;

	/* Matrix and array holder */
	double *A;
	double *B;
	double *X;
	
	/* Check for double matrices */
	check_matrix(prhs[0]);
	check_matrix(prhs[1]);

	mwSize nrows, ncols, ndata;

	nrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	ndata = mxGetN(prhs[1]);
	
	printf("%d\n", ndata);

	/* Copy the data */
	A = mxGetPr(prhs[0]);
	B = mxGetPr(prhs[1]);
	
	/* Create output data pointer */
	plhs[0] = mxCreateDoubleMatrix(ncols, ndata, mxREAL);
	X = mxGetPr(plhs[0]);

	/* Now invoke the function */
	lstsq(A, B, X, nrows, ncols, ndata);

	return;
}
