/*
 * C routine to check the usage of BLAS functions for C.
 */

#include<mex.h>
#include<utils.h>
#include<linalg.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Need to make sure that the number of arguments are correct.*/
	if ((nlhs != 0) || (nrhs != 1))
		return;

	/* Matrix and array holder */
	double *A;
	
	/* Check for double matrices */
	check_matrix(prhs[0], "Matrix error");

	mwSize nrows1, ncols1;

	nrows1 = mxGetM(prhs[0]);
	ncols1 = mxGetN(prhs[0]);

	/* Copy the data */
	A = mxGetPr(prhs[0]);
	
	/* Now invoke the function */
	normalize_mat(A, nrows1, ncols1);

	return;
}
