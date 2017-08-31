/*
 * C routine to check the usage of BLAS functions for C.
 */

#include<mex.h>
#include<utils.h>
#include<linalg.h>
#include<sparse.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Need to make sure that the number of arguments are correct.*/
	if ((nlhs != 1) || (nrhs != 2))
		return;

	/* Matrix and array holder */
	double *B, *Pr, *X;

	size_t *Ir, *Jc, nzmax;
	
	/* Check for double matrices */
	check_matrix(prhs[0], "Matrix error");
	check_matrix(prhs[1], "Matrix error");

	mwSize nrows1, ncols1, ncols2;

	nrows1 = mxGetM(prhs[0]);
	ncols1 = mxGetN(prhs[0]);
	ncols2 = mxGetN(prhs[1]);
	

	/* Sparse matrix data. */
	Pr = mxGetPr(prhs[0]);
	Ir = mxGetIr(prhs[0]);
	Jc = mxGetJc(prhs[0]);

	nzmax = Jc[ncols1];

	B = mxGetPr(prhs[1]);
	
	/* Create output data pointer */
	plhs[0] = mxCreateDoubleMatrix(nrows1, ncols2, mxREAL);
	X = mxGetPr(plhs[0]);

	/* Now invoke the function */
	sp_dotm(Pr, Ir, Jc, B, X, nrows1, ncols1, ncols2);
	return;
}
