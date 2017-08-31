/*
 * C implementation of Orthogonal Matching pursuit using BLAS and
 * LAPACK routines for ultra fast operation.
 */

/* Includes */
#include <mex.h>
#include <utils.h>
#include <linalg.h>
#include <omp_core.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Need to make sure that the number of arguments are correct.*/
	if ((nlhs != 1) || (nrhs != 5))
		mexErrMsgIdAndTxt("OMP:NInputs", "Invalid number of inputs or outputs");

	/* Matrix and array holder */
	double *D1, *D0, *Y, *X1, *X0;
	int k0, l0_resolve;
	double thres;
	
	/* (D0, y, x1, k0, l0_resolve) */
	/* Sanity check */
	check_matrix(prhs[0], "OMP:Dictionary");
	check_matrix(prhs[1], "OMP:Signals");
	check_matrix(prhs[2], "OMP:Signals");
	check_scalar(prhs[3], "OMP:Sparsity");
	check_scalar(prhs[4], "OMP:L0Resolve");

	/* Get the sizes of Y and D0, X1 */
	mwSize x1_m1, x1_n1, d_m0, d_n0, y_m, y_n;

	x1_m1 = mxGetM(prhs[2]); x1_n1 = mxGetN(prhs[2]);
	d_m0 = mxGetM(prhs[0]); d_n0 = mxGetN(prhs[0]);
	y_m = mxGetM(prhs[1]); y_n = mxGetN(prhs[1]);
	
	/* Need number of rows in D and Y to be same */
	if(d_m0 != y_m)
		mexErrMsgIdAndTxt("OMP:NRows", "Number of rows in D and Y have to be same");

	/* Copy the matrix data */
	D0 = mxGetPr(prhs[0]);
	Y = mxGetPr(prhs[1]);
	X1 = mxGetPr(prhs[2]);
	
	/* Copy the scalar data */
	k0 = (int)mxGetScalar(prhs[3]);
	l0_resolve = mxGetScalar(prhs[4]);
	
	/* Create output data pointer */
	plhs[0] = mxCreateDoubleMatrix(d_n0, y_n, mxREAL);
	X0 = mxGetPr(plhs[0]);

	thres = 1e-4;	
	/* Now invoke the function */
	c_omp_std(Y, X1, l0_resolve, D0, X0, k0, d_m0, d_n0, y_n, x1_n1);
	
	return;
}

