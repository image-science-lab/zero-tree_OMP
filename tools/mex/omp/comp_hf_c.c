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
	if ((nlhs != 2) || (nrhs != 8))
		mexErrMsgIdAndTxt("OMP:NInputs", "Invalid number of inputs or outputs");

	/* Matrix and array holder */
	double *D1, *D0, *Y, *X1_Pr, *X0_Pr, *R_Pr, *U_Pr;
	size_t *R_Ir, *R_Jc, *U_Ir, *U_Jc, *X1_Ir, *X1_Jc, *X0_Ir, *X0_Jc;
	int k0, k1;
	double thres;

	/* Sanity check */
	check_matrix(prhs[0], "OMP:Dictionary");
	check_matrix(prhs[1], "OMP:Dictionary");
	check_matrix(prhs[2], "OMP:Signals");
	check_matrix(prhs[3], "OMP:Resizer");
	check_sparse(prhs[3], "OMP:Resizer");
	check_matrix(prhs[4], "OMP:Resizer");
	check_sparse(prhs[4], "OMP:Resizer");
	check_scalar(prhs[5], "OMP:Sparsity");
	check_scalar(prhs[6], "OMP:Sparsity");
	check_scalar(prhs[7], "OMP:ErrorTolerance");

	/* Get the sizes of Y and D */
	mwSize d_m1, d_n1, d_m0, d_n0, y_m, y_n;

	d_m1 = mxGetM(prhs[0]); d_n1 = mxGetN(prhs[0]);
	d_m0 = mxGetM(prhs[1]); d_n0 = mxGetN(prhs[1]);
	y_m = mxGetM(prhs[2]); y_n = mxGetN(prhs[2]);

	/* Need number of rows in D and Y to be same */
	if(d_m0 != y_m)
		mexErrMsgIdAndTxt("OMP:NRows", "Number of rows in D and Y have to be same");

	/* Copy the matrix data */
	D1 = mxGetPr(prhs[0]);
	D0 = mxGetPr(prhs[1]);
	Y = mxGetPr(prhs[2]);

	/* R has to be sparse */
	R_Pr = mxGetPr(prhs[3]);
	R_Ir = mxGetIr(prhs[3]);
	R_Jc = mxGetJc(prhs[3]);

	/* Upsampling operator */
	U_Pr = mxGetPr(prhs[4]);
	U_Ir = mxGetIr(prhs[4]);
	U_Jc = mxGetJc(prhs[4]);

	/* Copy the scalar data */
	k1 = (int)mxGetScalar(prhs[5]);
	k0 = (int)mxGetScalar(prhs[6]);
	thres = mxGetScalar(prhs[7]);

	/* Create sparse output data pointer */
	plhs[0] = mxCreateSparse(d_n1, y_n, k1*y_n, mxREAL);
	plhs[1] = mxCreateSparse(d_n0, y_n, k0*y_n, mxREAL);

	X1_Pr = mxGetPr(plhs[0]);
	X1_Ir = mxGetIr(plhs[0]);
	X1_Jc = mxGetJc(plhs[0]);

	X0_Pr = mxGetPr(plhs[1]);
	X0_Ir = mxGetIr(plhs[1]);
	X0_Jc = mxGetJc(plhs[1]);

	/* Now invoke the function */
	c_omp_hf(Y, D1, X1_Pr, X1_Ir, X1_Jc, k1, d_m1, d_n1, D0, X0_Pr, X0_Ir, X0_Jc,
		  k0, d_m0, d_n0, y_n, R_Pr, R_Ir, R_Jc, U_Pr, U_Ir, U_Jc, thres, 1);

	return;
}
