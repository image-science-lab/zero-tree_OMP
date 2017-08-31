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
	if ((nlhs != 2) || (nrhs != 9))
		mexErrMsgIdAndTxt("OMP:NInputs", "Invalid number of inputs or outputs");

	/* Matrix and array holder */
	double *D1, *D0, *Y, *X1_Pr, *X0_Pr, *U_Pr, *A_Pr;
	size_t *X1_Ir, *X1_Jc, *X0_Ir, *X0_Jc, *U_Ir, *U_Jc, *A_Ir, *A_Jc;
	int k0, k1, hf;
	double thres;
	
	/* Sanity check */
	check_matrix(prhs[0], "OMP:Dictionary");
	check_matrix(prhs[1], "OMP:Dictionary");
	check_matrix(prhs[2], "OMP:Signals");
	check_matrix(prhs[3], "OMP:Resizer");
	check_matrix(prhs[4], "OMP:Operator");
	check_scalar(prhs[5], "OMP:Sparsity");
	check_scalar(prhs[6], "OMP:Sparsity");
	check_scalar(prhs[7], "OMP:ErrorTolerance");
	check_scalar(prhs[8], "OMP:Mode");

	/*Sparse matrices check */
	check_sparse(prhs[3], "OMP:Resizer_sparse");
	check_sparse(prhs[4], "OMP:Operator_sparse");

	/* Get the sizes of Y, D1, D0, and A */
	mwSize d_m1, d_n1, d_m0, d_n0, y_m, y_n, a_m, a_n, u_m, u_n;

	d_m1 = mxGetM(prhs[0]); d_n1 = mxGetN(prhs[0]);
	d_m0 = mxGetM(prhs[1]); d_n0 = mxGetN(prhs[1]);
	y_m = mxGetM(prhs[2]); y_n = mxGetN(prhs[2]);
	u_m = mxGetM(prhs[3]); u_n = mxGetN(prhs[3]);
	a_m = mxGetM(prhs[4]); a_n = mxGetN(prhs[4]);
	
	/* Sanity check for upscale operator */
	if (u_m != d_m0 || u_n != d_m1)
		mexErrMsgIdAndTxt("OMP:Dim", "Upscaler dimension mismatch");

	/* Sanity check for CS operator */
	if (a_m != y_m || a_n != d_m0)
		mexErrMsgIdAndTxt("OMP:Dim", "CS operator dimension mismatch");

	/* Copy the matrix data */
	D1 = mxGetPr(prhs[0]);
	D0 = mxGetPr(prhs[1]);
	Y = mxGetPr(prhs[2]);

	U_Pr = mxGetPr(prhs[3]);
	U_Ir = mxGetIr(prhs[3]);
	U_Jc = mxGetJc(prhs[3]);

	A_Pr = mxGetPr(prhs[4]);
	A_Ir = mxGetIr(prhs[4]);
	A_Jc = mxGetJc(prhs[4]);

	/* Copy the scalar data */
	k1 = (int)mxGetScalar(prhs[5]);
	k0 = (int)mxGetScalar(prhs[6]);
	thres = mxGetScalar(prhs[7]);
	hf = (int)mxGetScalar(prhs[8]);
	
	/* Create output data pointer */
	plhs[0] = mxCreateSparse(d_n1, y_n, k1*y_n, mxREAL);
	plhs[1] = mxCreateSparse(d_n0, y_n, k0*y_n, mxREAL);
	
	X1_Pr = mxGetPr(plhs[0]);
	X1_Ir = mxGetIr(plhs[0]);
	X1_Jc = mxGetJc(plhs[0]);

	X0_Pr = mxGetPr(plhs[1]);
	X0_Ir = mxGetIr(plhs[1]);
	X0_Jc = mxGetJc(plhs[1]);

	/* Now invoke the function */
	c_omp_cs(Y, D1, X1_Pr, X1_Ir, X1_Jc, k1, d_m1, d_n1, D0, X0_Pr, X0_Ir,
		   	 X0_Jc, k0, d_m0, d_n0, A_Pr, A_Ir, A_Jc, a_m, a_n, y_n, U_Pr,
			 U_Ir, U_Jc, thres, hf);
	return;
}

