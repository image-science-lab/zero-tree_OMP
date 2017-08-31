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
	double *D, *Y, *X_Pr, *A_Pr;
	size_t *X_Ir, *X_Jc, *A_Ir, *A_Jc;
	int blksize;
	double thres;
	
	/* Sanity check */
	check_matrix(prhs[0], "OMP:Dictionary");
	check_matrix(prhs[1], "OMP:Signals");
	check_matrix(prhs[2], "OMP:Operator");
	check_scalar(prhs[3], "OMP:BlockSize");
	check_scalar(prhs[4], "OMP:ErrorTolerance");

	/* Get the sizes of Y and D */
	mwSize d_m, d_n, y_m, y_n, a_m, a_n;

	d_m = mxGetM(prhs[0]); d_n = mxGetN(prhs[0]);
	y_m = mxGetM(prhs[1]); y_n = mxGetN(prhs[1]);
	a_m = mxGetM(prhs[2]); a_n = mxGetN(prhs[2]);
	
	/* Need number of rows in D and Columns in A need to be same */
	if(d_m != a_n)
		mexErrMsgIdAndTxt("OMP:NRows", "Dimension mismatch between A and D");
	if(a_m != y_m)
		mexErrMsgIdAndTxt("OMP:NRows", "Dimension mismatch for data");

	/* Copy the matrix data */
	D = mxGetPr(prhs[0]);
	Y = mxGetPr(prhs[1]);

	/* Get the sparse operator */
	A_Pr = mxGetPr(prhs[2]);
	A_Ir = mxGetIr(prhs[2]);
	A_Jc = mxGetJc(prhs[2]);

	/* Copy the scalar data */
	blksize = (int)mxGetScalar(prhs[3]);
	thres = mxGetScalar(prhs[4]);
	
	/* Create output sparse data pointer */
	plhs[0] = mxCreateSparse(d_n, y_n, (d_n/blksize)*y_n, mxREAL);
	X_Pr = mxGetPr(plhs[0]);
	X_Ir = mxGetIr(plhs[0]);
	X_Jc = mxGetJc(plhs[0]);
	
	/* Now invoke the function */
	mld_omp_cs(Y, D, X_Pr, X_Ir, X_Jc, d_m, d_n, A_Pr, A_Ir, A_Jc, a_m, a_n,
		       y_n, blksize, thres);
	return;
}

