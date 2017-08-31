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
	if ((nlhs != 1) || (nrhs != 4))
		mexErrMsgIdAndTxt("OMP:NInputs", "Invalid number of inputs or outputs");

	/* Matrix and array holder */
	double *D, *Y, *X_Pr;
	size_t *X_Ir, *X_Jc;
	int blksize;
	double thres;
	
	/* Sanity check */
	check_matrix(prhs[0], "OMP:Dictionary");
	check_matrix(prhs[1], "OMP:Signals");
	check_scalar(prhs[2], "OMP:BlockSize");
	check_matrix(prhs[3], "OMP:ErrorTolerance");

	/* Get the sizes of Y and D */
	mwSize d_m, d_n, y_m, y_n;

	d_m = mxGetM(prhs[0]); d_n = mxGetN(prhs[0]);
	y_m = mxGetM(prhs[1]); y_n = mxGetN(prhs[1]);
	
	/* Need number of rows in D and Y to be same */
	if(d_m != y_m)
		mexErrMsgIdAndTxt("OMP:NRows", "Number of rows in D and Y have to be same");

	/* Copy the matrix data */
	D = mxGetPr(prhs[0]);
	Y = mxGetPr(prhs[1]);

	/* Copy the scalar data */
	blksize = (int)mxGetScalar(prhs[2]);
	thres = mxGetScalar(prhs[3]);
	
	/* Create output sparse data pointer */
	plhs[0] = mxCreateSparse(d_n, y_n, (d_n/blksize)*y_n, mxREAL);
	X_Pr = mxGetPr(plhs[0]);
	X_Ir = mxGetIr(plhs[0]);
	X_Jc = mxGetJc(plhs[0]);
	
	/* Now invoke the function */
	omp_mld(Y, D, X_Pr, X_Ir, X_Jc, thres, d_m, d_n, y_n, blksize);

	return;
}

