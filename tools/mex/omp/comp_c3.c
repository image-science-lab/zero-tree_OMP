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
	if ((nlhs != 3) || (nrhs != 10))
		mexErrMsgIdAndTxt("OMP:NInputs", "Invalid number of inputs or outputs");

	/* Matrix and array holder */
	double *D2, *D1, *D0, *Y, *X2, *X1, *X0, *R1, *R0;
	int k0, k1, k2;
	double thres;
	
	/* Sanity check */
	check_matrix(prhs[0], "OMP:Dictionary");	/*D2*/
	check_matrix(prhs[1], "OMP:Dictionary");	/*D1*/
	check_matrix(prhs[2], "OMP:Dictionary");	/*D0*/
	check_matrix(prhs[3], "OMP:Signals");		/*Y*/
	check_matrix(prhs[4], "OMP:Resizer");		/*R1*/
	check_matrix(prhs[5], "OMP:Resizer");		/*R0*/
	check_scalar(prhs[6], "OMP:Sparsity");		/*k2*/
	check_scalar(prhs[7], "OMP:Sparsity");		/*k1*/
	check_matrix(prhs[8], "OMP:Sparsity");		/*k0*/
	check_scalar(prhs[9], "OMP:ErrorTolerance");

	/* Get the sizes of Y and D */
	mwSize d_m2, d_n2, d_m1, d_n1, d_m0, d_n0, y_m, y_n;

	d_m2 = mxGetM(prhs[0]); d_n2 = mxGetN(prhs[0]);
	d_m1 = mxGetM(prhs[1]); d_n1 = mxGetN(prhs[1]);
	d_m0 = mxGetM(prhs[2]); d_n0 = mxGetN(prhs[2]);
	y_m = mxGetM(prhs[3]); y_n = mxGetN(prhs[3]);
	
	/* Need number of rows in D and Y to be same */
	if(d_m0 != y_m)
		mexErrMsgIdAndTxt("OMP:NRows", "Number of rows in D and Y have to be same");

	/* Copy the matrix data */
	D2 = mxGetPr(prhs[0]);
	D1 = mxGetPr(prhs[1]);
	D0 = mxGetPr(prhs[2]);
	Y = mxGetPr(prhs[3]);
	R1 = mxGetPr(prhs[4]);
	R0 = mxGetPr(prhs[5]);

	/* Copy the scalar data */
	k2 = (int)mxGetScalar(prhs[6]);
	k1 = (int)mxGetScalar(prhs[7]);
	k0 = (int)mxGetScalar(prhs[8]);
	thres = mxGetScalar(prhs[9]);
	
	/* Create output data pointer */
	plhs[0] = mxCreateDoubleMatrix(d_n2, y_n, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(d_n1, y_n, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(d_n0, y_n, mxREAL);
	X2 = mxGetPr(plhs[0]);
	X1 = mxGetPr(plhs[1]);
	X0 = mxGetPr(plhs[2]);
	
	/* Now invoke the function */
	c_omp3(Y, D2, X2, k2, d_m2, d_n2, D1, X1, k1, d_m1, d_n1, D0, X0, k0, d_m0, d_n0, y_n, R1, R0, thres);

	return;
}

