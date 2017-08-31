/*
 * C implementation of Orthogonal Matching pursuit using BLAS and
 * LAPACK routines for ultra fast operation.
 */

/* Includes */
#include <mex.h>
#include <utils.h>
#include <linalg.h>
#include <sparse.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Need to make sure that the number of arguments are correct.*/
	if ((nlhs != 1) || (nrhs != 3))
		mexErrMsgIdAndTxt("OMP:NInputs", "Invalid number of inputs or outputs");

	/* Matrix and array holder */
	double *D, *x, *y;
	int *sup;
	int i;
	
	/* Sanity check */
	check_matrix(prhs[0], "MUL:Dictionary");
	check_matrix(prhs[1], "MUL:Coefficients");
	check_int_matrix(prhs[2], "MUL:Support");

	/* Get the sizes of Y and D */
	mwSize d_m, d_n, x_m, x_n, sup_m, sup_n;

	d_m = mxGetM(prhs[0]); d_n = mxGetN(prhs[0]);
	x_m = mxGetM(prhs[1]); x_n = mxGetN(prhs[1]);
	sup_m = mxGetM(prhs[2]); sup_n = mxGetN(prhs[2]);
	
	/* Need number of rows in D and Y to be same */
	if(x_m != sup_m)
		mexErrMsgIdAndTxt("MUL:NVals", "Number of elements in x and sup have to be same");

	/* Copy the matrix data */
	D = mxGetPr(prhs[0]);
	x = mxGetPr(prhs[1]);
	sup = (int *)mxGetPr(prhs[2]);

	/* Create output data pointer */
	plhs[0] = mxCreateDoubleMatrix(d_m, x_n, mxREAL);
	y = mxGetPr(plhs[0]);

	/* Matlab indexing and C indexing are off by 1 */
	for(i = 0; i < x_m; i++)
		sup[i] -= 1;

	/* Now invoke the function */
	Dxmul(y, x, sup, D, (int)d_m, (int)x_m);

	return;
}

