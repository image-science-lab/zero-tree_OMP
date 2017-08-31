/*
 * C implementation of Orthogonal Matching pursuit using BLAS and
 * LAPACK routines for ultra fast operation.
 */

/* Includes */
#include <mex.h>
#include <matrix.h>
#include <utils.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Need to make sure that the number of arguments are correct.*/
	if ((nlhs != 1) || (nrhs != 3))
		mexErrMsgIdAndTxt("OMP:NInputs", "Invalid number of inputs or outputs");

	/* Matrix and array holder */
	double *D, *y, *proxy;
	int *support;
	
	/* Sanity check */
	check_matrix(prhs[0], "SDOT:Dictionary");
	check_vector(prhs[1], "SDOT:Signals");
	check_int_matrix(prhs[2], "SDOT:Indices");

	/* Get the sizes of Y and D and support*/
	mwSize d_m, d_n, y_m, y_n, sup_m, sup_n, sup_size;

	d_m = mxGetM(prhs[0]); d_n = mxGetN(prhs[0]);
	y_m = mxGetM(prhs[1]); y_n = mxGetN(prhs[1]);
	sup_m = mxGetM(prhs[2]); sup_n = mxGetN(prhs[2]);
	
	/* Get vector size of support */
	sup_size = (sup_m > sup_n)?sup_m:sup_n;

	/* Need number of rows in D and Y to be same */
	if(d_m != y_m)
		mexErrMsgIdAndTxt("SDOT:NRows", "Number of rows in D and Y have to be same");

	if (y_n != 1)
		mexErrMsgIdAndTxt("SDOT:COLVECTOR", "Need a column vector");

	/* Copy the matrix data */
	D = mxGetPr(prhs[0]);
	y = mxGetPr(prhs[1]);

	/* Copy the scalar data */
	support = (int *)mxGetPr(prhs[2]);
		
	/* Create output sparse data pointer */
	plhs[0] = mxCreateDoubleMatrix(sup_size, 1, mxREAL);
	proxy = mxGetPr(plhs[0]);

	/* Helper variables */
	int i, j;

	for (i = 0; i < sup_size; i++)
	{
		proxy[i] = 0;
		for(j = 0; j < d_m; j++)
		{
			/* There is a 1 shift between matlab indices and C indices */
			proxy[i] += D[(support[i]-1)*d_m + j]*y[j];
		}
	}
	
	return;
}

