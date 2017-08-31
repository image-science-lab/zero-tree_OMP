/*
 * C routine to check the usage of BLAS functions for C.
 */

#include<mex.h>
#include<utils.h>
#include<linalg.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Need to make sure that the number of arguments are correct.*/
	if ((nlhs != 1) || (nrhs != 2)){
		return;
    }
	/* Matrix and array holder */
	double *dictionary;
	int *in_array;
	double *out_array;
	double *temp_in_array; /* Temporary array to store the double values*/

	mwSize nrows, ncols, nidx;

	nrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	nidx = mxGetM(prhs[1])==1?mxGetN(prhs[1]):mxGetM(prhs[1]);

	/* Copy the data */
	dictionary = mxGetPr(prhs[0]);
	temp_in_array = mxGetPr(prhs[1]);
	
	int i = 0;
	in_array = (int *)malloc(nidx*sizeof(int));
	for(i = 0; i < nidx; i++)
		in_array[i] = (int)temp_in_array[i];

	/* Create output data pointer */
	plhs[0] = mxCreateDoubleMatrix(nrows, nidx, mxREAL);
	out_array = mxGetPr(plhs[0]);

	/* Now invoke the function */
	data_copy(out_array, dictionary, nrows, nidx, in_array);

	return;
}
