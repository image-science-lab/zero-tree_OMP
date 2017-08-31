/*
 * C routine to test the usage of C backend for MATLAB functions.
 */

#include<mex.h>
#include<stdio.h>
void test(double in_array[], double out_array[], int rsize, int csize)
{
	/* Return a matrix formed by adding 1 to thie current matrix. */
	int i, j;
	printf("nrows -- %d, ncols -- %d\n", rsize, csize);
	for (i = 0; i < rsize; i++)
		for (j = 0; j < csize; j++){
			printf("%d.", i*rsize + j);
			out_array[i*csize + j] = in_array[i*csize + j] + 1;
		}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* I think I should write a routine called check_arguments(), which
	 will check for the sanity of the arguments. But for now, pack. */
	
	double *in_matrix; 		/* Input matrix. */
	double *out_matrix; 	/* Output matrix. */

	/* Create size holder. */
	mwSize nrows, ncols;

	nrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);

	/* Copy data. */
	in_matrix = mxGetPr(prhs[0]);

	/* create data pointer for output. */
	plhs[0] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
	out_matrix = mxGetPr(plhs[0]);

	/* Now invoke our function. */
	test(in_matrix, out_matrix, nrows, ncols);

	return;
}
