/*
 * C routine to check the usage of BLAS functions for C.
 */

#include<mex.h>

/* Declare the blas matrix multiplication primitive */
extern void dgemv_(char *trans, int *m, int *n,
				   double *alpha, double a[], int *lda,
				   double *x, int *incx, double *beta,
				   double *y, int *incy);

void test_blas(double dictionary[], double *in_array, double *out_array, int rsize, int csize)
{
	/*Function to find the dot product between a dictionary and vector*/
	char trans = 'N'; 	/* Don't take transpose */
	double alpha = 1;
	double beta = 0;
	int incx = 1;
	int incy = 1;

	/* Call the BLAS function now */
	dgemv_(&trans, &rsize, &csize, &alpha, dictionary, &rsize, 
		   in_array, &incx, &beta, out_array, &incy);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Need to make sure that the number of arguments are correct.*/
	if ((nlhs != 1) || (nrhs != 2)){
		return;
    }
	/* Matrix and array holder */
	double *dictionary;
	double *in_array;
	double *out_array;
	
	mwSize nrows, ncols;

	nrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);

	/* Copy the data */
	dictionary = mxGetPr(prhs[0]);
	in_array = mxGetPr(prhs[1]);

	/* Create output data pointer */
	plhs[0] = mxCreateDoubleMatrix(nrows, 1, mxREAL);
	out_array = mxGetPr(plhs[0]);

	/* Now invoke the function */
	test_blas(dictionary, in_array, out_array, nrows, ncols);

	return;
}
