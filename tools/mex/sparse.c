/*
 * Sparse matrix multiplication routines. Uses Compressed Sparse Column
 * format for all operations.
 */

/* Understanding CSC format:
 * The CSC format stores data as pr, ir and jc, where pr is an nzmax
 * array of data, ir is an nzmax array of row indices and jc is an ncols+1
 * array of column indices.
 *
 * If the jth column has non-zero rows, then jc[j] to jc[j+1] is the indices
 * in ir where the rows have non-zero values.
 */

#include<sparse.h>

/* Sparse matrix, vector multiplication routine */
void sp_dotv(double *pr, size_t *ir, size_t *jc, double *y, double *output,
	   	     int nrows, int ncols)
{
	int i, j;

	/* Make sure that output is initialized to zero */
	for(i = 0; i < nrows; i++)
		output[i] = 0;

	/* Now run the matrix multiplication. */
	for(i = 0; i < ncols; i++)
	{
		for(j = jc[i]; j < jc[i+1]; j++)
			output[ir[j]] += pr[j]*y[i];
	}
}

/* Sparse matrix and dense matrix multiplication routine */
void sp_dotm(double *pr, size_t *ir, size_t *jc, double *mat, double *output,
	     	 int nrows1, int ncols1, int ncols2)
{
	int i, j, k;

	/* Make sure output is zero */
	for(i = 0; i < ncols2; i++)
		for(j = 0; j < nrows1; j++)
			output[j + i*nrows1] = 0;

	/* Now run the matrix multiplication */
	for(i = 0; i < ncols2; i++)
		for(j = 0; j < ncols1; j++)
			for(k = jc[j]; k < jc[j+1]; k++)
				output[i*nrows1 + ir[k]] += pr[k]*mat[i*ncols1 + j];
}

/* Dictionary and sparse vector multiplication */
void Dxmul(double *y, double *x, int *sup, double *D, int nrows, int k)
{
	int i, j;
	for(i = 0; i < k; i++)
	{
		if(sup[i] == -1)
			return;
		for(j = 0; j < nrows; j++)
		{
			y[j] += D[nrows*sup[i] + j]*x[i];
		}
	}
}
