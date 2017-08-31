/*
 * C backend routines for some common operations
 */

#include<utils.h>

/* 2D double matrix checking function */
void check_matrix(const mxArray *matrix, const char *context)
{
	/* Double or quit */
	if (!mxIsDouble(matrix) || mxIsComplex(matrix))
		mexErrMsgIdAndTxt(context, "Need a double valued matrix");
	/* We are dealing with 2D matrices only */
	if (mxGetNumberOfDimensions(matrix) >2)
		mexErrMsgIdAndTxt(context, "Need a 2D matrix only");
}

/* 2D UINT16 matrix checking function */
void check_int_matrix(const mxArray *matrix, const char *context)
{
	/* Double or quit */
	if ((!mxIsUint8(matrix) && !mxIsUint16(matrix) && !mxIsUint32(matrix)) ||
		mxIsComplex(matrix))
		mexErrMsgIdAndTxt(context, "Need an integer valued matrix");
	/* We are dealing with 2D matrices only */
	if (mxGetNumberOfDimensions(matrix) >2)
		mexErrMsgIdAndTxt(context, "Need a 2D matrix only");
}

/* Sparse matrix checking function */
void check_sparse(const mxArray *matrix, const char *context)
{
	if(!mxIsSparse(matrix))
		mexErrMsgIdAndTxt(context, "Input should be a sparse matrix");
}
/* 1D column vector checking function */
void check_vector(const mxArray *vector, const char *context)
{
	/* Double or quit */
	if (!mxIsDouble(vector) || mxIsComplex(vector))
		mexErrMsgIdAndTxt(context, "Need a double valued matrix");
	/* Single dimensional vectors */
	if (mxGetNumberOfDimensions(vector) > 2)
		mexErrMsgIdAndTxt(context, "Need atmost two dimensional matrix");
	if (mxGetM(vector) != 1 && mxGetN(vector) != 1)
		mexErrMsgIdAndTxt(context, "Need a 1D vector input");
}

/* Scalar checking function */
void check_scalar(const mxArray *scalar, const char *context)
{
	/* Double or quit */
	if (!mxIsDouble(scalar) || mxIsComplex(scalar))
		mexErrMsgIdAndTxt(context, "Need a double valued matrix");
	/* Single dimensional */
	if (mxGetNumberOfDimensions(scalar) > 2)
		mexErrMsgIdAndTxt(context, "Need atmost two dimensional matrix");
	if (mxGetM(scalar) != 1 || mxGetN(scalar) != 1)
		mexErrMsgIdAndTxt(context, "Need a scalar input");
}

/* Matrix printing function. */
void print_matrix(double *matrix, int leading_dimension, int nrows, int ncols)
{
	int i, j;

	/* Now print the matrix */
	for(i = 0; i < nrows; i++)
	{
		for(j = 0; j < ncols; j++)
		{
			/* Row major printing. This is C default */
			if(leading_dimension == ROW_MAJOR){
				printf("%f ", (double)matrix[i*ncols + j]);
			}
			/* Column major printing. This is Matlab and Fortran default */
			else{
				printf("%f ", (double)matrix[i + j*nrows]);
			}
		}
		printf("\n");
	}
}

/* Matrix printing function. */
void print_int_matrix(int *matrix, int leading_dimension, int nrows, int ncols)
{
	int i, j;

	/* Now print the matrix */
	for(i = 0; i < nrows; i++)
	{
		for(j = 0; j < ncols; j++)
		{
			/* Row major printing. This is C default */
			if(leading_dimension == ROW_MAJOR){
				printf("%d ", (size_t)matrix[i*ncols + j]);
			}
			/* Column major printing. This is Matlab and Fortran default */
			else{
				printf("%d ", (size_t)matrix[i + j*nrows]);
			}
		}
		printf("\n");
	}
}

/* Matlab matrix printing function. */
void print_mx_matrix(const mxArray *mx_matrix)
{
	int nrows, ncols;
	int i, j;

	double *matrix;

	/* Get the data from the mxArray struct */
	nrows = mxGetM(mx_matrix);
	ncols = mxGetN(mx_matrix);
	matrix = mxGetPr(mx_matrix);

	/* Now print the matrix */
	for(i = 0; i < nrows; i++)
	{
		for(j = 0; j < ncols; j++)
		{
			/* Column major printing. This is Matlab and Fortran default */
			printf("%f \n", matrix[i + j*nrows]);
		}
		printf("\n");
	}
}

/* Non contiguous data copy utility */
void data_copy(double *dest, double *src, int nrows, int ncols, int *col_idx)
{
	int i, j;

	for(i = 0; i < nrows; i++)
		for(j = 0; j < ncols; j++)
		{
			/* Copy the col_idx[j]th column to jth column */
			dest[i + j*nrows] = src[i + col_idx[j]*nrows];
		}
}

int max(double *array, int nelem)
{
	int i;
	int max_idx = 0;
	double max_val = 0;

	for(i = 0; i < nelem; i++){
		if(max_val < fabs(array[i]))
		{
			max_val = fabs(array[i]);
			max_idx = i;
		}
	}
	return max_idx;
}
