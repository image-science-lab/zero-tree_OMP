/*
 * C backend routines for some common linear algebra operations
 */

#include<linalg.h>

/* Linear least squares fit */
void lstsq(double A[], double B[], double X[], int nrows, int ncols, int ndata)
{
	/* Some constants for passing to the Fortran subroutine. */
	char trans = 'N'; 			/* Don't take transpose */
	int nrhs = ndata;  			/* NRHS is the number of data points */
	int lda = nrows;
	int ldb = nrows; 			/* Leading array dimension of A and B */
	int lwork = -1; 			/* Used for quering the optimal work size */
	double wkopt; 				/* Dummy variable for work query */
	double *work; 				/* Work array */
	int info; 					/* Information of success */

	/* The fortran subroutine modifies the A and B array. Hence, to avoid
	 * the problem, we create copies of A and B. */
	double *copy_A, *copy_B;

	copy_A = (double *)mxMalloc(nrows*ncols*sizeof(double));
	copy_B = (double *)mxMalloc(nrows*ndata*sizeof(double));

	memcpy(copy_A, A, nrows*ncols*sizeof(double));
	memcpy(copy_B, B, nrows*ndata*sizeof(double));

	/* Make a query for optimal work size */
	dgels_(&trans, &nrows, &ncols, &nrhs, copy_A, &lda, copy_B, &ldb,
			&wkopt, &lwork, &info);

	/* Got the optimal work size. Now allocate the memory */
	lwork = (int)wkopt;
	work = (double *)mxMalloc(lwork*sizeof(double));

	/* Call the routine to solve the problem this time.*/
	dgels_(&trans, &nrows, &ncols, &nrhs, copy_A, &lda, copy_B, &ldb,
			work, &lwork, &info);

	/* Solution is done. Since we have an over-determined problem, the
	 * solution is contained in the first ncols rows of copy_B. We will
	 * manually copy this data to X */
	int i, j = 0;
	for(i = 0; i < ncols; i++)
		for(j = 0; j < ndata; j++)
			/* Matlab stores data in column-major fasion. */
			X[i + j*ncols] = copy_B[i + j*nrows];

	/* All done. Free the temporary data */
	mxFree(work);
	mxFree(copy_A);
	mxFree(copy_B);
}

/* Matrix and matrix dot product */
/* TODO: Make it generic for transposes also */
void dotmm(double matrix1[], double matrix2[], double output[],
		  int nrows1, int ncols1, int nrows2, int ncols2)
{
	/* Data to send to the fortran subroutine */
	double alpha = 1;
	double beta = 0;
	char transa = 'n';
	char transb = 'n';

	/* Now call the BLAS subroutine */
	dgemm_(&transa, &transb, &nrows1, &ncols2, &ncols1, &alpha,
		   matrix1, &nrows1, matrix2, &nrows2, &beta, output, &nrows1);
}

/* Matrix and vector dot product */
void dotm(char trans, double matrix[], double *vector, double *output, int nrows, int ncols, double beta)
{
	/* Data to send to the fortran subroutine. */
	double alpha = 1; 			/* No scalar multiplication */
	int incx = 1, incy = 1; 	/* Increments. */

	/* Now call the BLAS subroutine. */
	dgemv_(&trans, &nrows, &ncols, &alpha, matrix, &nrows, vector,
			&incx, &beta, output, &incy);
}

/* Vector dot product */
double dotv(double *vector1, double*vector2, int n)
{
	/* This function is at best a wrapper to the actual dot function */
	int incx = 1, incy = 1;
	double vec_product = 0;

	/* Call the subroutine */
	vec_product = ddot_(&n, vector1, &incx, vector2, &incy);
	return vec_product;
}

/* Normalization of matrix columns */
/* TODO: Add a switch for normalizing rows as well. */
double normalize_mat(double matrix[], double *norms, int nrows, int ncols)
{
	int row_idx, col_idx;
	double norm, temp;

	for(col_idx=0; col_idx < ncols; col_idx++)
	{
		/* First calculate the norm */
		norm = 0;

		for(row_idx = 0; row_idx < nrows; row_idx++)
		{
			temp = matrix[row_idx + col_idx*nrows];
			norm += temp*temp;
		}
		norm = sqrt(norm);
		norms[col_idx] = norm;

		/* Second pass, normalize the column */
		for(row_idx = 0; row_idx < nrows; row_idx++)
		{
			temp = matrix[row_idx + col_idx*nrows];
			matrix[row_idx + col_idx*nrows] = temp/norm;
		}
	}
}
