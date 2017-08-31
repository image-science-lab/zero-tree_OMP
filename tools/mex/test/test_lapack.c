/*
 * C routine to check the usage of BLAS functions for C.
 */

#include<mex.h>
#include<stdio.h>
#include<string.h>

/* Declare the lapack linear least square primitive */
extern void dgels_(char *trans, int *m, int *n, int *nrhs, 
				   double A[], int *lda, double *B, int *ldb,
				   double *work, int *lwork, int *info);

void test_lapack(double A[], double B[], double *X, int rsize, int csize, int ndata)
{
	/*Function to find the linear least square solutions min.||B - AX||^2*/
	char trans = 'N'; 	/* Don't take transpose */
	int nrhs = ndata;
	int lda = rsize, ldb = rsize;
	int lwork = -1;
	double wkopt;
	double *work;
	int info;

	/* Create ghosts to store A, B */
	double *ghost_A, *ghost_B;
	ghost_A = (double *)mxMalloc(rsize*csize*sizeof(double));
	ghost_B = (double *)mxMalloc(rsize*ndata*sizeof(double));

	/* Copy the data */
	memcpy(ghost_A, A, rsize*csize*sizeof(double));
	memcpy(ghost_B, B, rsize*ndata*sizeof(double));

	/* Need to make a query for work size */
	dgels_(&trans, &rsize, &csize, &nrhs, ghost_A, &lda, ghost_B, &ldb, &wkopt, &lwork, &info);
	/* work now has the optimum work size */
	lwork = (int)wkopt;
	work = (double *)mxMalloc(lwork*sizeof(double));
	/* Call the routine again, this time to solve the solution */
	dgels_(&trans, &rsize, &csize, &nrhs, ghost_A, &lda, ghost_B, &ldb, work, &lwork, &info);
	/* Free the work area */
	mxFree(work);

	/* Turns out, B is the required solution! */
	printf("nrows = %d; ncols = %d; ndata = %d\n", rsize, csize, ndata);	
	int i = 0, j = 0;
	for(i = 0; i < csize; i++){
		for (j = 0; j < ndata; j++){
			X[i + j*csize] = ghost_B[i + j*rsize];
			printf("%d -- %f ",i + j*rsize, X[i + j*csize]);
		}
		printf("\n");
	}
	printf("info = %d\n", info);
	/*memcpy(X, B, ndata*csize*sizeof(double));*/

	/* Free the ghost arrays */
	mxFree(ghost_A);
	mxFree(ghost_B);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Need to make sure that the number of arguments are correct.*/
	if ((nlhs != 1) || (nrhs != 2))
		return;

	/* Matrix and array holder */
	double *A;
	double *B;
	double *X;
	
	mwSize nrows, ncols, ndata;

	nrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	ndata = mxGetN(prhs[1]);
	
	printf("%d\n", ndata);

	/* Copy the data */
	A = mxGetPr(prhs[0]);
	B = mxGetPr(prhs[1]);
	
	/* Create output data pointer */
	plhs[0] = mxCreateDoubleMatrix(ncols, ndata, mxREAL);
	X = mxGetPr(plhs[0]);

	/* Now invoke the function */
	test_lapack(A, B, X, nrows, ncols, ndata);

	return;
}
