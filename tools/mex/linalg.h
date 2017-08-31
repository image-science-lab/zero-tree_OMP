/*
 * Header file for C backend routines for some common linear algebra opearations.
*/

#ifndef __LINALG_H__
#define __LINALG_H__

/* Includes */
#include<mex.h>
#include<math.h>
#include<string.h>

/* In case debugging is necessary */
#define DEBUG 0
#include<stdio.h>

/* Function prototypes */

/* External Fortran declarations first */

/* DGELS -- Linear least squares */
extern void dgels_(char *trans, int *m, int *n, int *nrhs, double A[],
		int *lda, double B[], int *ldb, double *work, int *lwork, int *info);

extern void dgemv_(char *trans, int *m, int *n, double *alpha, double a[],
		int *lda, double *x, int *incx, double *beta, double *y, int *incy);

extern double ddot_(int *n, double *x, int *incx, double *y, int *incy);

extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
				   double *alpha, double A[], int *lda, double B[],
				   int *ldb, double *beta, double C[], int *ldc);

/* Least squares fitting solution using Fortran DGELS */
void lstsq(double A[], double B[], double X[], int nrows, int ncols, int ndata);

/* Matrix matrix product */
void dotmm(double matrix1[], double matrix2[], double output[],
		   int nrows1, int ncols1, int nrows2, int ncols2);

/* Dot product of a matrix and a vector */
void dotm(char trans, double matrix[], double *vector, double *output,
		  int nrows, int ncols, double beta);

/* Dot product of two vectors */
double dotv(double *vector1, double *vector2, int n);

/* Matrix column normalization */
double normalize_mat(double matrix[], double *norms, int nrows, int ncols);

#endif
