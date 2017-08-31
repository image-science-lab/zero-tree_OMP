/*
 * Header file for C backend routines for some sparse matrix functions
 */

#ifndef __SPARSE_H__
#define __SPARSE_H__

/* Includes */
#include<mex.h>
#include<math.h>
#include<string.h>
#include<stdio.h>

/* In case debugging is necessary */
#define DEBUG 0

/* Sparse matrix, vector multiplication routine */
void sp_dotv(double *pr, size_t *ir, size_t *jc, double *y, double *output,
	   	     int nrows, int ncols);

/* Sparse matrix and dense matrix multiplication routine */
void sp_dotm(double *pr, size_t *ir, size_t *jc, double *mat, double *output,
	     	 int nrows1, int ncols1, int ncols2);

/* Dictionary and sparse vector multiplication */
void Dxmul(double *y, double *x, int *sup, double *D, int nrows, int k);

#endif
