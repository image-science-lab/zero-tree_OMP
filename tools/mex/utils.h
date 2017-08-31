/*
 * Header file for C backend routines for some common operations.
 */

#ifndef __UTILS_H__
#define __UTILS_H__

/* Includes */
#include<mex.h>
#include<string.h>
#include<math.h>

/* In case debugging is necessary */
#define DEBUG 0
#define COL_MAJOR 10
#define ROW_MAJOR 11
#include<stdio.h>

/* Function prototypes */

/* 2D double matrix checking function */
void check_matrix(const mxArray *matrix, const char *context);

/* Sparse matrix checking function */
void check_sparse(const mxArray *matrix, const char *context);

/* 1D column vector checking function */
void check_vector(const mxArray *vector, const char *context);

/* Scalar checking function */
void check_scalar(const mxArray *scalar, const char *context);

/* Matrix printing function. Will work for vectors as well */
void print_matrix(double *matrix, int leading_dimension, int nrows, int ncols);

/* Integer matrix printing function. Will work for vectors as well */
void print_int_matrix(int *matrix, int leading_dimension, int nrows, int ncols);

/* Matlab matrix array printing function */
void print_mx_matrix(const mxArray *mx_matrix);

/* Non contiguous data copy utility */
void data_copy(double *dest, double *src, int nrows, int ncols, int *col_idx);

/* Find the maximum in an array */
int max(double *array, int nelem);

#endif
