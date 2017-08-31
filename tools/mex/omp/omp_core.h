/*
 * Header file for Core implementation of OMP
 */

#ifndef __OMP_CORE_H__
#define __OMP_CORE_H__

#include<mex.h>
#include<utils.h>
#include<linalg.h>
#include<sparse.h>

/* Multi level dictionary OMP CS routine */
void mld_omp_cs(double Y[], double D[], double *X_Pr, size_t *X_Ir,
	   			size_t *X_Jc, int nrows, int ncols, double *A_Pr,
			   	size_t *A_Ir, size_t *A_Jc, int nrows_a, int ncols_a,
			   	int ndata, int blksize,	double thres);

/* Multi level dictionary OMP Routine */
void omp_mld(double Y[], double D[], double *X_Pr, size_t *X_Ir, size_t *X_Jc,
		 	 double thres, int nrows, int ncols, int ndata, int blksize);

/* Standalone Multi level dictionary OMP */
void _omp_mld(double *y, double *x, double D[], int *support, int L,
	   		  double thres, int nrows, int ncols, int blksize);

/* Constrained OMP CS routine */
void c_omp_cs(double Y[], double D1[], double *X1_Pr, size_t *X1_Ir,
	   		  size_t *X1_Jc, int k1, int nrows1, int ncols1, double D0[],
			  double *X0_Pr, size_t *X0_Ir, size_t *X0_Jc, int k0, int nrows0,
			  int ncols0, double *A_Pr, size_t *A_Ir, size_t *A_Jc, int nrows_a,
			  int ncols_a, int ndata, double *U_Pr, size_t *U_Ir, size_t *U_Jc,
			  double thres, int hf);

/* Constrained OMP routine */
void c_omp(double Y[], double D1[], double *X1_Pr, size_t *X1_Ir, size_t *X1_Jc,
	   		int k1, int nrows1, int ncols1, double D0[], double *X0_Pr,
		   	size_t *X0_Ir, size_t *X0_Jc, int k0, int nrows0, int ncols0,
		   	int ndata, double *R_Pr, size_t *R_Ir, size_t *R_Jc, double thres,
		   	int use_resizer);

/* Constrained OMP routine with high frequency dictionary*/
void c_omp_hf(double Y[], double D1[], double *X1_Pr, size_t *X1_Ir,
	   	      size_t *X1_Jc, int k1, int nrows1, int ncols1, double D0[],
			  double *X0_Pr, size_t *X0_Ir, size_t *X0_Jc, int k0, int nrows0,
			  int ncols0, int ndata, double *R_Pr, size_t *R_Ir, size_t *R_Jc,
			  double *U_Pr, size_t *U_Ir, size_t *U_Jc,
			  double thres, int use_resizer);

/* Single signal constrained OMP routine */
void _c_omp(double *y, double *x, double D[], int k0, int *sup0, double thres,
   	int nrows, int ncols, int k1, int l0_resolve, int *sup1, int *csub);

/* OMP CS routine */
void omp_cs(double Y[], double D[], double *X_Pr, size_t *X_Ir, size_t *X_Jc,
	   			int k, int nrows, int ncols, double *A_Pr, size_t *A_Ir,
					size_t *A_Jc, int nrows_a, int ncols_a, int ndata, double thres);

/* OMP routine */
void omp(double Y[], double D[], double *X_Pr, size_t *X_Ir, size_t *X_Jc,
		 int k, double thres, int nrows, int ncols, int ndata);

/* Single signal omp */
void _omp(double *y, double *x, double D[], int *support, int k,
		double thres, int nrows, int ncols);

/* Standalone constrained OMP */
void c_omp_std(double Y[], double X1[], int l0_resolve,
		double D0[], double X0[], int k0, int nrows0, int ncols0, int ndata,
	   	int ncols1);

/* Create constrained subset */
int create_constraint_set(int *c_sub, int *x1, int k1, int l0_resolve);

#endif
