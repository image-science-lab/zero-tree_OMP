/*
 * C function for simple implementation of OMP
 */

#include<omp_core.h>

/* Multi level dictionary OMP CS routine */
void mld_omp_cs(double Y[], double D[], double *X_Pr, size_t *X_Ir,
	   			size_t *X_Jc, int nrows, int ncols, double *A_Pr,
			   	size_t *A_Ir, size_t *A_Jc, int nrows_a, int ncols_a,
			   	int ndata, int blksize,	double thres)
{
	/* Variables to hold renormalized dictionary. */
	double *Dmod, *norms;
	int i;

	/* Allocate memory */
	Dmod = (double *)mxMalloc(nrows_a*ncols*sizeof(double));
	norms = (double *)mxMalloc(ncols*sizeof(double));

	/* Modify the dictionary. */
	sp_dotm(A_Pr, A_Ir, A_Jc, D, Dmod, nrows_a, ncols_a, ncols);

	/* Normalize the modified dictionary */
	normalize_mat(Dmod, norms, nrows_a, ncols);

	/* Now we can call the OMP routine with the modified dictionary */
	omp_mld(Y, Dmod, X_Pr, X_Ir, X_Jc, thres, nrows_a, ncols, ndata, blksize);

	/* Renormalize. This is easy work, as we have values and row indices in
	 * X_Pr and X_Ir respectively */
	for(i = 0; i < X_Jc[ndata]; i++)
		X_Pr[i] /= norms[X_Ir[i]];
}

/* Multi level dictionary OMP Routine */
void omp_mld(double Y[], double D[], double *X_Pr, size_t *X_Ir, size_t *X_Jc,
		 	 double thres, int nrows, int ncols, int ndata, int blksize)
{
	/* Create temporary holder variables. */
	double *y; 				/* Holds individual signals */
	double *x; 				/* Holds the non zero coefficients */
	int *support; 			/* Support set for the sparse signal */
	int dcnt; 				/* Data index counter */
	int i; 					/* Other counter variables */
	int x_cnt = 0; 			/* Count of the non-zero values for sparse X */
	int L = ncols/blksize; 	/* Number of blocks in the dictionary */

	/* Allocate memory */
	y = (double *)mxMalloc(nrows*sizeof(double));
	x = (double *)mxMalloc(L*sizeof(double));
	support = (int *)mxMalloc(L*sizeof(int));

	/* Individually solve for the signals */
	for(dcnt = 0; dcnt < ndata; dcnt++)
	{
		/* Store the beginning index of rows corresponding to this column. */
		X_Jc[dcnt] = x_cnt;
		/* Copy the individual signal to y */
		for(i = 0; i < nrows; i++)
			y[i] = Y[dcnt*nrows + i];

		/* Now invoke individual omp */
		_omp_mld(y, x, D, support, L, thres, nrows, ncols, blksize);

		/* Now copy the non zero coefficients to X */
		for(i = 0; i < L; i++)
		{
			/* If support[i] = -1, break. */
			if(support[i] == -1)
				break;
			/* We are using CSC format for sparse representation */
			X_Ir[x_cnt] = support[i];
			X_Pr[x_cnt] = x[i];
			x_cnt++;
		}
	}
	/* Finally add nzmax information. */
	X_Jc[ndata] = x_cnt;

	/* Done. Free the temporary arrays */
	mxFree(y);
	mxFree(x);
	mxFree(support);
}

/* Stand alone Multi level dictionary OMP */
void _omp_mld(double *y, double *x, double D[], int *support, int L,
	   		  double thres, int nrows, int ncols, int blksize)
{
	/* Holder variables for internal operations */
	double *dotprods, *dict_blk, *residue;
	int i, max_idx, blk_idx;
	int *blk_indices;

	/* Allocate memory */
	residue = (double *)mxMalloc(nrows*sizeof(double));
	dotprods = (double *)mxMalloc(blksize*sizeof(double));
	dict_blk = (double *)mxMalloc(nrows*blksize*sizeof(double));
	blk_indices = (int *)mxMalloc(blksize*sizeof(int));

	/* Initialize residue */
	for(i = 0; i < nrows; i++)
		residue[i] = y[i];

	for(blk_idx = 0; blk_idx < L; blk_idx++)
	{
		/* Check stopping condition */
		if(dotv(residue, residue, nrows) <= thres)
		{
			support[blk_idx] = -1;
			return;
		}

		/* Copy the dictionary block */
		for(i = 0; i < blksize; i++)
			blk_indices[i] = blk_idx*blksize + i;

		data_copy(dict_blk, D, nrows, blksize, blk_indices);

		/* Correlate and max */
		dotm('T', dict_blk, residue, dotprods, nrows, blksize, 0);
		max_idx = max(dotprods, blksize);

		/* Add it to support and update coefficients */
		support[blk_idx] = blk_idx*blksize + max_idx;
		x[blk_idx] = dotprods[max_idx];

		/* Update residue */
		for(i = 0; i < nrows; i++)
			residue[i] -= dict_blk[nrows*max_idx + i]*x[blk_idx];
	}

	/* Free the memory */
	mxFree(residue);
	mxFree(dotprods);
	mxFree(dict_blk);
	mxFree(blk_indices);
}

/* Constrained OMP CS routine */
void c_omp_cs(double Y[], double D1[], double *X1_Pr, size_t *X1_Ir,
	   		  size_t *X1_Jc, int k1, int nrows1, int ncols1, double D0[],
			  double *X0_Pr, size_t *X0_Ir, size_t *X0_Jc, int k0, int nrows0,
			  int ncols0, double *A_Pr, size_t *A_Ir, size_t *A_Jc, int nrows_a,
			  int ncols_a, int ndata, double *U_Pr, size_t *U_Ir, size_t *U_Jc,
			  double thres, int hf)
{
	/* Variables to hold renormalized dictionary. */
	double *D1mod, *D1mod1, *D0mod, *norms1, *norms0;
	int i;

	/* Allocate memory */
	D1mod1 = (double *)mxMalloc(nrows0*ncols1*sizeof(double));
	D1mod = (double *)mxMalloc(nrows_a*ncols1*sizeof(double));
	D0mod = (double *)mxMalloc(nrows_a*ncols0*sizeof(double));
	norms1 = (double *)mxMalloc(ncols1*sizeof(double));
	norms0 = (double *)mxMalloc(ncols0*sizeof(double));

	/* Modify the low resolution dictionary in two steps */
	sp_dotm(U_Pr, U_Ir, U_Jc, D1, D1mod1, nrows0, nrows1, ncols1);
	sp_dotm(A_Pr, A_Ir, A_Jc, D1mod1, D1mod, nrows_a, ncols_a, ncols1);

	/* Modify the high resolution dictionary now*/
	sp_dotm(A_Pr, A_Ir, A_Jc, D0, D0mod, nrows_a, ncols_a, ncols0);

	/* Normalize the dictionary now */
	normalize_mat(D1mod, norms1, nrows_a, ncols1);
	normalize_mat(D0mod, norms0, nrows_a, ncols0);

	/* Call constrained OMP next */
	/* Note: R_Pr, R_Ir, R_Jc are not used when use_resizer=0. Hence, we send
	 * some pointer as dummy. We are sending U here */

	if(hf)
		c_omp_hf(Y, D1mod, X1_Pr, X1_Ir, X1_Jc, k1, nrows_a, ncols1, D0mod,
			   	 X0_Pr, X0_Ir, X0_Jc, k0, nrows_a, ncols0, ndata, U_Pr, U_Ir,
				 U_Jc, U_Pr, U_Ir, U_Jc, thres, 0);
	else
		c_omp(Y, D1mod, X1_Pr, X1_Ir, X1_Jc, k1, nrows_a, ncols1, D0mod, X0_Pr,
		  	  X0_Ir, X0_Jc, k0, nrows_a, ncols0, ndata, U_Pr, U_Ir, U_Jc, thres,
			  0);

	/* Renormalize sparse vector */
	for(i = 0; i < X1_Jc[ndata]; i++)
		X1_Pr[i] /= norms1[X1_Ir[i]];

	for(i = 0; i < X0_Jc[ndata]; i++)
		X0_Pr[i] /= norms0[X0_Ir[i]];
}

/* OMP CS routine */
void omp_cs(double Y[], double D[], double *X_Pr, size_t *X_Ir, size_t *X_Jc,
	   		int k, int nrows, int ncols, double *A_Pr, size_t *A_Ir,
			size_t *A_Jc, int nrows_a, int ncols_a, int ndata, double thres)
{
	/* Variables to hold renormalized dictionary. */
	double *Dmod, *norms;
	int i;

	/* Allocate memory */
	Dmod = (double *)mxMalloc(nrows_a*ncols*sizeof(double));
	norms = (double *)mxMalloc(ncols*sizeof(double));

	/* Modify the dictionary. */
	sp_dotm(A_Pr, A_Ir, A_Jc, D, Dmod, nrows_a, ncols_a, ncols);

	/* Normalize the modified dictionary */
	normalize_mat(Dmod, norms, nrows_a, ncols);

	/* Now we can call the OMP routine with the modified dictionary */
	omp(Y, Dmod, X_Pr, X_Ir, X_Jc, k, thres, nrows_a, ncols, ndata);

	/* Renormalize. This is easy work, as we have values and row indices in
	 * X_Pr and X_Ir respectively */
	for(i = 0; i < X_Jc[ndata]; i++)
		X_Pr[i] /= norms[X_Ir[i]];
}

/* Constrained OMP routine for three scales*/
void c_omp3(double Y[], double D2[], double X2[], int k2, int nrows2, int ncols2,
	double D1[], double X1[], int k1, int nrows1, int ncols1,
	double D0[], double X0[], int k0, int nrows0, int ncols0, int ndata,
	double *resizer1, double *resizer0, double thres)
{
	/* Temporary holder variables. */
	double *y, *y1, *y2, *x2, *x1, *x0; 	/* Signal and the sparse vector coefficients */
	int *sup2, *sup1, *sup0, *csub1, *csub0;/* Support sets for L1 and L0 */
	int dcnt, i, j; 						/* Data and dummy counter */

	/* Memory allocation */
	y = (double *)mxMalloc(nrows0*sizeof(double));
	y1 = (double *)mxMalloc(nrows1*sizeof(double));
	y2 = (double *)mxMalloc(nrows2*sizeof(double));

	x2 = (double *)mxMalloc(k2*sizeof(double));
	x1 = (double *)mxMalloc(k1*sizeof(double));
	x0 = (double *)mxMalloc(k0*sizeof(double));

	sup2 = (int *)mxMalloc(k2*sizeof(int));
	sup1 = (int *)mxMalloc(k1*sizeof(int));
	sup0 = (int *)mxMalloc(k0*sizeof(int));

	csub1 = (int *)mxMalloc(k2*(ncols1 / ncols2)*sizeof(int));
	csub0 = (int *)mxMalloc(k1*(ncols0 / ncols1)*sizeof(int));

	/* Solve for individual signals */
	for (dcnt = 0; dcnt < ndata; dcnt++)
	{
		/* Copy individual signal */
		for (i = 0; i < nrows0; i++)
			y[i] = Y[dcnt*nrows0 + i];

		/* Create a low resolution signal */
		dotm('N', resizer0, y, y1, nrows1, nrows0, 0);

		/* Create further low resolution signal */
		dotm('N', resizer1, y1, y2, nrows2, nrows1, 0);

		/* Execute OMP for the least resolution dictionary*/
		_omp(y2, x2, D2, sup2, k2, thres, nrows2, ncols2);

		/* Execute OMP for the next low resolution dictionary */
		_c_omp(y1, x1, D1, k1, sup1, thres, nrows1, ncols1, k2,
			ncols1 / ncols2, sup2, csub1);

		/* Tricky stuff. We need to transfer the support from constraint set*/
		for (i = 0; i < k1; i++)
		{
			sup1[i] = csub1[sup1[i]];
		}

		/* Now that we have  the constrain set x1, execute constrained omp */
		_c_omp(y, x0, D0, k0, sup0, thres, nrows0, ncols0, k1,
			ncols0 / ncols1, sup1, csub0);

		/* Done executing OMP and cOMP. Create the sparse vector */
		for (i = 0; i < k2; i++)
		{
			if (sup2[i] == -1)
				break;
			X2[dcnt*ncols2 + sup2[i]] = x2[i];
		}

		for (i = 0; i < k1; i++)
		{
			if (sup1[i] == -1)
				break;
			X1[dcnt*ncols1 + sup1[i]] = x1[i];
		}

		for (i = 0; i < k0; i++)
		{
			if (sup0[i] == -1)
				break;
			X0[dcnt*ncols0 + csub0[sup0[i]]] = x0[i];
		}
	}

	/* Free the temporary variables */
	mxFree(y); mxFree(y1); mxFree(y2);
	mxFree(x2);  mxFree(x1); mxFree(x0);
	mxFree(sup2);  mxFree(sup1); mxFree(sup0);
	mxFree(csub1); mxFree(csub0);
}

/* Constrained OMP routine */
void c_omp(double Y[], double D1[], double *X1_Pr, size_t *X1_Ir, size_t *X1_Jc,
	   	   int k1, int nrows1, int ncols1, double D0[], double *X0_Pr,
		   size_t *X0_Ir, size_t *X0_Jc, int k0, int nrows0, int ncols0,
		   int ndata, double *R_Pr, size_t *R_Ir, size_t *R_Jc, double thres,
		   int use_resizer)
{
	/* Temporary holder variables. */
	double *y, *y1, *x1, *x0; 	/* Signal and the sparse vector coefficients */
	int *sup1, *sup0, *csub; 	/* Support sets for L1 and L0 */
	int dcnt, i, j; 			/* Data and dummy counter */
	int x1_cnt, x0_cnt; 		/* Counters for sparse data */

	x1_cnt = 0;
	x0_cnt = 0;

	/* Memory allocation */
	y = (double *)mxMalloc(nrows0*sizeof(double));
	y1 = (double *)mxMalloc(nrows1*sizeof(double));

	x1 = (double *)mxMalloc(k1*sizeof(double));
	x0 = (double *)mxMalloc(k0*sizeof(double));

	sup1 = (int *)mxMalloc(k1*sizeof(int));
	sup0 = (int *)mxMalloc(k0*sizeof(int));
	csub = (int *)mxMalloc(k1*(ncols0/ncols1)*sizeof(int));

	/* Solve for individual signals */
	for(dcnt = 0; dcnt < ndata; dcnt++)
	{
		/* Update Jc for X1 and X0 */
		X1_Jc[dcnt] = x1_cnt;
		X0_Jc[dcnt] = x0_cnt;

		/* Copy individual signal */
		for(i = 0; i < nrows0; i++)
		{
			y[i] = Y[dcnt*nrows0 + i];
			if (!use_resizer)
				y1[i] = Y[dcnt*nrows0 + i];
		}

		/* Create a low resolution signal if necessary*/
		if (use_resizer)
			sp_dotv(R_Pr, R_Ir, R_Jc, y, y1, nrows1, nrows0);

		/* Else, just copy the output, which is done above */

		/* Execute OMP for the low resolution dictionary */
		_omp(y1, x1, D1, sup1, k1,  thres, nrows1, ncols1);

		/* Now that we have  the constrain set x1, execute constrained omp */
		_c_omp(y, x0, D0, k0, sup0, thres, nrows0, ncols0, k1,
				   ncols0/ncols1, sup1, csub);

		/* Done executing OMP and cOMP. Create the sparse vector */
		for(i = 0; i < k1; i++)
		{
			if (sup1[i] == -1)
				break;
			/* Using CSC format for sparse representation. */
			X1_Ir[x1_cnt] = sup1[i];
			X1_Pr[x1_cnt] = x1[i];
			x1_cnt++;
		}

		for(i = 0; i < k0; i++)
		{
			if(sup0[i] == -1)
				break;
			X0_Ir[x0_cnt] = csub[sup0[i]];
			X0_Pr[x0_cnt] = x0[i];
			x0_cnt++;
		}
	}
	/* Finally add nzmax information. */
	X1_Jc[ndata] = x1_cnt;
	X0_Jc[ndata] = x0_cnt;

	/* Free the temporary variables */
	mxFree(y); mxFree(y1);
	mxFree(x1); mxFree(x0);
	mxFree(sup1); mxFree(sup0); mxFree(csub);
}

/* Constrained OMP routine */
void c_omp_hf(double Y[], double D1[], double *X1_Pr, size_t *X1_Ir,
	   	      size_t *X1_Jc, int k1, int nrows1, int ncols1, double D0[],
			  double *X0_Pr, size_t *X0_Ir, size_t *X0_Jc, int k0, int nrows0,
			  int ncols0, int ndata, double *R_Pr, size_t *R_Ir, size_t *R_Jc,
			  double *U_Pr, size_t *U_Ir, size_t *U_Jc,
			  double thres, int use_resizer)
{
	/* Temporary holder variables. */
	double *y, *y1, *x1, *x0; 	/* Signal and the sparse vector coefficients */
	int *sup1, *sup0, *csub; 	/* Support sets for L1 and L0 */
	int dcnt, i, j; 			/* Data and dummy counter */
	int x1_cnt, x0_cnt; 		/* Counters for sparse data */
	double *y_lf; 				/* Low resolution approximation */

	x1_cnt = 0;
	x0_cnt = 0;

	/* Memory allocation */
	y = (double *)mxMalloc(nrows0*sizeof(double));
	y_lf = (double *)mxMalloc(nrows0*sizeof(double));
	y1 = (double *)mxMalloc(nrows1*sizeof(double));

	x1 = (double *)mxMalloc(k1*sizeof(double));
	x0 = (double *)mxMalloc(k0*sizeof(double));

	sup1 = (int *)mxMalloc(k1*sizeof(int));
	sup0 = (int *)mxMalloc(k0*sizeof(int));
	csub = (int *)mxMalloc(k1*(ncols0/ncols1)*sizeof(int));

	/* Solve for individual signals */
	for(dcnt = 0; dcnt < ndata; dcnt++)
	{
		/* Update Jc for X1 and X0 */
		X1_Jc[dcnt] = x1_cnt;
		X0_Jc[dcnt] = x0_cnt;

		/* Copy individual signal */
		for(i = 0; i < nrows0; i++)
		{
			y[i] = Y[dcnt*nrows0 + i];
			if (!use_resizer)
				y1[i] = Y[dcnt*nrows0 + i];
		}

		/* Create a low resolution signal if necessary*/
		if (use_resizer)
			sp_dotv(R_Pr, R_Ir, R_Jc, y, y1, nrows1, nrows0);

		/* Else, just copy the output, which is done above */

		/* Execute OMP for the low resolution dictionary */
		_omp(y1, x1, D1, sup1, k1,  thres, nrows1, ncols1);
		
		/* Compute high frequency difference */
		Dxmul(y1, x1, sup1, D1, nrows1, k1);
		print_matrix(y1, COL_MAJOR, nrows1, 1);

		/* Create a high resolution signal if necessary */
		if (use_resizer)
		{
			sp_dotv(U_Pr, U_Ir, U_Jc, y1, y_lf, nrows0, nrows1);
			for(i = 0; i < nrows0; i++)
				y[i] -= y_lf[i];
		}
		else
		{
			for(i = 0; i < nrows0; i++)
				y[i] -= y1[i];
		}

		/* Now that we have  the constrain set x1, execute constrained omp */
		_c_omp(y, x0, D0, k0, sup0, thres, nrows0, ncols0, k1,
				   ncols0/ncols1, sup1, csub);

		/* Done executing OMP and cOMP. Create the sparse vector */
		for(i = 0; i < k1; i++)
		{
			if (sup1[i] == -1)
				break;
			/* Using CSC format for sparse representation. */
			X1_Ir[x1_cnt] = sup1[i];
			X1_Pr[x1_cnt] = x1[i];
			x1_cnt++;
		}

		for(i = 0; i < k0; i++)
		{
			if(sup0[i] == -1)
				break;
			X0_Ir[x0_cnt] = csub[sup0[i]];
			X0_Pr[x0_cnt] = x0[i];
			x0_cnt++;
		}
	}
	/* Finally add nzmax information. */
	X1_Jc[ndata] = x1_cnt;
	X0_Jc[ndata] = x0_cnt;

	/* Free the temporary variables */
	mxFree(y); mxFree(y1); mxFree(y_lf);
	mxFree(x1); mxFree(x0);
	mxFree(sup1); mxFree(sup0); mxFree(csub);
}

/* Standalone constrained OMP */
void c_omp_std(double Y[], double X1[], int l0_resolve,	double D0[],
	   	double X0[], int k0, int nrows0, int ncols0, int ndata,
	   	int ncols1)
{
	/* Temporary holder variables. */
	double *y, *y1, *x1, *x0; 	/* Signal and the sparse vector coefficients */
	int *sup1, *sup0, *csub; 	/* Support sets for L1 and L0 */
	int dcnt, i, j; 			/* Data and dummy counter */
	int k1; 					/* L1 sparsity, which we will estimate */

	double thres = 1e-4;
	/* Memory allocation */
	y = (double *)mxMalloc(nrows0*sizeof(double));

	x0 = (double *)mxMalloc(k0*sizeof(double));

	sup1 = (int *)mxMalloc(20*sizeof(int));
	sup0 = (int *)mxMalloc(k0*sizeof(int));

	/* Solve for individual signals */
	for(dcnt = 0; dcnt < ndata; dcnt++)
	{
		/* Copy individual signal */
		for(i = 0; i < nrows0; i++)
			y[i] = Y[dcnt*nrows0 + i];

		/* Find the support set*/
		j = 0;
		for(i = 0; i < ncols1; i++)
		{
			if (X1[dcnt*ncols1 + j] != 0)
			{
				sup1[j] = i;
				j++;
			}
		}

		/* Now allocate space for subset of indices */
		k1 = j;
		csub = (int *)mxMalloc(k1*l0_resolve*sizeof(int));

		/* Now that we have  the constrain set x1, execute constrained omp */
		_c_omp(y, x0, D0, k0, sup0, thres, nrows0, ncols0, k1,
				   l0_resolve, sup1, csub);

		for(i = 0; i < k0; i++)
		{
			if(sup0[i] == -1)
				break;
			X0[dcnt*ncols0 + csub[sup0[i]]] = x0[i];
		}
	}

	/* Free the temporary variables */
	mxFree(y);
	mxFree(x0);
	mxFree(sup1); mxFree(sup0); mxFree(csub);
}

/* Single signal constrained OMP routine */
void _c_omp(double *y, double *x, double D[], int k0, int *sup0, double thres,
	   	int nrows, int ncols, int k1, int l0_resolve, int *sup1, int *csub)
{
	/* Variables to store sub dictionary and constrain support set*/
	double *d_sub;
	int n_idx;

	/* Create the sub-support set */
	n_idx = create_constraint_set(csub, sup1, k1, l0_resolve);

	/* Allocate memory for sub-dictionary */
	d_sub = (double *)mxMalloc(nrows*n_idx*sizeof(double));

	/* Copy the sub dictionary now */
	data_copy(d_sub, D, nrows, n_idx, csub);

	/* Now run OMP on the sub set */
	_omp(y, x, d_sub, sup0, k0, thres, nrows, n_idx);

	/* Free the sub dictionary */
	mxFree(d_sub);
}


/* OMP Routine */
void omp(double Y[], double D[], double *X_Pr, size_t *X_Ir, size_t *X_Jc,
		 int k, double thres, int nrows, int ncols, int ndata)
{
	/* Create temporary holder variables. */
	double *y; 				/* Holds individual signals */
	double *x; 				/* Holds the non zero coefficients */
	int *support; 			/* Support set for the sparse signal */
	int dcnt; 				/* Data index counter */
	int i; 					/* Other counter variables */
	int x_cnt = 0; 			/* Count of the non-zero values for sparse X */

	/* Allocate memory */
	y = (double *)mxMalloc(nrows*sizeof(double));
	x = (double *)mxMalloc(k*sizeof(double));
	support = (int *)mxMalloc(k*sizeof(int));

	/* Individually solve for the signals */
	for(dcnt = 0; dcnt < ndata; dcnt++)
	{
		/* Store the beginning index of rows corresponding to this column. */
		X_Jc[dcnt] = x_cnt;
		/* Copy the individual signal to y */
		for(i = 0; i < nrows; i++)
			y[i] = Y[dcnt*nrows + i];

		/* Now invoke individual omp */
		_omp(y, x, D, support, k, thres, nrows, ncols);

		/* Now copy the non zero coefficients to X */
		for(i = 0; i < k; i++)
		{
			/* If support[i] = -1, break. */
			if(support[i] == -1)
				break;
			/* We are using CSC format for sparse representation */
			X_Ir[x_cnt] = support[i];
			X_Pr[x_cnt] = x[i];
			x_cnt++;
		}
	}
	/* Finally add nzmax information. */
	X_Jc[ndata] = x_cnt;

	/* Done. Free the temporary arrays */
	mxFree(y);
	mxFree(x);
	mxFree(support);
}


/* Single signal OMP */
void _omp(double *y, double *x, double D[], int *support, int k,
	   	double thres, int nrows, int ncols)
{
	/* Holder variables for internal operations */
	double *residue, *dotprods, *temp_dict, *temp_x;
	double residue_norm;
	int i, j, max_idx;
	int *temp_support;

	/* Allocate memory */
	residue = (double *)mxMalloc(nrows*sizeof(double));
	dotprods = (double *)mxMalloc(ncols*sizeof(double));

	/* Initially, residue = y */
	memcpy(residue, y, nrows*sizeof(double));

	for(i = 0; i < k; i++)
	{
		/* If norm of residue vector is less than thres, break. */
		if (dotv(residue, residue, nrows) <= thres)
		{
			/* Set support set index = -1 to flag termination. */
			support[i] = -1;
			break;
		}

		/* Find D'*residue */
		dotm('T', D, residue, dotprods, nrows, ncols, 0);

		/* Find the index of the maximum of absolute of the dot products */
		max_idx = max(dotprods, ncols);

		support[i] = max_idx;

		/* Now solve the least squares. For this, we need to create a
		 * temporary subdictionary */
		temp_dict = (double *)mxMalloc(nrows*(i+1)*sizeof(double));
		temp_support = (int *)mxMalloc((i+1)*sizeof(int));
		temp_x = (double *)mxMalloc((i+1)*sizeof(double));

		/* copy support indices to temporary support */
		for (j = 0; j < i+1; j++)
		{
			temp_support[j] = support[j];
		}

		/* Copy the support set atoms to the temporary dictionary */
		data_copy(temp_dict, D, nrows, i+1, temp_support);

		/* Now get a least squares solution for the system */
		lstsq(temp_dict, y, temp_x, nrows, i+1, 1);

		/* update the residue vector */
		/* A bit roundabout way of doing stuff, but should work */
		memcpy(residue, y, nrows*sizeof(double));
		dotm('N', temp_dict, temp_x, residue, nrows, i+1, -1);

		/* Copy data from temp_x to x */
		for (j = 0; j < i+1; j++)
			x[j] = temp_x[j];

		/* Free all the temporary variables */
		mxFree(temp_dict);
		mxFree(temp_support);
		mxFree(temp_x);
	}

	/* Clear residue and dotprods */
	mxFree(residue);
	mxFree(dotprods);

}

/* Constrained set creation */
int create_constraint_set(int *c_sub, int *x1, int k1, int l0_resolve)
{
	int i, j;
	for(i = 0; i < k1; i++){
		if (x1[i] == -1)
		{
			c_sub[i*l0_resolve] = -1;
			return i*l0_resolve;
		}
		for(j = 0; j < l0_resolve; j++)
			c_sub[i*l0_resolve + j] = x1[i]*l0_resolve + j;
	}
	return k1*l0_resolve;
}
