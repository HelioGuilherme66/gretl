/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef GRETL_MATRIX_H
#define GRETL_MATRIX_H

/* #define LDEBUG 1 */

/* minimum value of diagonal element of R (as in X = QR) that counts
   as non-zero for the purpose of determining the rank of X */

#define R_DIAG_MIN 1.0e-8

typedef enum {
    GRETL_MOD_NONE = 0,
    GRETL_MOD_TRANSPOSE,
    GRETL_MOD_SQUARE,
    GRETL_MOD_CUMULATE
} GretlMatrixMod;

typedef struct _gretl_matrix gretl_matrix;
typedef struct _gretl_matrix gretl_vector;

struct _gretl_matrix {
    int rows;
    int cols;
    int t;
    double *val;
};

#define mdx(a,i,j)   ((j)*(a)->rows+(i))
#define mdxtr(a,i,j) ((i)*(a)->rows+(j))

/**
 * gretl_matrix_cols:
 * @m: matrix to query.
 * 
 * Gives the number of columns in @m. 
 */

#define gretl_matrix_cols(m) ((m == NULL)? 0 : m->cols)

/**
 * gretl_matrix_rows:
 * @m: matrix to query.
 * 
 * Gives the number of rows in @m. 
 */

#define gretl_matrix_rows(m) ((m == NULL)? 0 : m->rows)

/**
 * gretl_vector_get_length:
 * @v: vector to examine.
 * 
 * Gives the length of vector @v (without regard to whether
 * it is a row or column vector).
 */

#define gretl_vector_get_length(v) ((v == NULL)? 0 : \
                                    (v->cols == 1)? v->rows : \
                                    (v->rows == 1)? v->cols : 0)

/**
 * gretl_vector_alloc:
 * @i: number of columns.
 *
 * Allocates a new #gretl_vector with @i columns.
 */

#define gretl_vector_alloc(i) gretl_matrix_alloc(1,(i))

/**
 * gretl_column_vector_alloc:
 * @i: number of rows.
 *
 * Allocates a new column gretl_vector with @i rows.
 */

#define gretl_column_vector_alloc(i) gretl_matrix_alloc((i),1)

/**
 * gretl_vector_free:
 * @v: %gretl_vector to free.
 *
 * Frees the vector @v and its associated storage.
 */

#define gretl_vector_free(v) gretl_matrix_free(v)

/**
 * gretl_matrix_is_scalar:
 * @m: matrix to test.
 *
 * Gives 1 if @m is 1 x 1, else 0.
 */

#define gretl_matrix_is_scalar(m) ((m) != NULL && \
                                   (m)->rows == 1 && \
                                   (m)->cols == 1)

double gretl_matrix_get (const gretl_matrix *m, int i, int j);

double gretl_vector_get (const gretl_vector *v, int i);

int gretl_matrix_set (gretl_matrix *m, int i, int j, double x);

int gretl_vector_set (gretl_vector *v, int i, double x);

gretl_matrix *gretl_matrix_alloc (int rows, int cols);

gretl_matrix *gretl_matrix_reuse (gretl_matrix *m, int rows, int cols);

gretl_matrix *gretl_identity_matrix_new (int n);

gretl_matrix *gretl_zero_matrix_new (int r, int c);

gretl_matrix *gretl_unit_matrix_new (int r, int c);

gretl_matrix *gretl_null_matrix_new (void);

gretl_matrix *gretl_matrix_copy (const gretl_matrix *m);

int gretl_matrix_inscribe_I (gretl_matrix *m, int row, int col, int n);

gretl_matrix *gretl_matrix_copy_transpose (const gretl_matrix *m);

gretl_vector *gretl_column_vector_from_array (const double *x, 
					      int n, GretlMatrixMod mod);

gretl_vector *gretl_data_series_to_vector (const double **Z, int varno, 
					   int t1, int t2);

gretl_matrix *gretl_vector_from_array (const double *x, int n,
				       GretlMatrixMod mod);

gretl_matrix *gretl_matrix_from_2d_array (const double **X, 
					  int rows, int cols);

gretl_matrix *gretl_matrix_from_scalar (double x);

gretl_matrix *gretl_matrix_get_diagonal (const gretl_matrix *m, int *err);

double gretl_matrix_trace (const gretl_matrix *m, int *err);

int gretl_matrix_random_fill (gretl_matrix *m, int dist);

gretl_matrix *gretl_random_matrix_new (int r, int c, int dist);

double gretl_vector_mean (const gretl_vector *v);

double gretl_vector_variance (const gretl_vector *v);

void gretl_matrix_zero (gretl_matrix *m);

int gretl_matrix_zero_upper (gretl_matrix *m);

int gretl_matrix_zero_lower (gretl_matrix *m);

void gretl_matrix_fill (gretl_matrix *m, double x);

void gretl_matrix_multiply_by_scalar (gretl_matrix *m, double x);

int gretl_matrix_divide_by_scalar (gretl_matrix *m, double x);

void gretl_matrix_dot_pow (gretl_matrix *m, double x);

void gretl_matrix_free (gretl_matrix *m);

double *gretl_matrix_steal_data (gretl_matrix *m);

int gretl_matrix_copy_values (gretl_matrix *targ, 
			      const gretl_matrix *src);

int gretl_matrix_add_to (gretl_matrix *targ, const gretl_matrix *src);

int 
gretl_matrix_subtract_from (gretl_matrix *targ, const gretl_matrix *src);

int gretl_matrix_I_minus (gretl_matrix *m);

int gretl_matrix_transpose (gretl_matrix *targ, const gretl_matrix *src);

int gretl_square_matrix_transpose (gretl_matrix *m);

int gretl_matrix_add_self_transpose (gretl_matrix *m);

int 
gretl_matrix_vectorize (gretl_matrix *targ, const gretl_matrix *src);

int 
gretl_matrix_unvectorize (gretl_matrix *targ, const gretl_matrix *src);

int 
gretl_matrix_vectorize_h (gretl_matrix *targ, const gretl_matrix *src);

int 
gretl_matrix_unvectorize_h (gretl_matrix *targ, const gretl_matrix *src);

int gretl_matrix_inscribe_matrix (gretl_matrix *targ,
				  const gretl_matrix *src,
				  int row, int col,
				  GretlMatrixMod mod);

int gretl_matrix_extract_matrix (gretl_matrix *targ,
				 const gretl_matrix *src,
				 int row, int col,
				 GretlMatrixMod mod);

int gretl_matrix_multiply_mod (const gretl_matrix *a, GretlMatrixMod amod,
			       const gretl_matrix *b, GretlMatrixMod bmod,
			       gretl_matrix *c, GretlMatrixMod cmod);

int gretl_matrix_multiply (const gretl_matrix *a, const gretl_matrix *b,
			   gretl_matrix *c);

int
gretl_matrix_kronecker_product (const gretl_matrix *A, const gretl_matrix *B,
				gretl_matrix *K);

gretl_matrix *
gretl_matrix_kronecker_product_new (const gretl_matrix *A, 
				    const gretl_matrix *B);

double gretl_matrix_dot_product (const gretl_matrix *a, GretlMatrixMod amod,
				 const gretl_matrix *b, GretlMatrixMod bmod,
				 int *errp);

double gretl_vector_dot_product (const gretl_vector *a, const gretl_vector *b,
				 int *errp);

gretl_matrix *gretl_matrix_dot_multiply (const gretl_matrix *a, 
					 const gretl_matrix *b,
					 int *err);

gretl_matrix *gretl_matrix_dot_divide (const gretl_matrix *a, 
				       const gretl_matrix *b,
				       int *err);

gretl_matrix *gretl_matrix_row_sum (const gretl_matrix *m);

gretl_matrix *gretl_matrix_column_sum (const gretl_matrix *m);

gretl_matrix *gretl_matrix_row_mean (const gretl_matrix *m);

gretl_matrix *gretl_matrix_column_mean (const gretl_matrix *m);

double gretl_matrix_row_i_mean (const gretl_matrix *m, int row);

double gretl_matrix_column_j_mean (const gretl_matrix *m, int col);

void gretl_matrix_demean_by_row (gretl_matrix *m);

void gretl_matrix_demean_by_column (gretl_matrix *m);

gretl_matrix *gretl_matrix_vcv (gretl_matrix *m);

double gretl_matrix_determinant (gretl_matrix *a, int *err);

double gretl_matrix_log_determinant (gretl_matrix *a, int *err);

double gretl_matrix_log_abs_determinant (gretl_matrix *a, int *err);

double gretl_vcv_log_determinant (const gretl_matrix *m);

double gretl_matrix_one_norm (const gretl_matrix *m);

int gretl_LU_solve (gretl_matrix *a, gretl_vector *b);

int gretl_invert_general_matrix (gretl_matrix *a);

int gretl_invert_symmetric_indef_matrix (gretl_matrix *a);

int gretl_invert_symmetric_matrix (gretl_matrix *a);

int gretl_invert_symmetric_matrix2 (gretl_matrix *a, double *ldet);

int gretl_invert_packed_symmetric_matrix (gretl_matrix *v);

int gretl_invert_diagonal_matrix (gretl_matrix *a);

int gretl_invert_matrix (gretl_matrix *a);

int gretl_SVD_invert_matrix (gretl_matrix *a);

double gretl_symmetric_matrix_rcond (const gretl_matrix *m, int *err);

int gretl_eigen_sort (double *evals, gretl_matrix *evecs, int rank);

double *gretl_general_matrix_eigenvals (gretl_matrix *m,
					int eigenvecs, 
					int *err);

double *gretl_symmetric_matrix_eigenvals (gretl_matrix *m,
					  int eigenvecs, 
					  int *err);

gretl_matrix *gretl_matrix_right_nullspace (const gretl_matrix *M);

gretl_matrix *
gretl_matrix_col_concat (const gretl_matrix *a, const gretl_matrix *b,
			 int *err);

gretl_matrix *gretl_matrix_lag (gretl_matrix *m, int k, double missval);

int gretl_matrix_cholesky_decomp (gretl_matrix *a);

int gretl_matrix_QR_decomp (gretl_matrix *M, gretl_matrix *R);

int gretl_matrix_QR_rank (gretl_matrix *R, char **pmask, int *errp);

int gretl_matrix_rank (const gretl_matrix *a, int *err);

int gretl_matrix_ols (const gretl_vector *y, const gretl_matrix *X,
		      gretl_vector *b, gretl_matrix *vcv,
		      gretl_vector *uhat, double *s2);

int 
gretl_matrix_restricted_ols (const gretl_vector *y, const gretl_matrix *X,
			     const gretl_matrix *R, const gretl_vector *q,
			     gretl_vector *b, gretl_matrix *vcv,
			     double *s2);

int gretl_matrix_svd_ols (const gretl_vector *y, const gretl_matrix *X,
			  gretl_vector *b, gretl_matrix *vcv,
			  gretl_vector *uhat, double *s2);

double gretl_scalar_b_X_b (const gretl_vector *b, GretlMatrixMod bmod,
			   const gretl_matrix *X, int *errp);

gretl_matrix *
gretl_matrix_A_X_A (const gretl_matrix *A, GretlMatrixMod amod,
		    const gretl_matrix *X, int *errp);

int gretl_matrix_qform (const gretl_matrix *A, GretlMatrixMod amod,
			const gretl_matrix *X, gretl_matrix *C, 
			GretlMatrixMod cmod);

int
gretl_matrix_diagonal_sandwich (const gretl_vector *d, const gretl_matrix *X,
				gretl_matrix *DXD);

gretl_matrix *
gretl_vcv_matrix_from_model (MODEL *pmod, const char *select);

gretl_vector *
gretl_coeff_vector_from_model (const MODEL *pmod, const char *select);

void 
gretl_matrix_print_to_prn (const gretl_matrix *m, const char *msg, PRN *prn);

void gretl_matrix_print (const gretl_matrix *m, const char *msg);

void gretl_packed_matrix_print (const gretl_matrix *m, const char *msg);

void debug_print_matrix (const gretl_matrix *m, const char *msg);

void gretl_matrix_set_int (gretl_matrix *m, int t);

int gretl_matrix_get_int (const gretl_matrix *m);

int gretl_is_identity_matrix (const gretl_matrix *m);

int gretl_is_zero_matrix (const gretl_matrix *m);

int gretl_matrices_are_equal (const gretl_matrix *a, const gretl_matrix *b,
			      int *err);

gretl_matrix *
gretl_covariance_matrix_from_varlist (const int *list, const double **Z, 
				      const DATAINFO *pdinfo, 
				      gretl_matrix **means, int *errp);

int gretl_matrix_row_to_array (const gretl_matrix *m, int i, double *x);

gretl_matrix **gretl_matrix_array_alloc (int n);

gretl_matrix **
gretl_matrix_array_alloc_with_size (int n, int rows, int cols);

void gretl_matrix_array_free (gretl_matrix **A, int n);

gretl_matrix *gretl_matrix_data_subset (const int *list, const double **Z,
					int t1, int t2, const char *mask);

gretl_matrix *
gretl_matrix_data_subset_no_missing (const int *list, const double **Z,
				     int t1, int t2, int *err);

void lapack_mem_free (void);

#endif /* GRETL_MATRIX_H */
