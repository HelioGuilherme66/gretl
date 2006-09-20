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

#ifndef USERMAT_H_
#define USERMAT_H_

enum {
    SEL_RANGE,
    SEL_MATRIX,
    SEL_DIAG,
    SEL_ALL,
    SEL_NULL
};

typedef struct matrix_subspec_ matrix_subspec;

union msel {
    int range[2];
    gretl_matrix *m;
};

struct matrix_subspec_ {
    int type[2];
    union msel sel[2];
};

int n_user_matrices (void);

gretl_matrix *user_matrix_by_index (int i, const char **name);

const char *get_matrix_name_by_index (int idx);

gretl_matrix *get_matrix_by_name (const char *name);

gretl_matrix *get_matrix_by_name_at_level (const char *name, int level,
					   const DATAINFO *pdinfo);

int user_matrix_add (gretl_matrix *M, const char *name);

int user_matrix_destroy (const char *name, PRN *prn);

int user_matrix_replace_matrix (const char *name, gretl_matrix *M);

int user_matrix_replace_submatrix (const char *name, gretl_matrix *M,
				   matrix_subspec *spec);

int add_or_replace_user_matrix (gretl_matrix *M, const char *name);

int copy_named_matrix_as (const char *orig, const char *new);

int user_matrix_set_name_and_level (const gretl_matrix *M, char *name, 
				    int level);

void destroy_user_matrices (void);

int destroy_user_matrices_at_level (int level);

double user_matrix_get_determinant (const gretl_matrix *m, int *err);

double user_matrix_get_log_determinant (const gretl_matrix *m, int *err);

gretl_matrix *user_matrix_get_inverse (const gretl_matrix *m);

gretl_matrix *user_matrix_cholesky_decomp (const gretl_matrix *m);

gretl_matrix *user_matrix_column_demean (const gretl_matrix *m);

gretl_matrix *user_matrix_vec (const gretl_matrix *m);

gretl_matrix *user_matrix_vech (const gretl_matrix *m, int *err);

gretl_matrix *user_matrix_unvech (const gretl_matrix *m, int *err);

gretl_matrix *
user_matrix_QR_decomp (const gretl_matrix *m, const char *rname, 
		       int *err);

gretl_matrix *
user_matrix_eigen_analysis (const gretl_matrix *m, const char *rname, int symm,
			    int *err);

gretl_matrix *matrix_get_submatrix (const gretl_matrix *M, 
				    matrix_subspec *spec,
				    int *err);

gretl_matrix *user_matrix_get_submatrix (const char *name, 
					 matrix_subspec *spec,
					 int *err);

#endif /* USERMAT_H_ */
