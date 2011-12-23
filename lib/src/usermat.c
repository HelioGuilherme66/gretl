/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "gretl_matrix.h"
#include "matrix_extra.h"
#include "gretl_func.h"
#include "libset.h"
#include "usermat.h"
#include "gretl_xml.h"
#include "genparse.h"
#include "gretl_bundle.h"

#define MDEBUG 0
#define LEVEL_AUTO -1

typedef enum {
    UM_PRIVATE = 1 << 0,
    UM_SHELL   = 1 << 1
} UMFlags;

struct user_matrix_ {
    gretl_matrix *M;
    int level;
    UMFlags flags;
    char name[VNAMELEN];
};

static user_matrix **matrices;
static int n_matrices;

#define matrix_is_private(u) ((u->flags & UM_PRIVATE) || *u->name == '$')
#define matrix_is_shell(u) (u->flags & UM_SHELL)

#if MDEBUG > 1
static void print_matrix_stack (const char *msg)
{
    int i;

    fprintf(stderr, "\nmatrix stack, %s:\n", msg);
    for (i=0; i<n_matrices; i++) {
	if (matrices[i] == NULL) {
	    fprintf(stderr, " %d: NULL\n", i);
	} else {
	    fprintf(stderr, " %d: '%s' at %p (%d x %d)\n",
		    i, matrices[i]->name, (void *) matrices[i]->M,
		    matrices[i]->M->rows, matrices[i]->M->cols);
	}
    }
    fputc('\n', stderr);
}
#endif

int n_user_matrices (void)
{
    return n_matrices;
}

const char *get_matrix_name_by_index (int idx)
{
    if (idx >= 0 && idx < n_matrices) {
	return matrices[idx]->name;
    } else {
	return NULL;
    }
}

static user_matrix *user_matrix_new (gretl_matrix *M, const char *name)
{
    user_matrix *u = malloc(sizeof *u);

    if (u == NULL) {
	return NULL;
    }

    u->M = M;
    u->level = gretl_function_depth();
    u->flags = 0;
    *u->name = '\0';
    strncat(u->name, name, VNAMELEN - 1);

    return u;
}

static int matrix_is_saved (const gretl_matrix *m)
{
    int i;

    for (i=0; i<n_matrices; i++) {
	if (m == matrices[i]->M) {
	    return 1;
	}
    }

    if (data_is_bundled((void *) m)) {
	return 1;
    }

    return 0;
}

/* callback for adding or deleting icons representing 
   matrices in the GUI session window */

static void (*matrix_add_delete_callback)(const char *name,
					  int action);

/**
 * set_matrix_add_delete_callback:
 * @callback: function function to out in place.
 *
 * Sets the callback function to be invoked when a user-defined
 * matrix is added to or removed from the stack of saved objects.  
 * Intended for synchronizing the GUI program with the saved object
 * state.
 */

void set_matrix_add_delete_callback (void (*callback))
{
    matrix_add_delete_callback = callback; 
}

static user_matrix *real_user_matrix_add (gretl_matrix *M, 
					  const char *name,
					  int callback_ok)
{
    user_matrix **tmp;
    user_matrix *u;
    int n = n_matrices;

    if (M == NULL) {
	return 0;
    }

    tmp = realloc(matrices, (n + 1) * sizeof *tmp);
    if (tmp == NULL) {
	return NULL;
    } else {
	matrices = tmp;
    }

    if (matrix_is_saved(M)) {
	/* ensure uniqueness of matrix pointers */
	gretl_matrix *Mcpy = gretl_matrix_copy(M);

	if (Mcpy == NULL) {
	    return NULL;
	}
	M = Mcpy;
    }

    matrices[n] = u = user_matrix_new(M, name);

    if (matrices[n] == NULL) {
	return NULL;
    }

    n_matrices++;

#if MDEBUG
    fprintf(stderr, "add_user_matrix: allocated '%s' at %p (M at %p, %dx%d)\n",
	    name, (void *) matrices[n], (void *) M,
	    gretl_matrix_rows(M), gretl_matrix_cols(M));
#endif

    if (callback_ok && matrix_add_delete_callback != NULL && 
	gretl_function_depth() == 0) {
	(*matrix_add_delete_callback)(name, USER_MATRIX_ADD);
    }

    return u;
}

int user_matrix_add (gretl_matrix *M, const char *name)
{
    user_matrix *u = real_user_matrix_add(M, name, 1);

    return (u == NULL)? E_ALLOC : 0;
}

int private_matrix_add (gretl_matrix *M, const char *name)
{
    user_matrix *u = real_user_matrix_add(M, name, 0);

    if (u == NULL) {
	return E_ALLOC;
    } else {
	u->flags |= UM_PRIVATE;
	return 0;
    }
}

static int 
matrix_insert_diagonal (gretl_matrix *M, const gretl_matrix *S,
			int mr, int mc)
{
    int i, n = gretl_vector_get_length(S);
    int k = (mr < mc)? mr : mc;

    if (n != k) {
	return E_NONCONF;
    }

    for (i=0; i<n; i++) {
	gretl_matrix_set(M, i, i, S->val[i]);
    }
    
    return 0;
}

/* If level is LEVEL_AUTO, search at the current function execution
   depth, otherwise search at the function execution depth given by
   slevel.
*/

static gretl_matrix *real_get_matrix_by_name (const char *name, 
					      int level)
{
    int i;

    if (level == LEVEL_AUTO) {
	level = gretl_function_depth();
    } 

    for (i=0; i<n_matrices; i++) {
	if (matrices[i]->level == level && 
	    !strcmp(name, matrices[i]->name)) {
	    return matrices[i]->M;
	}
    }

    return NULL;
}

user_matrix *get_user_matrix_by_name (const char *name)
{
    int i, level = gretl_function_depth();

    for (i=0; i<n_matrices; i++) {
	if (matrices[i]->level == level && 
	    !strcmp(name, matrices[i]->name)) {
	    return matrices[i];
	}
    }

    return NULL;
}

user_matrix *get_user_matrix_by_index (int idx)
{
    if (idx >= 0 && idx < n_matrices) {
	return matrices[idx];
    } else {
	return NULL;
    }
}

static void usermat_destroy_matrix (user_matrix *u)
{
    if (!data_is_bundled(u->M)) {
	gretl_matrix_free(u->M);
    }
}

int user_matrix_replace_matrix (user_matrix *u, gretl_matrix *M)
{
    if (u == NULL) {
	return E_UNKVAR;
    }

    if (M != u->M) {
#if MDEBUG
	fprintf(stderr, " freeing u->M at %p, replacing with "
		"matrix at %p\n", u->M, M);
#endif
	usermat_destroy_matrix(u);
	u->M = M;
    }

    return 0;
}

int user_matrix_get_level (user_matrix *u)
{
    return (u == NULL)? -1 : u->level;
}

int user_matrix_adjust_level (user_matrix *u, int adj)
{
    if (u == NULL) {
	return E_UNKVAR;
    }

    u->level += adj;

#if MDEBUG
    fprintf(stderr, " user matrix at %p, new level = %d\n", (void *) u,
	    u->level);
#endif

    return 0;
}

int user_matrix_set_name (user_matrix *u, const char *name)
{
    if (u == NULL) {
	return E_UNKVAR;
    }

    *u->name = '\0';
    strncat(u->name, name, VNAMELEN - 1);

#if MDEBUG
    fprintf(stderr, " user matrix at %p, new name = '%s'\n", (void *) u,
	    name);
#endif

    return 0;
}

const char *user_matrix_get_name (user_matrix *u)
{
    return (u == NULL)? NULL : u->name;
}

int user_matrix_replace_matrix_by_name (const char *name, 
					gretl_matrix *M)
{
    user_matrix *u = get_user_matrix_by_name(name);  

    return user_matrix_replace_matrix(u, M);
}

user_matrix *get_user_matrix_by_data (const gretl_matrix *M)
{
    int level = gretl_function_depth();
    int i;

    for (i=0; i<n_matrices; i++) {
	if (matrices[i]->M == M && matrices[i]->level == level) {
	    return matrices[i];
	}
    }

    return NULL;
}

/**
 * get_matrix_by_name:
 * @name: name of the matrix.
 *
 * Looks up a user-defined matrix by name.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

gretl_matrix *get_matrix_by_name (const char *name)
{
    if (name == NULL || *name == '\0') {
	return NULL;
    } else {
	return real_get_matrix_by_name(name, LEVEL_AUTO);
    }
}

/**
 * steal_matrix_by_name:
 * @name: name of the matrix.
 *
 * Looks up a user-defined matrix by name and if found,
 * grabs the matrix, leaving the matrix pointer on the
 * named matrix as %NULL.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

gretl_matrix *steal_matrix_by_name (const char *name)
{
    gretl_matrix *ret = NULL;

    if (name != NULL && *name != '\0') {
	user_matrix *u = get_user_matrix_by_name(name);

	if (u != NULL) {
	    ret = u->M;
	    u->M = NULL;
	}
    }

    return ret;
}

/**
 * get_matrix_copy_by_name:
 * @name: name of the matrix.
 * @err: location to receive error code;
 *
 * Looks up a user-defined matrix by name.
 *
 * Returns: a copy of the named matrix, or %NULL if not found.
 */

gretl_matrix *get_matrix_copy_by_name (const char *name, int *err)
{
    gretl_matrix *m;

    m = real_get_matrix_by_name(name, LEVEL_AUTO);

    if (m == NULL) {
	*err = E_UNKVAR;
    } else {
	m = gretl_matrix_copy(m);
	if (m == NULL) {
	    *err = E_ALLOC;
	}
    }

    return m;
}

/**
 * get_matrix_by_name_at_level:
 * @name: name of the matrix.
 * @level: level of function execution at which to search.
 *
 * Looks up a user-defined matrix by @name, at the given
 * @level of function execution.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

gretl_matrix *get_matrix_by_name_at_level (const char *name, int level)
{
    return real_get_matrix_by_name(name, level);
}

/**
 * user_matrix_get_matrix:
 * @u: user-matrix pointer.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

gretl_matrix *user_matrix_get_matrix (user_matrix *u)
{
    int i;

    for (i=0; i<n_matrices; i++) {
	if (matrices[i] == u) {
	    return matrices[i]->M;
	}
    }

    return NULL;
}

/**
 * copy_named_matrix_as:
 * @orig: the name of the original matrix.
 * @newname: the name to be given to the copy.
 *
 * If a saved matrix is found by the name @orig, a copy of
 * this matrix is added to the stack of saved matrices under the
 * name @newname.  This is intended for use when a matrix is given
 * as the argument to a user-defined function: it is copied
 * under the name assigned by the function's parameter list.
 *
 * Returns: 0 on success, non-zero on error.
 */

int copy_named_matrix_as (const char *orig, const char *newname)
{
    user_matrix *u;
    int err = 0;

    u = get_user_matrix_by_name(orig);
    if (u == NULL) {
	err = 1;
    } else {
	gretl_matrix *M = gretl_matrix_copy(u->M);

	if (M == NULL) {
	    err = E_ALLOC;
	} else {
	    err = user_matrix_add(M, newname);
	}
	if (!err) {
	    /* for use in functions: increment level of last-added matrix */
	    u = matrices[n_matrices - 1];
	    u->level += 1;
	}
    }

    return err;
}

/**
 * copy_matrix_as:
 * @m: the original matrix.
 * @newname: the name to be given to the copy.
 *
 * A copy of matrix @m is added to the stack of saved matrices
 * under the name @newname.  This is intended for use when a matrix
 * is given as the argument to a user-defined function: it is copied
 * under the name assigned by the function's parameter list.
 *
 * Returns: 0 on success, non-zero on error.
 */

int copy_matrix_as (const gretl_matrix *m, const char *newname)
{
    gretl_matrix *m2 = gretl_matrix_copy(m);
    user_matrix *u;
    int err = 0;

    if (m2 == NULL) {
	err = E_ALLOC;
    } else {
	u = real_user_matrix_add(m2, newname, 0);
	if (u == NULL) {
	    err = E_ALLOC;
	} else {
	    u->level += 1;
	}
    }

    return err;
}

/**
 * matrix_add_as_shell:
 * @M: the matrix to add.
 * @name: the name to be given to the "shell".
 *
 * Matrix @M is added to the stack of saved matrices under
 * the name @name with the shell flag set.  This is used
 * when an anonymous matrix is given as a %const argument to a 
 * user-defined function: it is temporarily given user_matrix, 
 * status, so that it is accessible by name within the function,
 * but the content @M is protected from destruction on exit 
 * from the function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int matrix_add_as_shell (gretl_matrix *M, const char *name)
{
    user_matrix *u = real_user_matrix_add(M, name, 0);
    int err = 0;

    if (u == NULL) {
	err = E_ALLOC;
    } else {
	u->flags |= UM_SHELL;
	u->level += 1;
    }

    return err;
}

static int msel_out_of_bounds (int *range, int n)
{
    int *bad = NULL;

    if (range[0] < 1 || range[0] > n) {
	bad = &range[0];
    } else if (range[1] < 1 || range[1] > n) {
	bad = &range[1];
    }

    if (bad != NULL) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), 
			     *bad);
	return 1;
    } else {
	return 0;
    }
}

static int vec_is_const (const gretl_matrix *m, int n)
{
    int i;

    for (i=1; i<n; i++) {
	if (m->val[i] != m->val[i-1]) {
	    return 0;
	}
    }

    return 1;
}

/* convert a matrix subspec component into list of rows 
   or columns */

static int *mspec_make_list (int type, union msel *sel, int n,
			     int *err)
{
    int *slice = NULL;
    int i, ns = 0;

    if (type == SEL_ALL || type == SEL_NULL) {
	return NULL;
    }

    if (type == SEL_MATRIX) {
	if (sel->m == NULL) {
	    gretl_errmsg_set(_("Range is non-positive!"));
	    *err = E_DATA;
	} else {
	    ns = gretl_vector_get_length(sel->m);
	    if (ns == 0) {
		fprintf(stderr, "selection matrix is %d x %d\n", 
			sel->m->rows, sel->m->cols);
		*err = E_NONCONF;
	    } else if (vec_is_const(sel->m, ns)) {
		ns = 1;
	    }
	}
    } else {
	/* range or element */
	if (sel->range[1] == MSEL_MAX) {
	    sel->range[1] = n;
	}
	if (msel_out_of_bounds(sel->range, n)) {
	    *err = E_DATA;
	} else {
	    ns = sel->range[1] - sel->range[0] + 1;
	    if (ns <= 0) {
		gretl_errmsg_sprintf(_("Range %d to %d is non-positive!"),
				     sel->range[0], sel->range[1]); 
		*err = E_DATA;
	    }
	}
    } 

    if (*err) {
	return NULL;
    }

    slice = gretl_list_new(ns);
    if (slice == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<slice[0]; i++) {
	if (type == SEL_MATRIX) {
	    slice[i+1] = sel->m->val[i];
	} else {
	    slice[i+1] = sel->range[0] + i;
	} 
    }

    for (i=1; i<=slice[0] && !*err; i++) {
	if (slice[i] < 1 || slice[i] > n) {
	    gretl_errmsg_sprintf(_("Index value %d is out of bounds"), 
				 slice[i]);
	    *err = 1;
	}
    }

    if (*err) {
	free(slice);
	slice = NULL;
    }

    return slice;
}

/* catch the case of an implicit column or row specification for
   a sub-matrix of an (n x 1) or (1 x m) matrix; also catch the
   error of giving just one row/col spec for a matrix that has
   more than one row and more than one column
*/

int check_matrix_subspec (matrix_subspec *spec, const gretl_matrix *m)
{
    int err = 0;

    if (spec->type[1] == SEL_NULL) {
	/* we got only one row/col spec */
	if (m->cols == 1) {
	    /* OK: implicitly col = 1 */
	    spec->type[1] = SEL_RANGE;
	    mspec_set_col_index(spec, 1);
	} else if (m->rows == 1) {
	    /* OK: implicitly row = 1, and transfer the single 
	       given spec to the column dimension */
	    spec->type[1] = spec->type[0];
	    if (spec->type[1] == SEL_MATRIX) {
		spec->sel[1].m = spec->sel[0].m;
	    } else {
		spec->sel[1].range[0] = spec->sel[0].range[0];
		spec->sel[1].range[1] = spec->sel[0].range[1];
	    }		
	    spec->type[0] = SEL_RANGE;
	    mspec_set_row_index(spec, 1);
	} else {
	    gretl_errmsg_set(_("Ambiguous matrix index"));
	    err = E_DATA;
	}	    
    }

    if (spec->type[0] == SEL_RANGE && spec->type[1] == SEL_RANGE) {
	if (spec->sel[0].range[0] == spec->sel[0].range[1] &&
	    spec->sel[1].range[0] == spec->sel[1].range[1]) {
	    spec->type[0] = spec->type[1] = SEL_ELEMENT;
	}
    }

    return err;
}

static int get_slices (matrix_subspec *spec, 
		       const gretl_matrix *M)
{
    int err = 0;

    spec->rslice = mspec_make_list(spec->type[0], &spec->sel[0], 
				   M->rows, &err);

    if (!err) {
	spec->cslice = mspec_make_list(spec->type[1], &spec->sel[1], 
				       M->cols, &err);
    }

    return err;
}

int assign_scalar_to_submatrix (gretl_matrix *M, double x,
				matrix_subspec *spec)
{
    int mr = gretl_matrix_rows(M);
    int mc = gretl_matrix_cols(M);
    int i, err = 0;

    if (spec == NULL) {
	fprintf(stderr, "matrix_replace_submatrix: spec is NULL!\n");
	return E_DATA;
    }

    if (spec->type[0] == SEL_DIAG) {
	int n = (mr < mc)? mr : mc;

	for (i=0; i<n; i++) {
	    gretl_matrix_set(M, i, i, x);
	}
	return 0;
    }

    if (spec->rslice == NULL && spec->cslice == NULL) {
	/* parse mspec into lists of affected rows and columns */
	err = get_slices(spec, M);
    }

    if (!err) {
	int sr = (spec->rslice == NULL)? mr : spec->rslice[0];
	int sc = (spec->cslice == NULL)? mc : spec->cslice[0];
	int j, l, k = 0;
	int mi, mj;

	for (i=0; i<sr; i++) {
	    mi = (spec->rslice == NULL)? k++ : spec->rslice[i+1] - 1;
	    l = 0;
	    for (j=0; j<sc; j++) {
		mj = (spec->cslice == NULL)? l++ : spec->cslice[j+1] - 1;
		gretl_matrix_set(M, mi, mj, x);
	    }
	}
    }

    return err;
}

/* @M is the target for partial replacement, @S is the source to
   substitute, and @spec tells how/where to make the
   substitution.
*/

static int matrix_replace_submatrix (gretl_matrix *M,
				     const gretl_matrix *S,
				     matrix_subspec *spec)
{
    int mr = gretl_matrix_rows(M);
    int mc = gretl_matrix_cols(M);
    int sr = gretl_matrix_rows(S);
    int sc = gretl_matrix_cols(S);
    int sscalar = 0;
    int err = 0;

    if (spec == NULL) {
	fprintf(stderr, "matrix_replace_submatrix: spec is NULL!\n");
	return E_DATA;
    }

    if (sr > mr || sc > mc) {
	/* the replacement matrix won't fit into M */
	fprintf(stderr, "matrix_replace_submatrix: target is %d x %d but "
		"replacement part is %d x %d\n", mr, mc, sr, sc);
	return E_NONCONF;
    }

    if (spec->type[0] == SEL_DIAG) {
	return matrix_insert_diagonal(M, S, mr, mc);
    }

    if (spec->rslice == NULL && spec->cslice == NULL) {
	/* parse mspec into lists of affected rows and columns */
	err = get_slices(spec, M);
	if (err) {
	    return err;
	}
    }

#if MDEBUG
    printlist(spec->rslice, "rslice (rows list)");
    printlist(spec->cslice, "cslice (cols list)");
    fprintf(stderr, "orig M = %d x %d, S = %d x %d\n", mr, mc, sr, sc);
#endif

    if (sr == 1 && sc == 1) {
	/* selection matrix is a scalar */
	sscalar = 1;
	sr = (spec->rslice == NULL)? mr : spec->rslice[0];
	sc = (spec->cslice == NULL)? mc : spec->cslice[0];
    } else if (spec->rslice != NULL && spec->rslice[0] != sr) {
	fprintf(stderr, "mspec has %d rows but substitute matrix has %d\n", 
		spec->rslice[0], sr);
	err = E_NONCONF;
    } else if (spec->cslice != NULL && spec->cslice[0] != sc) {
	fprintf(stderr, "mspec has %d cols but substitute matrix has %d\n", 
		spec->cslice[0], sc);
	err = E_NONCONF;
    }

    if (!err) {
	int i, j, l, k = 0;
	int mi, mj;
	double x;

	x = (sscalar)? S->val[0] : 0.0;

	for (i=0; i<sr; i++) {
	    mi = (spec->rslice == NULL)? k++ : spec->rslice[i+1] - 1;
	    l = 0;
	    for (j=0; j<sc; j++) {
		mj = (spec->cslice == NULL)? l++ : spec->cslice[j+1] - 1;
		if (!sscalar) {
		    x = gretl_matrix_get(S, i, j);
		}
		gretl_matrix_set(M, mi, mj, x);
	    }
	}
    }

    return err;
}

gretl_matrix *matrix_get_submatrix (const gretl_matrix *M, 
				    matrix_subspec *spec,
				    int prechecked, 
				    int *err)
{
    gretl_matrix *S;
    int r, c;

    if (gretl_is_null_matrix(M)) {
	*err = E_DATA;
	return NULL;
    }

    if (!prechecked) {
	*err = check_matrix_subspec(spec, M);
	if (*err) {
	    return NULL;
	}
    }

    if (spec->type[0] == SEL_DIAG) {
	return gretl_matrix_get_diagonal(M, err);
    } else if (spec->type[0] == SEL_ELEMENT) {
	int i = mspec_get_row_index(spec);
	int j = mspec_get_col_index(spec);
	double x = matrix_get_element(M, i, j, err);

	return (*err)? NULL : gretl_matrix_from_scalar(x);
    }

    if (spec->rslice == NULL && spec->cslice == NULL) {
	*err = get_slices(spec, M);
	if (*err) {
	    return NULL;
	}
    }

#if MDEBUG
    printlist(spec->rslice, "rslice");
    printlist(spec->cslice, "cslice");
    fprintf(stderr, "M = %d x %d\n", M->rows, M->cols);
#endif

    r = (spec->rslice == NULL)? M->rows : spec->rslice[0];
    c = (spec->cslice == NULL)? M->cols : spec->cslice[0];

    S = gretl_matrix_alloc(r, c);
    if (S == NULL) {
	*err = E_ALLOC;	
    }

    if (S != NULL) {
	int i, j, k, l;
	int mi, mj;
	double x;

	k = 0;
	for (i=0; i<r && !*err; i++) {
	    mi = (spec->rslice == NULL)? k++ : spec->rslice[i+1] - 1;
	    l = 0;
	    for (j=0; j<c && !*err; j++) {
		mj = (spec->cslice == NULL)? l++ : spec->cslice[j+1] - 1;
		x = gretl_matrix_get(M, mi, mj);
		gretl_matrix_set(S, i, j, x);
	    }
	}
    }

    if (S != NULL && S->rows == M->rows && gretl_matrix_is_dated(M)) {
	int mt1 = gretl_matrix_get_t1(M);
	int mt2 = gretl_matrix_get_t2(M);

	gretl_matrix_set_t1(S, mt1);
	gretl_matrix_set_t2(S, mt2);
    }

    return S;
}

double matrix_get_element (const gretl_matrix *M, int i, int j,
			   int *err)
{
    double x = NADBL;

    /* The incoming i and j are from userspace, and will
       be 1-based.
    */
    i--;
    j--;

    if (M == NULL) {
	*err = E_DATA;
    } else if (i < 0 || i >= M->rows || j < 0 || j >= M->cols) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), 
			     (i < 0 || i >= M->rows)? (j+1) : (i+1));
	*err = 1;
    } else {
	x = gretl_matrix_get(M, i, j);
    }

    return x;
}

gretl_matrix *user_matrix_get_submatrix (const char *name, 
					 matrix_subspec *spec,
					 int *err)
{
    gretl_matrix *S = NULL;
    gretl_matrix *M;

    M = real_get_matrix_by_name(name, LEVEL_AUTO); 
    if (M == NULL) {
	*err = E_UNKVAR;
    } else {
	S = matrix_get_submatrix(M, spec, 0, err);
    }

    return S;
}

/* Look up the existing matrix called @name, and substitute
   the matrix @S for part of the original, as specified by
   @spec.
*/

int user_matrix_replace_submatrix (const char *mname, 
				   const gretl_matrix *S,
				   matrix_subspec *spec)
{
    gretl_matrix *M;

    M = real_get_matrix_by_name(mname, LEVEL_AUTO); 
    if (M == NULL) {
	return E_UNKVAR;
    }

    return matrix_replace_submatrix(M, S, spec);
}

/**
 * add_or_replace_user_matrix:
 * @M: gretl matrix.
 * @name: name for the matrix.
 *
 * Checks whether a matrix of the given @name already exists.
 * If so, the original matrix is replaced by @M; if not, the
 * the matrix @M is added to the stack of user-defined
 * matrices.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int add_or_replace_user_matrix (gretl_matrix *M, const char *name)
{
    int err = 0;

    if (get_user_matrix_by_name(name) != NULL) {
	err = user_matrix_replace_matrix_by_name(name, M);
    } else {
	err = user_matrix_add(M, name);
    }

    return err;
}

static void destroy_user_matrix (user_matrix *u)
{
    if (u == NULL) {
	return;
    }

    if (matrix_is_shell(u)) {
	/* don't free content */
	free(u);
	return;
    }

#if MDEBUG
    fprintf(stderr, "destroy_user_matrix: freeing matrix at %p...", 
	    (void *) u->M);
#endif

    usermat_destroy_matrix(u);
    free(u);

#if MDEBUG
    fprintf(stderr, " done\n");
#endif
}

#define LEV_PRIVATE -1

static int matrix_levels_match (user_matrix *u, int lev)
{
    int ret = 0;

    if (u->level == lev) {
	ret = 1;
    } else if (lev == LEV_PRIVATE && matrix_is_private(u)) {
	ret = 1;
    }

    return ret;
}

/**
 * destroy_user_matrices_at_level:
 * @level: stack level of function execution.
 *
 * Destroys and removes from the stack of user matrices all
 * matrices that were created at the given @level.  This is 
 * part of the cleanup that is performed when a user-defined
 * function terminates.
 *
 * Returns: 0 on success, non-zero on error.
 */

int destroy_user_matrices_at_level (int level)
{
    user_matrix **tmp;
    int i, j, nm = 0;
    int err = 0;

#if MDEBUG
    fprintf(stderr, "destroy_user_matrices_at_level: level = %d, "
	    "total n_matrices = %d\n", level, n_matrices);
#endif
#if MDEBUG > 1
    print_matrix_stack("at top of destroy_user_matrices_at_level");
#endif

    for (i=0; i<n_matrices; i++) {
	if (matrices[i] == NULL) {
	    break;
	}
	if (matrix_levels_match(matrices[i], level)) {
#if MDEBUG
	    fprintf(stderr, "destroying matrix[%d] ('%s' at %p)\n",
		    i, matrices[i]->name, (void *) matrices[i]->M);
#endif
	    destroy_user_matrix(matrices[i]);
	    for (j=i; j<n_matrices - 1; j++) {
		matrices[j] = matrices[j+1];
	    }
	    matrices[n_matrices - 1] = NULL;
	    i--;
	} else {
	    nm++;
	}
    }

    if (nm < n_matrices) {
	n_matrices = nm;
	if (nm == 0) {
	    free(matrices);
	    matrices = NULL;
	} else {
	    tmp = realloc(matrices, nm * sizeof *tmp);
	    if (tmp == NULL) {
		err = E_ALLOC;
	    } else {
		matrices = tmp;
	    }
	}
    }

#if MDEBUG > 1
    print_matrix_stack("at end of destroy_user_matrices_at_level");
#endif

    return err;
}

int destroy_private_matrices (void)
{
    return destroy_user_matrices_at_level(LEV_PRIVATE);    
}

/**
 * destroy_user_matrices:
 *
 * Frees all resources associated with the stack of user-
 * defined matrices.
 */

void destroy_user_matrices (void)
{
    int i;

#if MDEBUG
    fprintf(stderr, "destroy_user_matrices called, n_matrices = %d\n",
	    n_matrices);
#endif

    if (matrices == NULL) {
	return;
    }

    for (i=0; i<n_matrices; i++) {
#if MDEBUG
	fprintf(stderr, "destroying user_matrix %d (%s) at %p\n", i,
		matrices[i]->name, (void *) matrices[i]);
#endif
	destroy_user_matrix(matrices[i]);
    }

    free(matrices);
    matrices = NULL;
    n_matrices = 0;
}

int user_matrix_destroy (user_matrix *u)
{
    int err = 0;

    if (u == NULL) {
	err = E_UNKVAR;
    } else {
	int i, j, nm = n_matrices - 1;

	for (i=0; i<n_matrices; i++) {
	    if (matrices[i] == u) {
		destroy_user_matrix(matrices[i]);
		for (j=i; j<nm; j++) {
		    matrices[j] = matrices[j+1];
		}
		matrices[nm] = NULL;
		break;
	    }
	} 

	if (nm == 0) {
	    free(matrices);
	    matrices = NULL;
	} else {
	    user_matrix **tmp = realloc(matrices, nm * sizeof *tmp);

	    if (tmp == NULL) {
		err = E_ALLOC;
	    } else {
		matrices = tmp;
	    }
	}

	n_matrices = nm;
    }

    return err;
}

int user_matrix_destroy_by_name (const char *name, PRN *prn)
{
    user_matrix *u = get_user_matrix_by_name(name);
    int err;

    err = user_matrix_destroy(u);

    if (!err && prn != NULL && gretl_messages_on()) {
	pprintf(prn, _("Deleted matrix %s"), name);
	pputc(prn, '\n');
    }

    if (matrix_add_delete_callback != NULL && 
	gretl_function_depth() == 0) {
	(*matrix_add_delete_callback)(name, USER_MATRIX_DELETE);
    }    

    return err;
}

int umatrix_set_names_from_string (gretl_matrix *M, 
				   const char *s,
				   int byrow)
{
    user_matrix *u = get_user_matrix_by_data(M);
    int n, err = 0;

    if (u == NULL) {
	return E_UNKVAR;
    }

    n = (byrow)? M->rows : M->cols;

    if (s == NULL || *s == '\0') {
	if (byrow) {
	    gretl_matrix_set_rownames(M, NULL);
	} else {
	    gretl_matrix_set_colnames(M, NULL);
	}
    } else {
	char **S;
	int ns;

	S = gretl_string_split(s, &ns);
	if (S == NULL) {
	    err = E_ALLOC;
	} else if (ns != n) {
	    err = E_NONCONF;
	    free_strings_array(S, ns);
	} else if (byrow) {
	    gretl_matrix_set_rownames(M, S);
	} else {
	    gretl_matrix_set_colnames(M, S);
	} 
    }

    return err;
}

int umatrix_set_names_from_list (gretl_matrix *M, 
				 const int *list,
				 const DATASET *dset,
				 int byrow)
{
    user_matrix *u = get_user_matrix_by_data(M);
    int i, n, err = 0;

    if (u == NULL) {
	return E_UNKVAR;
    }

    n = (byrow)? M->rows : M->cols;

    if (list == NULL || list[0] == 0) {
	if (byrow) {
	    gretl_matrix_set_rownames(M, NULL);
	} else {
	    gretl_matrix_set_colnames(M, NULL);
	}
    } else if (list[0] != n) {
	err = E_NONCONF;
    } else {
	char **S = strings_array_new(n);

	if (S == NULL) {
	    err = E_ALLOC;
	}

	for (i=0; i<n && !err; i++) {
	    S[i] = gretl_strndup(dset->varname[list[i+1]], 12);
	    if (S[i] == NULL) {
		err = E_ALLOC;
	    }
	}

	if (err) {
	    free_strings_array(S, n);
	} else if (byrow) {
	    gretl_matrix_set_rownames(M, S);
	} else {
	    gretl_matrix_set_colnames(M, S);
	}
    }

    return err;
}

char *user_matrix_get_column_name (const gretl_matrix *M, int col,
				   int *err)
{
    char *ret = NULL;

    if (M == NULL || col < 1 || col > M->cols) {
	*err = E_DATA;
    } else {
	const char **S = gretl_matrix_get_colnames(M);

	if (S == NULL) {
	    ret = gretl_strdup("");
	} else {
	    ret = gretl_strdup(S[col-1]);
	}
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

double 
user_matrix_get_determinant (gretl_matrix *m, int f, int *err)
{
    gretl_matrix *tmp = NULL;
    double d = NADBL;

    if (gretl_is_null_matrix(m)) {
	return d;
    } else if (matrix_is_saved(m)) {
	tmp = gretl_matrix_copy(m);
    } else {
	tmp = m;
    }

    if (tmp != NULL) {
	if (f == F_LDET) {
	    d = gretl_matrix_log_determinant(tmp, err);
	} else {
	    d = gretl_matrix_determinant(tmp, err);
	}
	if (tmp != m) {
	    gretl_matrix_free(tmp);
	}
    }

    return d;
}

gretl_matrix *user_matrix_matrix_func (gretl_matrix *m, int f, 
				       int *err)
{
    gretl_matrix *R = NULL;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
    } else if (matrix_is_saved(m)) {
	/* don't mess with the original matrix! */
	R = gretl_matrix_copy(m);
	if (R == NULL) {
	    *err = E_ALLOC;
	}	
    } else {
	/* use (overwrite) the input matrix */
	R = m;
    } 

    if (R != NULL) {
	if (f == F_CDEMEAN) {
	    gretl_matrix_demean_by_column(R);
	} else if (f == F_CHOL) {
	    *err = gretl_matrix_cholesky_decomp(R);
	} else if (f == F_PSDROOT) {
	    *err = gretl_matrix_psd_root(R);
	} else if (f == F_INVPD) {
	    *err = gretl_invpd(R);
	} else if (f == F_GINV) {
	    *err = gretl_matrix_moore_penrose(R);
	} else if (f == F_INV) {
	    *err = gretl_invert_matrix(R);
	} else if (f == F_UPPER) {
	    *err = gretl_matrix_zero_lower(R);
	} else if (f == F_LOWER) {
	    *err = gretl_matrix_zero_upper(R);
	} else {
	    *err = E_DATA;
	}
	if (*err && R != m) {
	    gretl_matrix_free(R);
	    R = NULL;
	}
    } 
   
    return R;
}

static void matrix_cannibalize (gretl_matrix *targ, gretl_matrix *src)
{
    targ->rows = src->rows;
    targ->cols = src->cols;

    free(targ->val);
    targ->val = src->val;
    src->val = NULL;
}

int matrix_invert_in_place (gretl_matrix *m)
{
    gretl_matrix *R = gretl_matrix_copy(m);
    int err = 0;

    if (R == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_invert_matrix(R);
	if (!err) {
	    matrix_cannibalize(m, R);
	}
	gretl_matrix_free(R);
    } 

    return err;
}

int matrix_cholesky_in_place (gretl_matrix *m)
{
    gretl_matrix *R = gretl_matrix_copy(m);
    int err = 0;

    if (R == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_matrix_cholesky_decomp(R);
	if (!err) {
	    matrix_cannibalize(m, R);
	}
	gretl_matrix_free(R);
    } 

    return err;
}

int matrix_transpose_in_place (gretl_matrix *m)
{
    gretl_matrix *R = gretl_matrix_copy_transpose(m);
    int err = 0;

    if (R == NULL) {
	err = E_ALLOC;
    } else {
	matrix_cannibalize(m, R);
	gretl_matrix_free(R);
    }

    return err;
}

int matrix_XTX_in_place (gretl_matrix *m)
{
    gretl_matrix *R = gretl_matrix_alloc(m->cols, m->cols);
    int err;

    if (R == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_matrix_multiply_mod(m, GRETL_MOD_TRANSPOSE,
					m, GRETL_MOD_NONE,
					R, GRETL_MOD_NONE);
    }

    if (!err) {
	matrix_cannibalize(m, R);
    }

    gretl_matrix_free(R);

    return err;
}

gretl_matrix *user_matrix_vec (const gretl_matrix *m, int *err)
{
    gretl_matrix *R = NULL;

    if (gretl_is_null_matrix(m)) {
	R = gretl_null_matrix_new();
    } else {
	R = gretl_matrix_alloc(m->rows * m->cols, 1);
	if (R != NULL) {
	    gretl_matrix_vectorize(R, m);
	} 
    }

    if (R == NULL) {
	*err = E_ALLOC;
    }

    return R;
}

gretl_matrix *user_matrix_vech (const gretl_matrix *m, int *err)
{
    gretl_matrix *R = NULL;

    if (gretl_is_null_matrix(m)) {
	R = gretl_null_matrix_new();
    } else if (m->rows != m->cols) {
	*err = E_NONCONF;
    } else {
	int n = m->rows;
	int k = n * (n + 1) / 2;

	R = gretl_matrix_alloc(k, 1);
	if (R != NULL) {
	    *err = gretl_matrix_vectorize_h(R, m);
	}
    } 

    if (R == NULL && !*err) {
	*err = E_ALLOC;
    }

    return R;
}

gretl_matrix *user_matrix_unvech (const gretl_matrix *m, int *err)
{
    gretl_matrix *R = NULL;

    if (gretl_is_null_matrix(m)) {
	R = gretl_null_matrix_new();
    } else if (m->cols != 1) {
	*err = E_NONCONF;
    } else {
	int n = (int) ((sqrt(1.0 + 8.0 * m->rows) - 1.0) / 2.0);

	R = gretl_matrix_alloc(n, n);
	if (R != NULL) {
	    *err = gretl_matrix_unvectorize_h(R, m);
	} 
    }

    if (R == NULL && !*err) {
	*err = E_ALLOC;
    }

    return R;
}

static int 
real_user_matrix_QR_decomp (const gretl_matrix *m, gretl_matrix **Q, 
			    gretl_matrix **R)
{
    int mc = gretl_matrix_cols(m);
    int err = 0;

    *Q = gretl_matrix_copy(m);

    if (*Q == NULL) {
	err = E_ALLOC;
    } else if (R != NULL) {
	*R = gretl_matrix_alloc(mc, mc);
	if (*R == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = gretl_matrix_QR_decomp(*Q, (R == NULL)? NULL : *R);
    }

    if (err) {
	gretl_errmsg_set(_("Matrix decomposition failed"));
	gretl_matrix_free(*Q);
	*Q = NULL;
	if (R != NULL) {
	    gretl_matrix_free(*R);
	    *R = NULL;
	}
    }

    return err;
}

#define nullarg(s) (s == NULL || !strcmp(s, "null"))

gretl_matrix *
user_matrix_QR_decomp (const gretl_matrix *m, const char *rname, int *err)
{
    gretl_matrix *Q = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix **pR = NULL;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (!nullarg(rname)) {
	if (get_matrix_by_name(rname) == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix"), rname);
	    *err = E_UNKVAR;
	} else {
	    pR = &R;
	}
    }

    if (!*err) {
	*err = real_user_matrix_QR_decomp(m, &Q, pR);
    }

    if (!*err && R != NULL) {
	user_matrix_replace_matrix_by_name(rname, R);
    }

    return Q;
}

static int revise_SVD_V (gretl_matrix **pV, int r, int c)
{
    gretl_matrix *V;
    double x;
    int i, j;

    V = gretl_matrix_alloc(r, c);
    if (V == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<r; i++) {
	for (j=0; j<c; j++) {
	    x = gretl_matrix_get((*pV), i, j);
	    gretl_matrix_set(V, i, j, x);
	}
    }

    gretl_matrix_free(*pV);
    *pV = V;

    return 0;
}

gretl_matrix *user_matrix_SVD (const gretl_matrix *m, 
			       const char *uname, 
			       const char *vname, 
			       int *err)
{
    gretl_matrix *U = NULL;
    gretl_matrix *S = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix **pU = NULL;
    gretl_matrix **pV = NULL;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (!nullarg(uname)) {
	if (get_matrix_by_name(uname) == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix"), uname);
	    *err = E_UNKVAR;
	} else {
	    pU = &U;
	}
    }

    if (!*err && !nullarg(vname)) {
	if (get_matrix_by_name(vname) == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix"), vname);
	    *err = E_UNKVAR;
	} else {
	    pV = &V;
	}
    }

    if (!*err) {
	*err = gretl_matrix_SVD(m, pU, &S, pV);
    }

    if (!*err && (U != NULL || V != NULL)) {
	int tall = m->rows - m->cols;
	int minrc = (m->rows > m->cols)? m->cols : m->rows;

	if (U != NULL) {
	    if (tall > 0) {
		*err = gretl_matrix_realloc(U, m->rows, minrc);
	    }
	    if (!*err) {
		user_matrix_replace_matrix_by_name(uname, U);
	    }
	}
	if (V != NULL) {
	    if (tall < 0) {
		*err = revise_SVD_V(&V, minrc, m->cols);
	    } 
	    if (!*err) {
		user_matrix_replace_matrix_by_name(vname, V);
	    }
	}
    }

    return S;
}

static gretl_matrix *get_ols_matrix (const char *mname, 
				     int r, int c,
				     int *newmat, int *err)
{
    gretl_matrix *m = get_matrix_by_name(mname);

    if (m == NULL) {
	gretl_errmsg_sprintf(_("'%s': no such matrix"), mname);
	*err = E_UNKVAR;
    } else if (newmat != NULL && (m->rows != r || m->cols != c)) {
	m = gretl_matrix_alloc(r, c);
	if (m == NULL) {
	    *err = E_ALLOC;
	} else {
	    *newmat = 1;
	}
    }

    return m;
}

gretl_matrix *user_matrix_ols (const gretl_matrix *Y, 
			       const gretl_matrix *X, 
			       const char *Uname, 
			       const char *Vname, 
			       gretlopt opt,
			       int *err)
{
    gretl_matrix *B = NULL;
    gretl_matrix *U = NULL;
    gretl_matrix *V = NULL;
    int newU = 0, newV = 0;
    double s2, *ps2 = NULL;
    int g, k, T;

    if (gretl_is_null_matrix(Y) || gretl_is_null_matrix(X)) {
	*err = E_DATA;
	return NULL;
    }

    T = Y->rows;
    k = X->cols;
    g = Y->cols;

    if (X->rows != T) {
	*err = E_NONCONF;
	return NULL;
    }

    if (g > 1 && (opt & OPT_M)) {
	/* multiple precision: only one y var wanted */
	*err = E_DATA;
	return NULL;
    }

    if (!nullarg(Uname)) {
	U = get_ols_matrix(Uname, T, g, &newU, err);
	if (*err) {
	    return NULL;
	}
    } 

    if (!nullarg(Vname)) {
	if (g > 1) {
	    /* multiple dependent variables */
	    get_ols_matrix(Vname, 0, 0, NULL, err);
	    if (!*err) {
		newV = 1;
	    }
	} else {
	    /* single dependent variable */
	    int nv = g * k;

	    V = get_ols_matrix(Vname, nv, nv, &newV, err);
	    if (!*err) {
		ps2 = &s2;
	    }
	}
    }

    if (!*err) {
	B = gretl_matrix_alloc(k, g);
	if (B == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	if (g == 1) {
	    if (opt & OPT_M) {
		/* use multiple precision */
		*err = gretl_matrix_mp_ols(Y, X, B, V, U, ps2);
	    } else {
		*err = gretl_matrix_ols(Y, X, B, V, U, ps2);
	    }
	} else {
	    if (newV) {
		/* note: "V" will actually be (X'X)^{-1} */
		*err = gretl_matrix_multi_ols(Y, X, B, U, &V);
	    } else {
		*err = gretl_matrix_multi_ols(Y, X, B, U, NULL);
	    }
	}
    }

    if (*err) {
	gretl_matrix_free(B);
	B = NULL;
	if (newU) gretl_matrix_free(U);
	if (newV) gretl_matrix_free(V);
    } else {
	if (newU) {
	    user_matrix_replace_matrix_by_name(Uname, U);
	}
	if (newV) {
	    user_matrix_replace_matrix_by_name(Vname, V);
	}	
    }

    return B;
}

gretl_matrix *user_matrix_rls (const gretl_matrix *Y, 
			       const gretl_matrix *X,
			       const gretl_matrix *R,
			       const gretl_matrix *Q,
			       const char *Uname, 
			       const char *Vname, 
			       int *err)
{
    gretl_matrix *B = NULL;
    gretl_matrix *U = NULL;
    gretl_matrix *V = NULL;
    int newU = 0, newV = 0;
    int g, k, T;

    if (gretl_is_null_matrix(Y) || gretl_is_null_matrix(X)) {
	*err = E_DATA;
	return NULL;
    }

    T = Y->rows;
    k = X->cols;
    g = Y->cols;

    if (X->rows != T) {
	*err = E_NONCONF;
	return NULL;
    }

    if (!nullarg(Uname)) {
	U = get_ols_matrix(Uname, T, g, &newU, err);
	if (*err) {
	    return NULL;
	}
    } 

    if (!nullarg(Vname)) {
	get_ols_matrix(Vname, 0, 0, NULL, err);
	if (!*err) {
	    newV = 1;
	}
    }

    if (!*err) {
	B = gretl_matrix_alloc(k, g);
	if (B == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	if (newV) {
	    /* note: "V" will actually be M (X'X)^{-1} M' */
	    *err = gretl_matrix_restricted_multi_ols(Y, X, R, Q, B, U, &V);
	} else {
	    *err = gretl_matrix_restricted_multi_ols(Y, X, R, Q, B, U, NULL);
	}
    }

    if (*err) {
	gretl_matrix_free(B);
	B = NULL;
	if (newU) gretl_matrix_free(U);
	if (newV) gretl_matrix_free(V);
    } else {
	if (newU) {
	    user_matrix_replace_matrix_by_name(Uname, U);
	}
	if (newV) {
	    user_matrix_replace_matrix_by_name(Vname, V);
	}	
    }

    return B;
}

static void maybe_eigen_trim (gretl_matrix *E)
{
    double x;
    int i, allreal = 1;

    for (i=0; i<E->rows; i++) {
	x = gretl_matrix_get(E, i, 1);
	if (x != 0.0) {
	    allreal = 0;
	    break;
	}
    }

    if (allreal) {
	gretl_matrix_reuse(E, -1, 1);
    }
}

gretl_matrix *
user_matrix_eigen_analysis (const gretl_matrix *m, const char *rname, int symm,
			    int *err)
{
    gretl_matrix *C = NULL;
    gretl_matrix *E = NULL;
    int vecs = 0;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (gretl_matrix_xna_check(m)) {
	*err = E_NAN;
	return NULL;
    }

    if (!nullarg(rname)) {
	vecs = 1;
	if (get_matrix_by_name(rname) == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix"), rname);
	    *err = E_UNKVAR;
	    return NULL;
	}
    }

    C = gretl_matrix_copy(m);
    if (C == NULL) {
	*err = E_ALLOC;
    }

    if (!*err) {
	if (symm) {
	    E = gretl_symmetric_matrix_eigenvals(C, vecs, err);
	} else {
	    E = gretl_general_matrix_eigenvals(C, vecs, err);
	    if (E != NULL && E->cols == 2) {
		maybe_eigen_trim(E);
	    }
	}
    }

    if (!*err && vecs) {
	user_matrix_replace_matrix_by_name(rname, C);
    }

    if (!vecs) {
	gretl_matrix_free(C);
    }

    return E;
}

gretl_matrix *user_gensymm_eigenvals (const gretl_matrix *A, 
				      const gretl_matrix *B,
				      const char *rname,
				      int *err)
{
    gretl_matrix *E = NULL, *V = NULL;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
	*err = E_DATA;
	return NULL;
    }

    if (gretl_matrix_xna_check(A) || gretl_matrix_xna_check(B)) {
	*err = E_NAN;
	return NULL;
    }

    if (!nullarg(rname)) {
	if (get_matrix_by_name(rname) == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix"), rname);
	    *err = E_UNKVAR;
	    return NULL;
	} else {
	    V = gretl_matrix_alloc(B->cols, A->rows);
	    if (V == NULL) {
		*err = E_ALLOC;
		return NULL;
	    }
	}
    }

    E = gretl_gensymm_eigenvals(A, B, V, err);

    if (V != NULL) {
	if (*err) {
	    gretl_matrix_free(V);
	} else {
	    user_matrix_replace_matrix_by_name(rname, V);
	}
    }

    return E;
}

static void xml_put_user_matrix (user_matrix *u, FILE *fp)
{
    gretl_matrix *M;
    const char **S;
    int i, j;

    if (u == NULL || u->M == NULL) {
	return;
    }

    M = u->M;

    fprintf(fp, "<gretl-matrix name=\"%s\" rows=\"%d\" cols=\"%d\"", 
	    u->name, M->rows, M->cols);

    S = gretl_matrix_get_colnames(M);

    if (S != NULL) {
	fputs(" colnames=\"", fp);
	for (j=0; j<M->cols; j++) {
	    fputs(S[j], fp);
	    fputc((j < M->cols - 1)? ' ' : '"', fp);
	}
    } 

    S = gretl_matrix_get_rownames(M);

    if (S != NULL) {
	fputs(" rownames=\"", fp);
	for (j=0; j<M->rows; j++) {
	    fputs(S[j], fp);
	    fputc((j < M->rows - 1)? ' ' : '"', fp);
	}
    }     

    fputs(">\n", fp);

    for (i=0; i<M->rows; i++) {
	for (j=0; j<M->cols; j++) {
	    fprintf(fp, "%.16g ", gretl_matrix_get(M, i, j));
	}
	fputc('\n', fp);
    }

    fputs("</gretl-matrix>\n", fp); 
}

void write_matrices_to_file (FILE *fp)
{
    int i;

    gretl_xml_header(fp);
    fprintf(fp, "<gretl-matrices count=\"%d\">\n", n_matrices);

    gretl_push_c_numeric_locale();

    for (i=0; i<n_matrices; i++) {
	if (matrices[i]->M != NULL) {
	    xml_put_user_matrix(matrices[i], fp);
	}
    }

    gretl_pop_c_numeric_locale();

    fputs("</gretl-matrices>\n", fp);
}
