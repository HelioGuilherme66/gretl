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

/* syntax tree evaluator for "genr" */

#include "genparse.h"
#include "monte_carlo.h"
#include "gretl_string_table.h"
#include "gretl_typemap.h"
#include "matrix_extra.h"
#include "usermat.h"
#include "uservar.h"
#include "gretl_bfgs.h"
#include "gretl_normal.h"
#include "gretl_panel.h"
#include "kalman.h"
#include "libset.h"
#include "version.h"
#include "csvdata.h"
#include "uservar_priv.h"
#include "genr_optim.h"
#include "gretl_cmatrix.h"

#ifdef USE_CURL
# include "gretl_www.h"
#endif

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#ifdef USE_RLIB
# include "gretl_foreign.h"
#endif

#include <errno.h>

#if GENDEBUG
# define EDEBUG GENDEBUG
# define LHDEBUG GENDEBUG
#else
# define EDEBUG 0
# define LHDEBUG 0
#endif

#define AUX_NODES_DEBUG 0

#if AUX_NODES_DEBUG
# include <stdarg.h>
static void real_rndebug (const char *format, ...)
{
    va_list args;

    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
}
# define rndebug(x) real_rndebug x
#else
# define rndebug(x)
#endif

#define MATRIX_NA_CHECK 1
#define ONE_BY_ONE_CAST 1

#define SCALARS_ENSURE_FINITE 1 /* debatable, but watch out for read/write */
#define SERIES_ENSURE_FINITE 1  /* debatable */

enum {
    FR_TREE = 1,
    FR_RET,
    FR_LHTREE,
    FR_LHRES,
    FR_ARET
};

#define is_aux_node(n) (n != NULL && (n->flags & AUX_NODE))
#define is_tmp_node(n) (n != NULL && (n->flags & TMP_NODE))
#define is_proxy_node(n) (n != NULL && (n->flags & PRX_NODE))

#define nullmat_ok(f) (f == F_ROWS || f == F_COLS || f == F_DET || \
		       f == F_LDET || f == F_DIAG || f == F_TRANSP || \
		       f == F_VEC || f == F_VECH || f == F_UNVECH || \
		       f == F_DATAOK)

#define dataset_dum(n) (n->t == DUM && n->v.idnum == DUM_DATASET)

#define postfix_node(n) (n->t == NUM_P || n->t == NUM_M)

#define ok_list_node(n) (n->t == LIST || n->t == NUM || \
			 n->t == MAT || n->t == EMPTY || \
			 (n->t == SERIES && n->vnum >= 0))

#define uscalar_node(n) ((n->t == NUM && n->vname != NULL) || postfix_node(n))

#define umatrix_node(n) (n->t == MAT && n->vname != NULL)
#define ubundle_node(n) (n->t == BUNDLE && n->vname != NULL)
#define uarray_node(n)  (n->t == ARRAY && n->vname != NULL)
#define ulist_node(n)   (n->t == LIST && n->vname != NULL)
#define ustring_node(n) (n->t == STR && n->vname != NULL)
#define useries_node(n) (n->t == SERIES && n->vnum >= 0)
#define uvar_node(n)    (n->vname != NULL)

#define scalar_matrix_node(n) (n->t == MAT && gretl_matrix_is_scalar(n->v.m))
#define scalar_node(n) (n->t == NUM || scalar_matrix_node(n))

#define stringvec_node(n) (n->flags & SVL_NODE)

#define empty_or_num(n) (n == NULL || n->t == EMPTY || n->t == NUM)
#define empty_or_string(n) (n == NULL || n->t == EMPTY || n->t == STR)
#define null_or_empty(n) (n == NULL || n->t == EMPTY)

#define ok_bundled_type(t) (t == NUM || t == STR || t == MAT || t == LIST || \
			    t == SERIES || t == BUNDLE || t == ARRAY)

#define compiled(p) (p->flags & P_EXEC)

static void parser_init (parser *p, const char *str, DATASET *dset,
			 PRN *prn, int flags, int targtype);
static void parser_reinit (parser *p, DATASET *dset, PRN *prn);
static NODE *eval (NODE *t, parser *p);
static void node_type_error (int ntype, int argnum, int goodt,
			     NODE *bad, parser *p);
static int node_is_true (NODE *n, parser *p);
static gretl_matrix *series_to_matrix (const double *x, parser *p,
				       int *prechecked);
static NODE *object_var_node (NODE *t, parser *p);
static void printnode (NODE *t, parser *p, int value);
static inline int attach_aux_node (NODE *t, NODE *ret, parser *p);
static char *get_opstr (int op);

static const char *typestr (int t)
{
    switch (t) {
    case NUM:
	return "scalar";
    case SERIES:
	return "series";
    case MAT:
	return "matrix";
    case STR:
	return "string";
    case U_ADDR:
	return "address";
    case LIST:
	return "list";
    case BUNDLE:
	return "bundle";
    case DBUNDLE:
	return "$-bundle";
    case ARRAY:
	return "array";
    case USERIES:
	return "named series";
    case EMPTY:
	return "empty";
    default:
	return "?";
    }
}

static void free_mspec (matrix_subspec *spec, parser *p)
{
    if (spec != NULL) {
	free(spec->rslice);
	free(spec->cslice);
	free(spec);
    }
}

static void clear_mspec (matrix_subspec *spec, parser *p)
{
    if (spec != NULL) {
	/* the slice elements may not be reusable as is */
	if (spec->rslice != NULL) {
	    free(spec->rslice);
	    spec->rslice = NULL;
	}
	if (spec->cslice != NULL) {
	    free(spec->cslice);
	    spec->cslice = NULL;
	}
    }
}

#if EDEBUG || LHDEBUG

static void print_tree (NODE *t, parser *p, int level)
{
    if (t == NULL) {
	fputs("null\n", stderr);
	return;
    }

    if (bnsym(t->t)) {
	int i;

	for (i=0; i<t->v.bn.n_nodes; i++) {
	    print_tree(t->v.bn.n[i], p, level+1);
	}
    } else if (b3sym(t->t)) {
	print_tree(t->v.b3.l, p, level+1);
	print_tree(t->v.b3.m, p, level+1);
	print_tree(t->v.b3.r, p, level+1);
    } else if (b2sym(t->t)) {
	print_tree(t->v.b2.l, p, level+1);
	print_tree(t->v.b2.r, p, level+1);
    } else if (b1sym(t->t)) {
	print_tree(t->v.b1.b, p, level+1);
    }

    if (t->vname != NULL) {
	fprintf(stderr, " %d: node at %p (type %03d, %s, flags %d), vname='%s'",
		level, (void *) t, t->t, getsymb(t->t), t->flags, t->vname);
	if (t->t == NUM) {
	    fprintf(stderr, ", val %g\n", t->v.xval);
	} else {
	    fputs("\n", stderr);
	}
    } else if (t->t == STR) {
	fprintf(stderr, " %d: node at %p (type %03d, %s, flags %d, val '%s')\n",
		level, (void *) t, t->t, getsymb(t->t), t->flags, t->v.str);
    } else if (t->t == NUM) {
	fprintf(stderr, " %d: node at %p (type %03d, %s, flags %d, val %g)\n",
		level, (void *) t, t->t, getsymb(t->t), t->flags, t->v.xval);
    } else {
	fprintf(stderr, " %d: node at %p (type %03d, %s, flags %d)\n",
		level, (void *) t, t->t, getsymb(t->t), t->flags);
    }

    if (t->aux != NULL) {
	fprintf(stderr, "  aux node at %p (type %03d, %s, flags %d)\n",
		(void *) t->aux, t->aux->t, getsymb(t->aux->t), t->aux->flags);
    }	
}

#endif /* EDEBUG */

#if EDEBUG

static const char *free_tree_tag (int t)
{
    if (t == FR_TREE) {
	return "free tree";
    } else if (t == FR_RET) {
	return "free ret";
    } else if (t == FR_LHTREE) {
	return "free lhtree";
    } else if (t == FR_LHRES) {
	return "free lhres";
    } else {
	return "free other";
    }
}

#endif /* EDEBUG */

/* used when we know that @t is a terminal node: skip
   the tests for attached tree */

static void free_node (NODE *t, parser *p)
{
    if (t->refcount > 1) {
	rndebug(("free node %p (%s): decrement refcount to %d\n",
		 (void *) t, getsymb(t->t), t->refcount - 1));
	t->refcount -= 1;
	return;
    }

    if (is_tmp_node(t)) {
#if EDEBUG
	fprintf(stderr, " tmp node: freeing attached data\n");
#endif
	if (t->t == SERIES) {
	    free(t->v.xvec);
	} else if (t->t == LIST || t->t == IVEC) {
	    free(t->v.ivec);
	} else if (t->t == MAT) {
	    gretl_matrix_free(t->v.m);
	} else if (t->t == MSPEC) {
	    free_mspec(t->v.mspec, p);
	} else if (t->t == BUNDLE) {
	    gretl_bundle_destroy(t->v.b);
	} else if (t->t == ARRAY) {
	    gretl_array_destroy(t->v.a);
	} else if (t->t == STR) {
	    free(t->v.str);
	}
    }

    if (t->t == UOBJ || t->t == WLIST) {
	free(t->v.str);
    }

    if (t->vname != NULL) {
	free(t->vname);
    }

    if (p != NULL && t == p->ret) {
	p->ret = NULL;
    }

    free(t);
}

/* A word on "aux" nodes. These come in two sorts, which
   might be described as "robust" and "fragile" respectively.

   A robust node (identified by the TMP_NODE flag) is one
   whose data pointer is independently allocated. With such
   a node it's OK simply to "pass on" the pointer in
   assignment, and if it's not passed on it should be freed
   on completion of "genr". (So nota bene: if it's assigned
   elsewhere, the pointer on the aux node itself must then
   be set to NULL to avoid double-freeing.)

   A fragile node is one whose data pointer is not
   independently allocated; it actually "belongs to someone
   else". In assignment, then, it must be deeply copied,
   and it must _not_ be freed on completion of genr.

   Obviously, it's necessary to be careful in handling
   fragile nodes, but the advantage of allowing them
   is that they cut down on wasteful deep-copying of objects
   that may be used in calculation, without being modified,
   on the fly.
*/

static void free_tree (NODE *t, parser *p, int code)
{
    if (t == NULL) {
	return;
    }

#if EDEBUG
    fprintf(stderr, "%-11s: starting at %p (type %03d, %s)\n",
	    free_tree_tag(code), (void *) t, t->t,
	    getsymb(t->t));
#endif

    /* free recursively */
    if (bnsym(t->t)) {
	int i;

	for (i=0; i<t->v.bn.n_nodes; i++) {
	    free_tree(t->v.bn.n[i], p, code);
	}
	free(t->v.bn.n);
    } else if (b3sym(t->t)) {
	free_tree(t->v.b3.l, p, code);
	free_tree(t->v.b3.m, p, code);
	free_tree(t->v.b3.r, p, code);
    } else if (b2sym(t->t)) {
	if (!(t->flags & LHT_NODE)) {
	    free_tree(t->v.b2.l, p, code);
	    free_tree(t->v.b2.r, p, code);
	}
    } else if (b1sym(t->t)) {
	free_tree(t->v.b1.b, p, code);
    }

    if (t->aux != NULL && t->aux != p->ret && t->aux != p->lhres) {
	rndebug(("freeing aux node at %p\n", (void *) t->aux));
	free_node(t->aux, p);
    } else if (t->aux != NULL) {
	rndebug(("NOT freeing aux at %p (= p->ret)\n", (void *) t->aux));
	t->aux->refcount -= 1;
    }

#if EDEBUG
    fprintf(stderr, "%-11s: freeing node at %p (type %03d, %s, flags = %d)\n",
	    free_tree_tag(code), (void *) t, t->t, getsymb(t->t),
	    t->flags);
#endif

    free_node(t, p);
}

static void clear_uvnodes (NODE *t)
{
    if (t == NULL) {
	return;
    }

    if (bnsym(t->t)) {
	int i;

	for (i=0; i<t->v.bn.n_nodes; i++) {
	    clear_uvnodes(t->v.bn.n[i]);
	}
    } else if (b3sym(t->t)) {
	clear_uvnodes(t->v.b3.l);
	clear_uvnodes(t->v.b3.m);
	clear_uvnodes(t->v.b3.r);
    } else if (b2sym(t->t)) {
	clear_uvnodes(t->v.b2.l);
	clear_uvnodes(t->v.b2.r);
    } else if (b1sym(t->t)) {
	clear_uvnodes(t->v.b1.b);
    }

    if (t->t == SERIES) {
	if (t->vnum >= 0 || t->vname != NULL) {
#if EDEBUG
	    fprintf(stderr, " clear_uvnode: series at %p\n", (void *) t);
#endif
	    t->v.xvec = NULL;
	}
    } else if (t->uv != NULL) {
#if EDEBUG
	fprintf(stderr, " clear_uvnode: uvar '%s' at %p, uv %p\n", t->vname,
		(void *) t, (void *) t->uv);
#endif
	t->uv = NULL;
    }
}

#if 0 /* may still be useful for debugging? */

static int in_tree (NODE *t, NODE *n)
{
    if (t == NULL || n == NULL) {
	return 0;
    }

    if (n == t || n == t->aux) {
	return 1;
    }

    if (bnsym(t->t)) {
	int i;

	for (i=0; i<t->v.bn.n_nodes; i++) {
	    if (in_tree(t->v.bn.n[i], n)) {
		return 1;
	    }
	}
    } else if (b3sym(t->t)) {
	if (in_tree(t->v.b3.l, n)) return 1;
	if (in_tree(t->v.b3.m, n)) return 1;
	if (in_tree(t->v.b3.r, n)) return 1;
    } else if (b2sym(t->t)) {
	if (in_tree(t->v.b2.l, n)) {
	    return 1;
	}
	if (in_tree(t->v.b2.r, n)) {
	}
    } else if (b1sym(t->t)) {
	if (in_tree(t->v.b1.b, n)) return 1;
    }

    return 0;
}

#endif

#if AUX_NODES_DEBUG
static void reset_p_aux (parser *p, NODE *n)
{
    fprintf(stderr, "resetting p->aux = %p\n", (void *) n);
    p->aux = n;
}
#else
# define reset_p_aux(p, n) (p->aux = n)
#endif

static NODE *newmdef (int k)
{
    NODE *n = new_node(MDEF);

    if (n != NULL) {
	int i;

	if (k > 0) {
	    n->v.bn.n = malloc(k * sizeof n);
	    if (n->v.bn.n != NULL) {
		for (i=0; i<k; i++) {
		    n->v.bn.n[i] = NULL;
		}
	    } else {
		free(n);
		n = NULL;
	    }
	} else {
	    n->v.bn.n = NULL;
	}
	if (n != NULL) {
	    n->v.bn.n_nodes = k;
	}
    }

    return n;
}

static double *na_array (int n)
{
    double *x = malloc(n * sizeof *x);
    int i;

    if (x != NULL) {
	for (i=0; i<n; i++) {
	    x[i] = NADBL;
	}
    }

    return x;
}

/* new node to hold array of doubles */

static NODE *newseries (int n, int flags)
{
    NODE *b = new_node(SERIES);

    if (b != NULL) {
	b->flags = flags;
	if (n > 0) {
	    b->v.xvec = na_array(n);
	    if (b->v.xvec == NULL) {
		free(b);
		b = NULL;
	    }
	} else {
	    b->v.xvec = NULL;
	}
    }

    return b;
}

/* new node to hold array of @n ints */

static NODE *newivec (int n)
{
    NODE *b = new_node(IVEC);

    if (b != NULL) {
	b->flags = TMP_NODE;
	if (n > 0) {
	    b->v.ivec = malloc(n * sizeof(int));
	    if (b->v.ivec == NULL) {
		free(b);
		b = NULL;
	    }
	} else {
	    b->v.ivec = NULL;
	}
    }

    return b;
}

/* new node to hold a gretl_matrix */

static NODE *newmat (int flags)
{
    NODE *n = new_node(MAT);

    if (n != NULL) {
	n->flags = flags;
	n->v.m = NULL;
    }

    return n;
}

/* new node to hold a matrix specification */

static NODE *newmspec (void)
{
    NODE *n = new_node(MSPEC);

    if (n != NULL) {
	n->flags = TMP_NODE;
	n->v.mspec = NULL;
    }

    return n;
}

/* new node to hold a list */

static NODE *newlist (int flags)
{
    NODE *n = new_node(LIST);

    if (n != NULL) {
	n->flags = flags;
	n->v.ivec = NULL;
    }

    return n;
}

static int node_allocate_matrix (NODE *t, int m, int n, parser *p)
{
#if AUX_NODES_DEBUG
    void *in = t->v.m;
#endif

    if (m == 0 || n == 0) {
	t->v.m = gretl_null_matrix_new();
    } else {
	t->v.m = gretl_matrix_alloc(m, n);
    }

#if AUX_NODES_DEBUG
    fprintf(stderr, "%p: node_allocate_matrix: in %p, out %p\n",
	    (void *) t, in, (void *) t->v.m);
#endif

    if (t->v.m == NULL) {
	p->err = E_ALLOC;
    }

    return p->err;
}

static NODE *newstring (int flags)
{
    NODE *n = new_node(STR);

    if (n != NULL) {
	n->flags = flags;
	n->v.str = NULL;
    }

    return n;
}

static NODE *newbundle (int flags)
{
    NODE *n = new_node(BUNDLE);

    if (n != NULL) {
	n->flags = flags;
	n->v.b = NULL;
    }

    return n;
}

static NODE *newarray (int flags)
{
    NODE *n = new_node(ARRAY);

    if (n != NULL) {
	n->flags = flags;
	n->v.a = NULL;
    }

    return n;
}

static void clear_tmp_node_data (NODE *n, parser *p)
{
    int nullify = 1;

    if (n->t == LIST) {
	free(n->v.ivec);
    } else if (n->t == MAT) {
	gretl_matrix_free(n->v.m);
    } else if (n->t == MSPEC) {
	clear_mspec(n->v.mspec, p);
	nullify = 0;
    } else if (n->t == BUNDLE) {
	gretl_bundle_destroy(n->v.b);
    } else if (n->t == ARRAY) {
	gretl_array_destroy(n->v.a);
    } else if (n->t == STR) {
	free(n->v.str);
    } else if (n->t == SERIES) {
	/* preserve any existing tmp series, unless the
	   dataset series length has changed
	*/
	if (p->flags & P_DELTAN) {
	    free(n->v.xvec);
	    n->v.xvec = NULL;
	    if (p->dset_n > 0) {
		n->v.xvec = na_array(p->dset_n);
		if (n->v.xvec == NULL) {
		    p->err = E_ALLOC;
		}
	    }
	} else {
	    /* scrub any pre-existing values in the current
	       sample range */
	    int t;

	    for (t=p->dset->t1; t<=p->dset->t2; t++) {
		n->v.xvec[t] = NADBL;
	    }
	}
	nullify = 0;
    } else {
	nullify = 0;
    }

    if (nullify) {
	n->v.ptr = NULL;
    }
}

/* We allow here for some equivocation in type between
   1 x 1 matrices and scalars in the course of executing
   a compiled parser with saved aux nodes.
*/

static void maybe_switch_node_type (NODE *n, int type,
				    int flags, parser *p)
{
    if (n->t == MAT && type == NUM) {
	/* switch aux node @n from matrix to scalar */
	if (is_tmp_node(n)) {
	    gretl_matrix_free(n->v.m);
	}
	n->t = NUM;
	n->v.xval = NADBL;
	n->flags = 0;
	n->vnum = NO_VNUM;
	n->vname = NULL;
    } else if (n->t == NUM && type == MAT) {
	/* switch @n from scalar to matrix */
	n->t = MAT;
	n->v.m = NULL;
	n->flags = flags;
    } else if (type == EMPTY) {
	; /* LHS mechanism: OK */
    } else {
	/* any other discrepancy presumably means that
	   things have gone badly wrong
	*/
	fprintf(stderr, "aux node mismatch: n->t = %d (%s), type = %d (%s), tmp = %d\n",
		n->t, getsymb(n->t), type, getsymb(type), (flags == TMP_NODE));
	gretl_errmsg_set("internal genr error: aux node mismatch");
	p->err = E_DATA;
    }
}

/* get an auxiliary node: if starting from scratch we allocate
   a new node, otherwise we look up an existing one */

static NODE *get_aux_node (parser *p, int t, int n, int flags)
{
    NODE *ret = p->aux;

#if EDEBUG
    fprintf(stderr, "get_aux_node: p=%p, t=%s, tmp=%d, starting=%d, "
	    "p->aux=%p\n", (void *) p, getsymb(t),
	    (flags & TMP_NODE)? 1 : 0, starting(p) ? 1 : 0,
	    (void *) p->aux);
#endif

    if (is_proxy_node(ret)) {
	/* this node will get freed later */
	ret = NULL;
    }

    if (ret != NULL) {
	/* got a pre-existing aux node */
	if (starting(p)) {
	    if (ret->t != t) {
		maybe_switch_node_type(ret, t, flags, p);
	    } else if (is_tmp_node(ret) && !(p->flags & P_MSAVE)) {
		clear_tmp_node_data(ret, p);
	    }
	}
    } else {
	/* we need to create a new aux node */
	if (t == NUM) {
	    ret = newdbl(NADBL);
	} else if (t == SERIES) {
	    ret = newseries(n, flags);
	} else if (t == IVEC) {
	    ret = newivec(n);
	} else if (t == LIST) {
	    ret = newlist(flags);
	} else if (t == MAT) {
	    ret = newmat(flags);
	} else if (t == MSPEC) {
	    ret = newmspec();
	} else if (t == MDEF) {
	    ret = newmdef(n);
	} else if (t == STR) {
	    ret = newstring(flags);
	} else if (t == BUNDLE) {
	    ret = newbundle(flags);
	} else if (t == ARRAY) {
	    ret = newarray(flags);
	} else if (t == EMPTY) {
	    ret = newempty();
	} else {
	    /* invalid aux node spec */
	    p->err = E_DATA;
	}

	if (!p->err && ret == NULL) {
	    p->err = E_ALLOC;
	}

	if (!p->err) {
	    ret->flags |= AUX_NODE;
	}
    }

    return ret;
}

static NODE *aux_scalar_node (parser *p)
{
    return get_aux_node(p, NUM, 0, 0);
}

static int no_data_error (parser *p)
{
    p->err = E_NODATA;
    return E_NODATA;
}

static NODE *aux_series_node (parser *p)
{
    if (p->dset == NULL || p->dset->n == 0) {
	no_data_error(p);
	return NULL;
    } else {
	return get_aux_node(p, SERIES, p->dset->n, TMP_NODE);
    }
}

static NODE *aux_empty_series_node (parser *p)
{
    if (p->dset == NULL || p->dset->n == 0) {
	no_data_error(p);
	return NULL;
    } else {
	return get_aux_node(p, SERIES, 0, TMP_NODE);
    }
}

static NODE *aux_ivec_node (parser *p, int n)
{
    return get_aux_node(p, IVEC, n, TMP_NODE);
}

static NODE *aux_list_node (parser *p)
{
    if (p->dset == NULL || p->dset->n == 0) {
	no_data_error(p);
	return NULL;
    } else {
	return get_aux_node(p, LIST, 0, TMP_NODE);
    }
}

static NODE *list_pointer_node (parser *p)
{
    if (p->dset == NULL || p->dset->n == 0) {
	no_data_error(p);
	return NULL;
    } else {
	return get_aux_node(p, LIST, 0, 0);
    }
}

static NODE *aux_matrix_node (parser *p)
{
    return get_aux_node(p, MAT, 0, TMP_NODE);
}

static NODE *matrix_pointer_node (parser *p)
{
    return get_aux_node(p, MAT, 0, 0);
}

static NODE *aux_mspec_node (parser *p)
{
    return get_aux_node(p, MSPEC, 0, 0);
}

static NODE *aux_string_node (parser *p)
{
    return get_aux_node(p, STR, 0, TMP_NODE);
}

static NODE *string_pointer_node (parser *p)
{
    return get_aux_node(p, STR, 0, 0);
}

static NODE *aux_bundle_node (parser *p)
{
    return get_aux_node(p, BUNDLE, 0, TMP_NODE);
}

static NODE *bundle_pointer_node (parser *p)
{
    return get_aux_node(p, BUNDLE, 0, 0);
}

static NODE *aux_array_node (parser *p)
{
    return get_aux_node(p, ARRAY, 0, TMP_NODE);
}

static NODE *array_pointer_node (parser *p)
{
    return get_aux_node(p, ARRAY, 0, 0);
}

static NODE *aux_b2_node (parser *p)
{
    return get_aux_node(p, EMPTY, 0, 0);
}

static NODE *aux_any_node (parser *p)
{
    return get_aux_node(p, 0, 0, 0);
}

/* Start of functions that probably should not be needed in
   their present full form, but testing is required before
   they're slimmed down. The general idea is that we should
   already have the user_var pointer we're in need of without
   having to look up it by name, again. However, we'll fall
   back to name look-up (and squawk about it on stderr) if
   need be.
*/

static void *gen_get_lhs_var (parser *p, GretlType type)
{
    void *data = NULL;

    if (p->lh.uv != NULL && p->lh.uv->type == type) {
	data = p->lh.uv->ptr;
    } else {
	if (p->lh.uv == NULL) {
	    fprintf(stderr, "*** get: LHS %s '%s' is NULL!\n",
		    gretl_type_get_name(type), p->lh.name);
	} else {
	    fprintf(stderr, "*** get: LHS uv '%s' of wrong type!\n",
		    p->lh.name);
	}
	if (type == GRETL_TYPE_BUNDLE) {
	    data = get_bundle_by_name(p->lh.name);
	} else if (type == GRETL_TYPE_ARRAY) {
	    data = get_array_by_name(p->lh.name);
	} else if (type == GRETL_TYPE_MATRIX) {
	    data = get_matrix_by_name(p->lh.name);
	} else if (type == GRETL_TYPE_STRING) {
	    data = get_string_by_name(p->lh.name);
	} else if (type == GRETL_TYPE_LIST) {
	    data = get_list_by_name(p->lh.name);
	}
    }

    return data;
}

static int gen_type_from_gretl_type (GretlType t)
{
    switch (t) {
    case GRETL_TYPE_DOUBLE:
	return NUM;
	break;
    case GRETL_TYPE_MATRIX:
	return MAT;
	break;
    case GRETL_TYPE_LIST:
	return LIST;
	break;
    case GRETL_TYPE_STRING:
	return STR;
	break;
    case GRETL_TYPE_BUNDLE:
	return BUNDLE;
	break;
    case GRETL_TYPE_ARRAY:
	return ARRAY;
	break;
    default:
	return UNDEF;
	break;
    }
}

static GretlType gretl_type_from_gen_type (int t)
{
    switch (t) {
    case NUM:
	return GRETL_TYPE_DOUBLE;
	break;
    case MAT:
	return GRETL_TYPE_MATRIX;
	break;
    case LIST:
	return GRETL_TYPE_LIST;
	break;
    case STR:
	return GRETL_TYPE_STRING;
	break;
    case BUNDLE:
	return GRETL_TYPE_BUNDLE;
	break;
    case ARRAY:
	return GRETL_TYPE_ARRAY;
	break;
    default:
	return GRETL_TYPE_NONE;
	break;
    }
}

static NODE *maybe_rescue_undef_node (NODE *n, parser *p)
{
    int v = current_series_index(p->dset, n->vname);
    user_var *uv = NULL;

    if (v >= 0) {
	n->t = SERIES;
	n->vnum = v;
	n->v.xvec = p->dset->Z[v];
	if (is_string_valued(p->dset, n->vnum)) {
	    n->flags |= SVL_NODE;
	}
    } else if ((uv = get_user_var_by_name(n->vname)) != NULL) {
	GretlType type = user_var_get_type(uv);

	n->t = gen_type_from_gretl_type(type);
	n->uv = uv;
	if (type == GRETL_TYPE_DOUBLE) {
	    n->v.xval = *(double *) uv->ptr;
	} else {
	    n->v.ptr = uv->ptr;
	}
    } else {
	undefined_symbol_error(n->vname, p);
    }

    return n;
}

static int gen_add_or_replace (parser *p, GretlType type, void *data)
{
    int err;

    if (p->lh.uv != NULL) {
	err = user_var_replace_value(p->lh.uv, data, type);
    } else {
	err = user_var_add_or_replace(p->lh.name, type, data);
    }

    return err;
}

static int gen_replace_lhs (parser *p, GretlType type, void *data)
{
    if (p->lh.uv == NULL) {
	fputs("*** gen_replace_lhs: lhs user_var is NULL ***\n", stderr);
	fprintf(stderr, " (type is specified as %s)\n",
		gretl_type_get_name(type));
	return E_DATA;
    } else {
	return user_var_replace_value(p->lh.uv, data, type);
    }
}

static int gen_add_uvar (parser *p, GretlType type, void *data)
{
    int err;
    
    err = user_var_add(p->lh.name, type, data);

    /* FIXME attach lh.uv pointer? */
    return err;
}

static int gen_edit_list (parser *p, int *list, int op)
{
    user_var *u;
    int err;

    if (p->lh.uv != NULL && p->lh.uv->type == GRETL_TYPE_LIST) {
	u = p->lh.uv;
    } else {
	if (p->lh.uv == NULL) {
	    fprintf(stderr, "*** replace list: LHS uv is NULL!\n");
	} else {
	    fprintf(stderr, "*** replace list: LHS uv of wrong type!\n");
	}
	u = get_user_var_of_type_by_name(p->lh.name, GRETL_TYPE_LIST);
    }

    if (op == B_ASN) {
	err = replace_list_by_data(u, list);
    } else if (op == B_ADD) {
	err = append_to_list_by_data(u, list);
    } else { /* must be B_SUB */
	err = subtract_from_list_by_data(u, list);
    }

    return err;
}

#ifdef HAVE_MPI

static int node_replace_matrix (NODE *n, gretl_matrix *m)
{
    int err;

    if (n->uv != NULL) {
	err = user_var_replace_value(n->uv, m, GRETL_TYPE_MATRIX);
    } else {
	fprintf(stderr, "*** replace matrix: node uv is NULL!\n");
	err = user_matrix_replace_matrix_by_name(n->vname, m);
    }

    return err;
}

static int node_replace_bundle (NODE *n, gretl_bundle *b)
{
    int err;

    if (n->uv != NULL) {
	err = user_var_replace_value(n->uv, b, GRETL_TYPE_BUNDLE);
    } else {
	fprintf(stderr, "*** replace bundle: node uv is NULL!\n");
	err = E_DATA;
    }

    return err;
}

#endif

static int node_replace_scalar (NODE *n, double x)
{
    int err;

    if (n->uv != NULL && n->uv->type == GRETL_TYPE_DOUBLE) {
	err = user_var_set_scalar_value(n->uv, x);
    } else {
	if (n->uv == NULL) {
	    fprintf(stderr, "*** replace scalar: node uv is NULL!\n");
	} else {
	    fprintf(stderr, "*** replace scalar: node uv of wrong type!\n");
	}
	err = gretl_scalar_set_value(n->vname, x);
    }

    return err;
}

static int gen_replace_scalar (parser *p, double x)
{
    int err;

    if (p->lh.uv != NULL && p->lh.uv->type == GRETL_TYPE_DOUBLE) {
	err = user_var_set_scalar_value(p->lh.uv, x);
    } else {
	if (p->lh.uv == NULL) {
	    fprintf(stderr, "*** replace scalar: LHS uv is NULL!\n");
	} else {
	    fprintf(stderr, "*** replace scalar: LHS uv of wrong type!\n");
	}
	err = gretl_scalar_set_value(p->lh.name, x);
    }

    return err;
}

/* end of functions that can probably be slimmed down */

static void eval_warning (parser *p, int op, int errnum)
{
    if (!check_gretl_warning()) {
	const char *w = (op == B_POW)? "pow" : getsymb(op);
	const char *s = (errnum)? strerror(errnum) : NULL;

	if (s != NULL) {
	    gretl_warnmsg_sprintf("%s: %s", w, s);
	} else {
	    gretl_warnmsg_set(w);
	}
    }
}

/* evaluation of binary operators (yielding x op y) for
   scalar operands (also increment/decrement operators)
*/

static double xy_calc (double x, double y, int op, int targ, parser *p)
{
    double z = NADBL;

#if EDEBUG > 1
    fprintf(stderr, "xy_calc: x = %g, y = %g, op = %d ('%s')\n",
	    x, y, op, getsymb(op));
#endif

#if SERIES_ENSURE_FINITE
    if (targ == SERIES) {
	/* this may be questionable */
	if (isnan(x)) x = NADBL;
	if (isnan(y)) y = NADBL;
    }
#endif

    /* assignment */
    if (op == B_ASN) {
	return y;
    }

    /* testing for presence of NAs? */
    if ((p->flags & P_NATEST) && (na(x) || na(y))) {
	return NADBL;
    }

    /* 0 * NA = NA * 0 = 0 */
    if (op == B_MUL && ((na(x) && y == 0) || (na(y) && x == 0))) {
	return 0;
    }

    /* logical OR: if x or y is valid and non-zero, ignore NA for
       the other term */
    if (op == B_OR && ((!na(x) && x != 0) || (!na(y) && y != 0))) {
	return 1.0;
    }

    /* logical AND: if either x or y is false, the logical product
       should be false, even if the other term is NA */
    if (op == B_AND && (x == 0 || y == 0)) {
	return 0;
    }

    /* otherwise NA propagates to the result */
    if (na(x) || na(y)) {
	return NADBL;
    }

    errno = 0;

    switch (op) {
    case B_ADD:
	return x + y;
    case B_SUB:
	return x - y;
    case B_MUL:
	return x * y;
    case B_DIV:
	return x / y;
    case B_MOD:
	return fmod(x, y);
    case B_AND:
	return x != 0 && y != 0;
    case B_OR:
	return x != 0 || y != 0;
    case B_EQ:
	return x == y;
    case B_NEQ:
	return x != y;
    case B_GT:
	return x > y;
    case B_LT:
	return x < y;
    case B_GTE:
	return x >= y;
    case B_LTE:
	return x <= y;
    case B_POW:
	z = pow(x, y);
	if (errno) {
	    eval_warning(p, op, errno);
	}
	return z;
    default:
	return z;
    }
}

#define randgen(f) (f == F_RANDGEN || f == F_MRANDGEN || f == F_RANDGEN1)

static int check_dist_count (int d, int f, int *np, int *argc)
{
    int err = 0;

    *np = *argc = 0;

    if (d == D_NC_T || d == D_NC_F || d == D_NC_CHISQ) {
	/* non-central t, chisq and F: only CDF, PDF and INVCDF supported */
	if (f == F_PDF || f == F_CDF || f == F_INVCDF) {
	    *np = (d == D_NC_F) ? 3 : 2;
	} else {
	    err = E_INVARG;
	}
    } else if (d == D_UNIFORM || d == D_UDISCRT) {
	/* only RANDGEN is supported */
	if (randgen(f)) {
	    *np = 2; /* min, max */
	} else {
	    err = E_INVARG;
	}
    } else if (d == D_NORMAL) {
	/* all functions supported */
	if (randgen(f)) {
	    *np = 2; /* mu, sigma */
	} else {
	    *np = 0; /* N(0,1) is assumed */
	}
    } else if (d == D_STUDENT) {
	/* Student t: all functions supported */
	*np = 1; /* df */
    } else if (d == D_CHISQ) {
	/* chi-square: all functions supported */
	*np = 1; /* df */
    } else if (d == D_SNEDECOR) {
	/* all functions supported */
	*np = 2; /* dfn, dfd */
    } else if (d == D_GAMMA) {
	/* partial support */
	if (f == F_CRIT) {
	    err = 1;
	} else {
	    *np = 2; /* shape, scale */
	}
    } else if (d == D_BINOMIAL) {
	*np = 2; /* prob, trials */
    } else if (d == D_BINORM) {
	/* bivariate normal: cdf only */
	if (f == F_CDF) {
	    *np = 1; /* rho */
	    *argc = 2; /* note: special */
	} else {
	    err = E_INVARG;
	}
    } else if (d == D_POISSON) {
	*np = 1;
    } else if (d == D_WEIBULL) {
	/* inverse cdf not supported */
	if (f == F_INVCDF) {
	    err = E_INVARG;
	} else {
	    *np = 2; /* shape, scale */
	}
    } else if (d == D_GED) {
	/* GED: critical values not supported */
	if (f == F_CRIT) {
	    err = E_INVARG;
	} else {
	    *np = 1; /* shape */
	}
    } else if (d == D_DW) {
	/* Durbin-Watson: only critical value */
	if (f == F_CRIT) {
	    *np = 2; /* n, k */
	} else {
	    err = E_INVARG;
	}
    } else if (d == D_JOHANSEN) {
	/* Johansen trace test: only p-value */
	if (f == F_PVAL) {
	    *np = 3;
	} else {
	    err = E_INVARG;
	}
    } else if (d == D_BETA) {
	/* randgen only */
	if (randgen(f)) {
	    *np = 2; /* shape1, shape2 */
	} else {
	    err = E_INVARG;
	}
    } else if (d == D_BETABIN) {
	/* randgen only */
	if (randgen(f)) {
	    *np = 3; /* n, shape1, shape2 */
	} else {
	    err = E_INVARG;
	}
    } else {
	err = E_INVARG;
    }

    if (!err && !randgen(f) && *argc == 0) {
	*argc = 1;
    }

    return err;
}

static double scalar_pdist (int t, int d, const double *parm,
			    int np, double arg, parser *p)
{
    double x = NADBL;
    int i;

    for (i=0; i<np; i++) {
	if (na(parm[i])) {
	    return NADBL;
	}
    }

    if (t == F_PVAL) {
	x = gretl_get_pvalue(d, parm, arg);
    } else if (t == F_PDF) {
	x = gretl_get_pdf(d, parm, arg);
    } else if (t == F_CDF) {
	x = gretl_get_cdf(d, parm, arg);
    } else if (t == F_INVCDF) {
	x = gretl_get_cdf_inverse(d, parm, arg);
    } else if (t == F_CRIT) {
	x = gretl_get_critval(d, parm, arg);
    } else {
	p->err = E_PARSE;
    }

    return x;
}

/* @parm contains an array of scalar parameters;
   @argvec contains a series of argument values.
*/

static int series_pdist (double *x, int f, int d,
			 double *parm, int np,
			 const double *argvec,
			 parser *p)
{
    int t;

    if (f == F_PDF) {
	/* fast treatment, for pdf only at this point */
	int n = sample_size(p->dset);

	for (t=p->dset->t1; t<=p->dset->t2; t++) {
	    x[t] = argvec[t];
	}
	gretl_fill_pdf_array(d, parm, x + p->dset->t1, n);
    } else {
	for (t=p->dset->t1; t<=p->dset->t2; t++) {
	    x[t] = scalar_pdist(f, d, parm, np, argvec[t], p);
	}
    }

    return 0;
}

/* @parm contains an array of zero to two scalar parameters;
   @argmat contains an array of argument values.
*/

static gretl_matrix *matrix_pdist (int f, int d,
				   double *parm, int np,
				   gretl_matrix *argmat,
				   parser *p)
{
    gretl_matrix *m;
    double x;
    int i, n;

    if (gretl_is_null_matrix(argmat)) {
	return gretl_null_matrix_new();
    }

    m = gretl_matrix_alloc(argmat->rows, argmat->cols);
    if (m == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

    n = m->rows * m->cols;

    for (i=0; i<n && !p->err; i++) {
	x = scalar_pdist(f, d, parm, np, argmat->val[i], p);
	if (na(x)) {
	    p->err = E_MISSDATA;
	} else {
	    m->val[i] = x;
	}
    }

    if (p->err) {
	gretl_matrix_free(m);
	m = NULL;
    }

    return m;
}

static double node_get_scalar (NODE *n, parser *p)
{
    if (n->t == NUM) {
	return n->v.xval;
    } else if (scalar_matrix_node(n)) {
	return n->v.m->val[0];
    } else {
	p->err = E_INVARG;
	return NADBL;
    }
}

static int node_get_int (NODE *n, parser *p)
{
    double x = node_get_scalar(n, p);

    if (p->err) {
	return -1;
    } else {
	return gretl_int_from_double(x, &p->err);
    }
}

static int node_get_bool (NODE *n, parser *p, int deflt)
{
    int ret = -1;

    if (!null_or_empty(n)) {
	int k = node_get_int(n, p);

	if (!p->err) {
	    ret = (k != 0);
	}
    } else if (deflt == 0 || deflt == 1) {
	ret = deflt;
    } else {
	p->err = E_ARGS;
    }

    return ret;
}

static NODE *DW_node (NODE *r, parser *p)
{
    NODE *s, *e, *ret = NULL;
    NODE *save_aux = p->aux;
    int i, parm[2] = {0};

    for (i=0; i<2 && !p->err; i++) {
	s = r->v.bn.n[i+1];
	if (scalar_node(s)) {
	    parm[i] = node_get_int(s, p);
	} else {
	    e = eval(s, p);
	    if (!p->err) {
		if (scalar_node(e)) {
		    parm[i] = node_get_int(e, p);
		} else {
		    p->err = E_INVARG;
		}
	    }
	}
    }

    if (!p->err && (parm[0] < 6 || parm[1] < 0)) {
	p->err = E_INVARG;
    }

    if (!p->err) {
	reset_p_aux(p, save_aux);
	ret = aux_matrix_node(p);
	if (ret != NULL) {
	    ret->v.m = gretl_get_DW(parm[0], parm[1], &p->err);
	}
    }

    return ret;
}

static NODE *eval_urcpval (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	NODE *save_aux = p->aux;
	NODE *s, *e, *r = n->v.b1.b;
	int i, m = r->v.bn.n_nodes;
	int iargs[3] = {0};
	double tau = NADBL;

	if (m != 4) {
	    p->err = E_INVARG;
	}

	/* need double, int, int, int */

	for (i=0; i<4 && !p->err; i++) {
	    s = r->v.bn.n[i];
	    e = eval(s, p);
	    if (!p->err) {
		if (scalar_node(e)) {
		    if (i == 0) {
			tau = node_get_scalar(e, p);
		    } else {
			iargs[i-1] = node_get_int(e, p);
		    }
		} else {
		    p->err = E_TYPES;
		}
	    }
	}

	if (!p->err) {
	    int nobs = iargs[0];
	    int niv = iargs[1];
	    int itv = iargs[2];

	    reset_p_aux(p, save_aux);
	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		ret->v.xval = get_urc_pvalue(tau, nobs, niv,
					     itv, OPT_NONE);
	    }
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static int get_matrix_size (gretl_matrix *a, gretl_matrix *b,
			    int *r, int *c)
{
    int err = 0;

    /* if both matrices are present, they must be the
       same size */

    if (a != NULL) {
	*r = a->rows;
	*c = b->cols;
	if (b != NULL && (b->rows != *r || b->cols != *c)) {
	    err = E_NONCONF;
	}
    } else if (b != NULL) {
	*r = b->rows;
	*c = b->cols;
    } else {
	*r = *c = 0;
    }

    return err;
}

static NODE *bvnorm_node (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	NODE *save_aux = p->aux;
	double *avec = NULL, *bvec = NULL;
	gretl_matrix *amat = NULL, *bmat = NULL;
	double a, b, args[2];
	double rho = NADBL;
	NODE *e;
	int i, mode = 0;

	for (i=0; i<3 && !p->err; i++) {
	    e = eval(n->v.bn.n[i+1], p);
	    if (p->err) {
		break;
	    }
	    if (scalar_node(e)) {
		if (i == 0) {
		    rho = node_get_scalar(e, p);
		} else {
		    args[i-1] = node_get_scalar(e, p);
		}
	    } else if (i == 1) {
		if (e->t == SERIES) {
		    avec = e->v.xvec;
		} else if (e->t == MAT) {
		    amat = e->v.m;
		}
	    } else if (i == 2) {
		if (e->t == SERIES) {
		    bvec = e->v.xvec;
		} else if (e->t == MAT) {
		    bmat = e->v.m;
		}
	    } else {
		node_type_error(F_CDF, i+1, NUM, e, p);
	    }
	}

	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    if ((avec != NULL && bmat != NULL) ||
		(bvec != NULL && amat != NULL)) {
		p->err = E_INVARG;
	    } else if (avec != NULL || bvec != NULL) {
		mode = 1;
		ret = aux_series_node(p);
	    } else if (amat != NULL || bmat != NULL) {
		mode = 2;
		ret = aux_matrix_node(p);
	    } else {
		mode = 0;
		ret = aux_scalar_node(p);
	    }
	}

	if (p->err) {
	    return ret;
	}

	if (mode == 0) {
	    /* a, b are both scalars */
	    ret->v.xval = bvnorm_cdf(rho, args[0], args[1]);
	} else if (mode == 1) {
	    /* a and/or b are series */
	    int t;

	    for (t=p->dset->t1; t<=p->dset->t2; t++) {
		a = (avec != NULL)? avec[t] : args[0];
		b = (bvec != NULL)? bvec[t] : args[1];
		if (na(a) || na(b)) {
		    ret->v.xvec[t] = NADBL;
		} else {
		    ret->v.xvec[t] = bvnorm_cdf(rho, a, b);
		}
	    }
	} else if (mode == 2) {
	    /* a and/or b are matrices */
	    gretl_matrix *m = NULL;
	    int r, c;

	    p->err = get_matrix_size(amat, bmat, &r, &c);

	    if (!p->err && r > 0 && c > 0) {
		m = gretl_matrix_alloc(r, c);
	    }

	    if (m != NULL) {
		int i, n = r * c;

		for (i=0; i<n && !p->err; i++) {
		    a = (amat != NULL)? amat->val[i] : args[0];
		    b = (bmat != NULL)? bmat->val[i] : args[1];
		    m->val[i] = bvnorm_cdf(rho, a, b);
		    if (na(m->val[i])) {
			/* matrix: change NAs to NaNs */
			m->val[i] = 0.0/0.0;
		    }
		}
	    }

	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }

	    ret->v.m = m;
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

/* return a node containing the evaluated result of a
   probability distribution function */

static NODE *eval_pdist (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	NODE *save_aux = p->aux;
	NODE *e, *s, *r = n->v.b1.b;
	int i, k, m = r->v.bn.n_nodes;
	int rgen = (n->t == F_RANDGEN);
	int mrgen = (n->t == F_MRANDGEN);
	int rgen1 = (n->t == F_RANDGEN1);
	double parm[3] = {0};
	double argval = NADBL;
	double *parmvec[2] = { NULL };
	double *argvec = NULL;
	gretl_matrix *argmat = NULL;
	int rows = 0, cols = 0;
	int d, np, argc, bb;

	if (mrgen) {
	    if (m < 4 || m > 7) {
		p->err = E_INVARG;
		goto disterr;
	    }
	} else if (m < 2 || m > 5) {
	    p->err = E_INVARG;
	    goto disterr;
	}

	s = r->v.bn.n[0];
	if (s->t == STR) {
	    d = dist_code_from_string(s->v.str);
	    if (d == 0) {
		p->err = E_INVARG;
		goto disterr;
	    }
	} else {
	    node_type_error(n->t, 0, STR, s, p);
	    goto disterr;
	}

	p->err = check_dist_count(d, n->t, &np, &argc);
	k = np + argc + 2 * mrgen;
	if (!p->err && k != m - 1) {
	    p->err = E_INVARG;
	}
	if (p->err) {
	    goto disterr;
	}

	bb = (d == D_BETABIN);

	if (d == D_DW) {
	    /* special: Durbin-Watson */
	    return DW_node(r, p);
	} else if (d == D_BINORM) {
	    /* special: bivariate normal */
	    return bvnorm_node(r, p);
	}

	for (i=1; i<=k && !p->err; i++) {
	    s = r->v.bn.n[i];
	    e = eval(s, p);
	    if (p->err) {
		break;
	    }
	    if (scalar_node(e)) {
		/* scalars always acceptable */
		if (mrgen) {
		    if (i == k) {
			cols = node_get_int(e, p);
		    } else if (i == k-1) {
			rows = node_get_int(e, p);
		    } else {
			parm[i-1] = node_get_scalar(e, p);
		    }
		} else if (i == k && argc > 0) {
		    argval = node_get_scalar(e, p);
		} else {
		    parm[i-1] = node_get_scalar(e, p);
		}
	    } else if (i == k && e->t == SERIES) {
		/* a series in the last place? */
		if (bb) {
		    node_type_error(n->t, i, NUM, e, p);
		} else if (rgen) {
		    parmvec[i-1] = e->v.xvec;
		} else if (mrgen) {
		    node_type_error(n->t, i, NUM, e, p);
		} else {
		    argvec = e->v.xvec;
		}
	    } else if (i == k && e->t == MAT) {
		/* a matrix in the last place? */
		if (rgen || mrgen) {
		    node_type_error(n->t, i, NUM, e, p);
		} else {
		    argmat = e->v.m;
		}
	    } else if (e->t == SERIES) {
		/* a series param for randgen? */
		if (rgen && !bb) {
		    parmvec[i-1] = e->v.xvec;
		} else {
		    node_type_error(n->t, i, NUM, e, p);
		}
	    } else {
		p->err = E_INVARG;
		fprintf(stderr, "eval_pdist: arg %d, bad type %d\n", i+1, e->t);
	    }
	}

	if (p->err) {
	    goto disterr;
	}

	reset_p_aux(p, save_aux);

	if (mrgen) {
	    ret = aux_matrix_node(p);
	} else if (rgen || argvec != NULL) {
	    ret = aux_series_node(p);
	} else if (argmat != NULL) {
	    ret = aux_matrix_node(p);
	} else {
	    ret = aux_scalar_node(p);
	}

	if (ret == NULL) {
	    goto disterr;
	}

	if (rgen) {
	    p->err = gretl_fill_random_series(ret->v.xvec, d, parm,
					      parmvec[0], parmvec[1],
					      p->dset);
	} else if (mrgen) {
	    ret->v.m = gretl_get_random_matrix(d, parm, rows, cols,
					       &p->err);
	} else if (rgen1) {
	    ret->v.xval = gretl_get_random_scalar(d, parm, &p->err);
	} else if (argvec != NULL) {
	    p->err = series_pdist(ret->v.xvec, n->t, d, parm, np,
				  argvec, p);
	} else if (argmat != NULL) {
	    ret->v.m = matrix_pdist(n->t, d, parm, np, argmat, p);
	} else {
	    ret->v.xval = scalar_pdist(n->t, d, parm, np, argval, p);
	}
    } else {
	ret = aux_any_node(p);
    }

  disterr:

    return ret;
}

static int mpi_rank = -1;
static int mpi_size = 0;

void set_mpi_rank_and_size (int rank, int size)
{
    mpi_rank = rank;
    mpi_size = size;
}

static double get_const_by_id (int id)
{
    if (id == CONST_PI) {
	return M_PI;
    } else if (id == CONST_EPS) {
	/* see https://en.wikipedia.org/wiki/Machine_epsilon :
	   we now use the (b) definition, as per Matlab, Gauss,
	   R and others
	*/
	return pow(2.0, -52);
    } else if (id == CONST_INF) {
#ifdef INFINITY
	return INFINITY;
#else
	return 1.0/0.0;
#endif
    } else if (id == CONST_NAN) {
#ifdef NAN
	return NAN;
#else
	return 0.0/0.0;
#endif
    } else if (id == CONST_WIN32) {
#ifdef WIN32
	return 1;
#else
	return 0;
#endif
    } else if (id == CONST_HAVE_MPI) {
#ifdef HAVE_MPI
	return check_for_mpiexec();
#else
	return 0;
#endif
    } else if (id == CONST_MPI_RANK) {
	return mpi_rank;
    } else if (id == CONST_MPI_SIZE) {
	return mpi_size;
    } else if (id == CONST_N_PROC) {
	return gretl_n_processors();
    } else {
	return NADBL;
    }
}

/* look up and return numerical values of symbolic constants */

static NODE *retrieve_const (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	ret->v.xval = get_const_by_id(n->v.idnum);
    }

    return ret;
}

double get_const_by_name (const char *name, int *err)
{
    int id = const_lookup(name);

    if (id > 0) {
	return get_const_by_id(id);
    } else {
	if (err != NULL) {
	    *err = E_DATA;
	}
	return NADBL;
    }
}

#ifdef HAVE_MPI

static Gretl_MPI_Op reduce_op_from_string (const char *s)
{
    if (!strcmp(s, "sum")) {
	return GRETL_MPI_SUM;
    } else if (!strcmp(s, "prod")) {
	return GRETL_MPI_PROD;
    } else if (!strcmp(s, "max")) {
	return GRETL_MPI_MAX;
    } else if (!strcmp(s, "min")) {
	return GRETL_MPI_MIN;
    } else if (!strcmp(s, "hcat")) {
	return GRETL_MPI_HCAT;
    } else if (!strcmp(s, "vcat")) {
	return GRETL_MPI_VCAT;
    } else {
	return 0;
    }
}

static Gretl_MPI_Op scatter_op_from_string (const char *s)
{
    if (!strcmp(s, "bycols")) {
	return GRETL_MPI_HSPLIT;
    } else if (!strcmp(s, "byrows")) {
	return GRETL_MPI_VSPLIT;
    } else {
	return 0;
    }
}

#endif

static NODE *mpi_transfer_node (NODE *l, NODE *r, NODE *r2,
				int f, parser *p)
{
#ifndef HAVE_MPI
    gretl_errmsg_set(_("MPI is not supported in this gretl build"));
    p->err = 1;
    return NULL;
#else
    NODE *ret = NULL;
    GretlType type = 0;
    int root = 0;
    int id = 0;

    if (!gretl_mpi_initialized()) {
	gretl_errmsg_set(_("The MPI library is not loaded"));
	p->err = 1;
	return NULL;
    }

    if (f == F_MPI_SEND) {
	/* we support sending a matrix, scalar or bundle; we need
	   the destination id as second argument
	*/
	if (l->t == MAT) {
	    type = GRETL_TYPE_MATRIX;
	} else if (l->t == NUM) {
	    type = GRETL_TYPE_DOUBLE;
	} else if (l->t == BUNDLE) {
	    type = GRETL_TYPE_BUNDLE;
	} else {
	    p->err = E_TYPES;
	}
	if (!p->err) {
	    /* destination id */
	    id = node_get_int(r, p);
	}
    } else if (f == F_MPI_RECV) {
	/* the single argument is the source id */
	id = node_get_int(l, p);
    } else if (f == F_BCAST || f == F_REDUCE ||
	       f == F_ALLREDUCE || f == F_SCATTER) {
	/* we need a variable's address on the left */
	if (l->t != U_ADDR) {
	    p->err = E_TYPES;
	} else {
	    /* switch to 'content' sub-node */
	    l = l->v.b1.b;
	    if (umatrix_node(l)) {
		/* matrix: all operations OK */
		type = GRETL_TYPE_MATRIX;
	    } else if (ubundle_node(l) && f == F_BCAST) {
		/* bundle: only broadcast OK */
		type = GRETL_TYPE_BUNDLE;
	    } else if (uscalar_node(l) && f != F_SCATTER) {
		/* scalar: all ops OK apart from scatter */
		type = GRETL_TYPE_DOUBLE;
	    } else {
		p->err = E_TYPES;
	    }
	}
	if (!p->err && f != F_ALLREDUCE) {
	    /* optional root specification */
	    NODE *rootspec = (f == F_BCAST)? r : r2;

	    if (!null_or_empty(rootspec)) {
		root = node_get_int(rootspec, p);
	    }
	}
	if (!p->err) {
	    /* "self" id */
	    id = gretl_mpi_rank();
	}
    }

    if (p->err) {
	return NULL;
    } else if (f == F_MPI_SEND) {
	void *sendp;

	if (type == GRETL_TYPE_MATRIX) {
	    sendp = l->v.m;
	} else if (type == GRETL_TYPE_BUNDLE) {
	    sendp = l->v.b;
	} else {
	    sendp = &l->v.xval;
	}
	ret = aux_scalar_node(p);
	if (!p->err) {
	    p->err = ret->v.xval = gretl_mpi_send(sendp, type, id);
	}
    } else if (f == F_MPI_RECV) {
	gretl_matrix *m = NULL;
	gretl_bundle *b = NULL;
	double x = NADBL;
	int ival = 0;

	p->err = gretl_mpi_receive(id, &type, &m, &b, &x, &ival);

	if (!p->err) {
	    if (type == GRETL_TYPE_MATRIX) {
		ret = aux_matrix_node(p);
		if (!p->err) {
		    ret->v.m = m;
		}
	    } else if (type == GRETL_TYPE_BUNDLE) {
		ret = aux_bundle_node(p);
		if (!p->err) {
		    ret->v.b = b;
		}
	    } else if (type == GRETL_TYPE_DOUBLE) {
		ret = aux_scalar_node(p);
		if (!p->err) {
		    ret->v.xval = x;
		}
	    } else {
		ret = aux_scalar_node(p);
		if (!p->err) {
		    ret->v.xval = ival;
		}
	    }
	}
    } else if (f == F_BCAST) {
	gretl_matrix *m = NULL;
	gretl_bundle *b = NULL;
	double x = NADBL;
	void *bcastp;

	if (type == GRETL_TYPE_MATRIX) {
	    if (id == root) {
		m = l->v.m;
	    }
	    bcastp = &m;
	} else if (type == GRETL_TYPE_BUNDLE) {
	    if (id == root) {
		b = l->v.b;
	    }
	    bcastp = &b;
	} else {
	    x = l->v.xval;
	    bcastp = &x;
	}

	ret = aux_scalar_node(p);
	if (!p->err) {
	    p->err = gretl_mpi_bcast(bcastp, type, root);
	    if (!p->err && id != root) {
		if (type == GRETL_TYPE_MATRIX) {
		    p->err = node_replace_matrix(l, m);
		} else if (type == GRETL_TYPE_BUNDLE) {
		    p->err = node_replace_bundle(l, b);
		} else {
		    p->err = node_replace_scalar(l, x);
		}
	    }
	    ret->v.xval = p->err;
	}
    } else if (f == F_REDUCE || f == F_ALLREDUCE) {
	ret = aux_scalar_node(p);
	if (!p->err) {
	    Gretl_MPI_Op op = reduce_op_from_string(r->v.str);
	    gretlopt opt = (f == F_REDUCE)? OPT_NONE : OPT_A;
	    gretl_matrix *m = NULL;
	    double x = NADBL;

	    if (type == GRETL_TYPE_MATRIX) {
		p->err = gretl_matrix_mpi_reduce(l->v.m, &m, op, root, opt);
	    } else {
		p->err = gretl_scalar_mpi_reduce(l->v.xval, &x, op, root, opt);
	    }
	    if (!p->err && (id == root || f == F_ALLREDUCE)) {
		if (type == GRETL_TYPE_MATRIX) {
		    p->err = node_replace_matrix(l, m);
		} else {
		    p->err = node_replace_scalar(l, x);
		}
	    }
	    ret->v.xval = p->err;
	}
    } else if (f == F_SCATTER) {
	ret = aux_scalar_node(p);
	if (!p->err) {
	    Gretl_MPI_Op op = scatter_op_from_string(r->v.str);
	    gretl_matrix *m = NULL;

	    p->err = ret->v.xval = gretl_matrix_mpi_scatter(l->v.m, &m,
							    op, root);
	    if (!p->err) {
		p->err = node_replace_matrix(l, m);
	    }
	}
    } else {
	gretl_errmsg_set("MPI function not yet supported");
	p->err = 1;
    }

    return ret;
#endif /* HAVE_MPI */
}

static NODE *scalar_calc (NODE *x, NODE *y, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	ret->v.xval = xy_calc(x->v.xval, y->v.xval, f, NUM, p);
    }

    return ret;
}

static NODE *string_offset (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	int n = strlen(l->v.str);
	int k = r->v.xval;

	if (k < 0) {
	    p->err = E_DATA;
	} else if (k >= n) {
	    ret->v.str = gretl_strdup("");
	} else {
	    ret->v.str = gretl_strdup(l->v.str + k);
	}

	if (!p->err && ret->v.str == NULL) {
	    p->err = E_ALLOC;
	}
    }

    return ret;
}

static NODE *compare_strings (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	int s = strcmp(l->v.str, r->v.str);

	ret->v.xval = (f == B_EQ)? (s == 0) : (s != 0);
    }

    return ret;
}

/*
   We're looking at a comparison, with either a series on the left and
   a string on the right or vice versa.  This can work if the series
   in question is string-valued, as in

     series foo = x=="strval"

   It can also work if the string is an observation marker, as in

     series foo = obs>="CA"
*/

static NODE *series_string_calc (NODE *l, NODE *r, int f, parser *p)
{
    double xt = NADBL, yt = NADBL;
    double *x = NULL, *y = NULL;
    double *alt;
    const char *strval;
    int vnum, t, t1, t2;
    NODE *ret;

    if (r->t == STR) {
	strval = r->v.str;
	vnum = l->vnum;
	x = l->v.xvec;
	alt = &yt;
    } else {
	strval = l->v.str;
	vnum = r->vnum;
	y = r->v.xvec;
	alt = &xt;
    }

    if (vnum > 0 && is_string_valued(p->dset, vnum)) {
	*alt = series_decode_string(p->dset, vnum, strval);
    }

    if (na(*alt)) {
	/* try for an observation string */
	if (annual_data(p->dset)) {
	    *alt = get_date_x(p->dset->pd, strval);
	} else {
	    t = dateton(strval, p->dset);
	    if (t >= 0) {
		*alt = t + 1;
	    }
	}
    }

    if (na(*alt)) {
	gretl_errmsg_sprintf(_("got invalid field '%s'"), strval);
	p->err = E_TYPES;
	return NULL;
    }

    ret = aux_series_node(p);
    if (ret == NULL) {
	return NULL;
    }

    t1 = (autoreg(p))? p->obs : p->dset->t1;
    t2 = (autoreg(p))? p->obs : p->dset->t2;

    for (t=t1; t<=t2; t++) {
	if (x != NULL) {
	    xt = x[t];
	} else if (y != NULL) {
	    yt = y[t];
	}
	ret->v.xvec[t] = xy_calc(xt, yt, f, SERIES, p);
    }

    return ret;
}

static double *list_node_get_series (NODE *n, parser *p)
{
    if (n->v.ivec[0] == 1) {
	int v = n->v.ivec[1];

	if (v >= 0 && v < p->dset->v) {
	    return p->dset->Z[v];
	}
    }

    p->err = E_INVARG;
    return NULL;
}

/* At least one of the nodes is a series; the other may be a
   scalar or 1 x 1 matrix */

static NODE *series_calc (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_series_node(p);
    const double *x = NULL, *y = NULL;
    double xt = 0, yt = 0;

    if (ret == NULL) {
	return NULL;
    }

    if (l->t == SERIES) {
	x = l->v.xvec;
    } else if (l->t == LIST) {
	x = list_node_get_series(l, p);
    } else if (l->t == NUM) {
	xt = l->v.xval;
    } else if (l->t == MAT) {
	xt = l->v.m->val[0];
    }

    if (r->t == SERIES) {
	y = r->v.xvec;
    } else if (r->t == LIST) {
	y = list_node_get_series(r, p);
    } else if (r->t == NUM) {
	yt = r->v.xval;
    } else if (r->t == MAT) {
	yt = r->v.m->val[0];
    }

    if (!p->err) {
	int t1 = (autoreg(p))? p->obs : p->dset->t1;
	int t2 = (autoreg(p))? p->obs : p->dset->t2;
	int t;

	for (t=t1; t<=t2; t++) {
	    if (x != NULL) {
		xt = x[t];
	    }
	    if (y != NULL) {
		yt = y[t];
	    }
	    ret->v.xvec[t] = xy_calc(xt, yt, f, SERIES, p);
	}
    }

    return ret;
}

/* Both nodes are string-valued series */

static NODE *stringvec_calc (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;
    const char *sl, *sr;
    int t, eq;

    if (f != B_EQ && f != B_NEQ) {
	p->err = E_TYPES;
	return NULL;
    }

    ret = aux_series_node(p);
    if (ret == NULL) {
	return NULL;
    }

    for (t=p->dset->t1; t<=p->dset->t2; t++) {
	sl = series_get_string_for_obs(p->dset, l->vnum, t);
	sr = series_get_string_for_obs(p->dset, r->vnum, t);
	if (sl == NULL || sr == NULL) {
	    ret->v.xvec[t] = NADBL;
	} else {
	    eq = strcmp(sl, sr) == 0;
	    ret->v.xvec[t] = (f == B_EQ)? eq : !eq;
	}
    }

    return ret;
}

static int op_symbol (int op)
{
    switch (op) {
    case B_DOTMULT: return '*';
    case B_DOTDIV:  return '/';
    case B_DOTPOW:  return '^';
    case B_DOTADD:  return '+';
    case B_DOTSUB:  return '-';
    case B_DOTEQ:   return '=';
    case B_DOTGT:   return '>';
    case B_DOTLT:   return '<';
    case B_DOTGTE:  return ']';
    case B_DOTLTE:  return '[';
    case B_DOTNEQ:  return '!';
    default: return 0;
    }
}

static gretl_matrix *nullmat_multiply (const gretl_matrix *A,
				       const gretl_matrix *B,
				       int op, int *err)
{
    gretl_matrix *C = NULL;

    if (A->rows == 0 && A->cols == 0 &&
	B->rows == 0 && B->cols == 0) {
	C = gretl_null_matrix_new();
    } else {
	int Lc = op == B_TRMUL ? A->rows : A->cols;
	int Cr = op == B_TRMUL ? A->cols : A->rows;
	int Cc = B->cols;

	if (Lc != B->rows) {
	    *err = E_NONCONF;
	} else {
	    if (Cr > 0 && Cc > 0) {
		C = gretl_zero_matrix_new(Cr, Cc);
	    } else {
		C = gretl_matrix_alloc(Cr, Cc);
	    }
	    if (C == NULL) {
		*err = E_ALLOC;
	    }
	}
    }

    return C;
}

static gretl_matrix *
matrix_add_sub_scalar (const gretl_matrix *A,
		       const gretl_matrix *B,
		       int op)
{
    gretl_matrix *C;
    double xval, *xvec;
    int r, c;

    if (gretl_matrix_is_scalar(A)) {
	r = B->rows;
	c = B->cols;
	xval = A->val[0];
	xvec = B->val;
    } else {
	r = A->rows;
	c = A->cols;
	xval = B->val[0];
	xvec = A->val;
    }

    C = gretl_matrix_alloc(r, c);

    if (C != NULL) {
	int i, n = r * c;

	if (op == B_ADD) {
	    for (i=0; i<n; i++) {
		C->val[i] = xvec[i] + xval;
	    }
	} else {
	    if (xvec == A->val) {
		for (i=0; i<n; i++) {
		    C->val[i] = xvec[i] - xval;
		}
	    } else {
		for (i=0; i<n; i++) {
		    C->val[i] = xval - xvec[i];
		}
	    }
	}
    }

    return C;
}

/* See if we can reuse an existing matrix on an
   auxiliary node. If so, return it; otherwise
   free it and return a newly allocated matrix.
*/

static gretl_matrix *calc_get_matrix (gretl_matrix **pM,
				      int r, int c)
{
    if (*pM == NULL) {
	/* just get on with the allocation */
	return gretl_matrix_alloc(r, c);
    } else if ((*pM)->rows == r && (*pM)->cols == c) {
	/* reusable as-is */
	return *pM;
    } else if ((*pM)->rows == c && (*pM)->cols == r) {
	/* reusable if reoriented */
	(*pM)->rows = r;
	(*pM)->cols = c;
	return *pM;
    } else {
	/* new matrix needed */
	gretl_matrix_free(*pM);
	*pM = NULL;
	return gretl_matrix_alloc(r, c);
    }
}

/* return allocated result of binary operation performed on
   two matrices */

static int real_matrix_calc (const gretl_matrix *A,
			     const gretl_matrix *B,
			     int op, gretl_matrix **pM)
{
    gretl_matrix *C = NULL;
    int ra, ca;
    int rb, cb;
    int r, c;
    int err = 0;

    if (gretl_is_null_matrix(A) ||
	gretl_is_null_matrix(B)) {
	if (op != B_HCAT && op != B_VCAT && op != F_DSUM &&
	    op != B_MUL && op != B_TRMUL) {
	    return E_NONCONF;
	}
	if (op == B_MUL || op == B_TRMUL) {
	    *pM = nullmat_multiply(A, B, op, &err);
	    return err;
	}
    }

    switch (op) {
    case B_ADD:
    case B_SUB:
	if (gretl_matrix_is_scalar(A) ||
	    gretl_matrix_is_scalar(B)) {
	    C = matrix_add_sub_scalar(A, B, op);
	    if (C == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    C = calc_get_matrix(pM, A->rows, A->cols);
	    if (C == NULL) {
		err = E_ALLOC;
	    } else if (op == B_ADD) {
		err = gretl_matrix_add(A, B, C);
	    } else {
		err = gretl_matrix_subtract(A, B, C);
	    }
	}
	break;
    case B_HCAT:
	C = gretl_matrix_col_concat(A, B, &err);
	break;
    case B_VCAT:
	C = gretl_matrix_row_concat(A, B, &err);
	break;
    case F_DSUM:
	C = gretl_matrix_direct_sum(A, B, &err);
	break;
    case B_MUL:
	ra = gretl_matrix_rows(A);
	ca = gretl_matrix_cols(A);
	rb = gretl_matrix_rows(B);
	cb = gretl_matrix_cols(B);

	r = (ra == 1 && ca == 1)? rb : ra;
	c = (rb == 1 && cb == 1)? ca : cb;

	C = calc_get_matrix(pM, r, c);
	if (C == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_matrix_multiply(A, B, C);
	    if (!err) {
		gretl_matrix_transcribe_obs_info(C, A);
	    }
	}
	break;
    case B_TRMUL:
	ra = gretl_matrix_cols(A);
	ca = gretl_matrix_rows(A);
	rb = gretl_matrix_rows(B);
	cb = gretl_matrix_cols(B);

	r = (ra == 1 && ca == 1)? rb : ra;
	c = (rb == 1 && cb == 1)? ca : cb;

	C = calc_get_matrix(pM, r, c);
	if (C == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
					    B, GRETL_MOD_NONE,
					    C, GRETL_MOD_NONE);
	}
	break;
    case F_QFORM:
	/* quadratic form, A * B * A', for symmetric B */
	ra = gretl_matrix_rows(A);
	ca = gretl_matrix_cols(A);
	rb = gretl_matrix_rows(B);
	cb = gretl_matrix_cols(B);

	if (ca != rb || cb != rb) {
	    err = E_NONCONF;
	} else {
	    gretl_matrix_set_equals_tolerance(1.0e-7);
	    if (!gretl_matrix_is_symmetric(B)) {
		gretl_errmsg_set(_("Matrix is not symmetric"));
		err = E_NONCONF;
	    }
	    gretl_matrix_unset_equals_tolerance();
	}
	if (!err) {
	    C = calc_get_matrix(pM, ra, ra);
	    if (C == NULL) {
		err = E_ALLOC;
	    } else {
		err = gretl_matrix_qform(A, GRETL_MOD_NONE, B,
					 C, GRETL_MOD_NONE);
	    }
	}
	break;
    case B_DIV:
    case B_LDIV:
	/* matrix right or left "division" */
	if (op == B_LDIV) {
	    C = gretl_matrix_divide(A, B, GRETL_MOD_NONE, &err);
	} else {
	    /* A/B = (B'\A')' */
	    C = gretl_matrix_divide(A, B, GRETL_MOD_TRANSPOSE, &err);
	}
	break;
    case B_DOTMULT:
    case B_DOTDIV:
    case B_DOTPOW:
    case B_DOTADD:
    case B_DOTSUB:
    case B_DOTEQ:
    case B_DOTGT:
    case B_DOTLT:
    case B_DOTGTE:
    case B_DOTLTE:
    case B_DOTNEQ:
	/* apply operator element-wise */
	C = gretl_matrix_dot_op(A, B, op_symbol(op), &err);
	break;
    case B_KRON:
	/* Kronecker product */
	C = gretl_matrix_kronecker_product_new(A, B, &err);
	break;
    case F_HDPROD:
	C = gretl_matrix_hdproduct_new(A, B, &err);
	break;
    case F_CMULT:
	C = gretl_matrix_complex_multiply(A, B, &err);
	break;
    case F_CDIV:
	C = gretl_matrix_complex_divide(A, B, &err);
	break;
    case F_MRSEL:
	C = gretl_matrix_bool_sel(A, B, 1, &err);
	break;
    case F_MCSEL:
	C = gretl_matrix_bool_sel(A, B, 0, &err);
	break;
    default:
	err = E_TYPES;
	break;
    }

    if (err) {
	if (C != NULL) {
	    if (pM != NULL && *pM == C) {
		*pM = NULL;
	    }
	    gretl_matrix_free(C);
	    C = NULL;
	}
    } else {
	/* preserve data-row info? */
	int At1 = gretl_matrix_get_t1(A);
	int At2 = gretl_matrix_get_t2(A);
	int Bt1 = gretl_matrix_get_t1(B);
	int Bt2 = gretl_matrix_get_t2(B);

	if (C->rows == A->rows && At1 >= 0 && At2 > At1) {
	    gretl_matrix_set_t1(C, At1);
	    gretl_matrix_set_t2(C, At2);
	} else if (C->rows == B->rows && Bt1 >= 0 && Bt2 > Bt1) {
	    gretl_matrix_set_t1(C, Bt1);
	    gretl_matrix_set_t2(C, Bt2);
	}
    }

    if (*pM != NULL && *pM != C) {
	/* we neither freed nor reused *pM */
	gretl_matrix_free(*pM);
    }

    *pM = C;

    return err;
}

static gretl_matrix *tmp_matrix_from_series (NODE *n, parser *p)
{
    int t, T = sample_size(p->dset);
    const double *x = n->v.xvec;
    gretl_matrix *m = NULL;

    m = gretl_column_vector_alloc(T);

    if (m == NULL) {
	p->err = E_ALLOC;
    } else {
	int i = 0;

	for (t=p->dset->t1; t<=p->dset->t2; t++) {
	    if (na(x[t])) {
		m->val[i++] = M_NA;
	    } else {
		m->val[i++] = x[t];
	    }
	}
    }

    return m;
}

const double *get_colvec_as_series (NODE *n, int f, parser *p)
{
    if (n->t != MAT) {
	node_type_error(f, 1, SERIES, n, p);
	return NULL;
    } else {
	const gretl_matrix *m = n->v.m;

	if (m->rows == p->dset->n && m->cols == 1) {
	    return m->val;
	} else {
	    node_type_error(f, 1, SERIES, n, p);
	    return NULL;
	}
    }
}

/* One of the operands is a matrix, the other a series: we
   try "casting" the series to a matrix.
*/

static NODE *matrix_series_calc (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *a, *b, *tmp;

	if (l->t == SERIES) {
	    tmp = a = tmp_matrix_from_series(l, p);
	    b = r->v.m;
	} else {
	    a = l->v.m;
	    tmp = b = tmp_matrix_from_series(r, p);
	}

	if (!p->err) {
	    p->err = real_matrix_calc(a, b, op, &ret->v.m);
	}

	gretl_matrix_free(tmp);
    }

    return ret;
}

/* Here we know have a scalar and a 1 x 1 matrix to work with,
   in either order */

static NODE *matrix_scalar_calc2 (NODE *l, NODE *r, int op,
				  parser *p)
{
    NODE *ret;

    if (l->t == NUM && (op == B_MOD || op == B_POW)) {
	/* the matrix on the right is functioning as
	   a scalar argument, so produce a scalar
	*/
	ret = aux_scalar_node(p);
    } else {
	/* one of the operands is a matrix, albeit 1 x 1,
	   so it's safer to produce a matrix result
	*/
	ret = aux_matrix_node(p);
    }

    if (!p->err) {
	double x, y;

	if (l->t == NUM) {
	    x = l->v.xval;
	    y = r->v.m->val[0];
	} else {
	    x = l->v.m->val[0];
	    y = r->v.xval;
	}

	if (ret->t == NUM) {
	    ret->v.xval = xy_calc(x, y, op, NUM, p);
	} else {
	    node_allocate_matrix(ret, 1, 1, p);
	    if (!p->err) {
		ret->v.m->val[0] = xy_calc(x, y, op, MAT, p);
	    }
	}
    }

    return ret;
}

#define comparison_op(o) (o == B_EQ  || o == B_NEQ || \
			  o == B_LT  || o == B_GT ||  \
			  o == B_LTE || o == B_GTE)

/* Mixed types: one of the operands is a matrix, the other a scalar,
   giving a matrix result unless we're looking at a comparison
   operator.
*/

static NODE *matrix_scalar_calc (NODE *l, NODE *r, int op, parser *p)
{
    const gretl_matrix *m = NULL;
    int comp = comparison_op(op);
    double x;
    NODE *ret = NULL;

    if (op != B_TRMUL) {
	/* Check for the simple case of scalar and
	   1 x 1 matrix, either way round
	*/
	if ((l->t == NUM && scalar_node(r)) ||
	    (r->t == NUM && scalar_node(l))) {
	    return matrix_scalar_calc2(l, r, op, p);
	}
    }

    /* get a scalar @x and matrix @m */
    x = (l->t == NUM)? l->v.xval : r->v.xval;
    m = (l->t == MAT)? l->v.m : r->v.m;

    if (gretl_is_null_matrix(m)) {
	p->err = E_DATA;
	return NULL;
    }

    /* mod, pow: the right-hand term must be scalar */
    if ((op == B_MOD || op == B_POW) && !scalar_node(r)) {
	p->err = E_TYPES;
	return NULL;
    }

    if (comp) {
	ret = aux_scalar_node(p);
    } else {
	ret = aux_matrix_node(p);
    }

    if (ret == NULL) {
	return NULL;
    }

    if (op == B_POW) {
	/* note: the (scalar, 1x1) and (1x1, scalar) cases are
	   handled above
	*/
	int s = node_get_int(r, p);

	if (!p->err) {
	    ret->v.m = gretl_matrix_pow(m, s, &p->err);
	}
	return ret;
    } else if (op == B_TRMUL) {
	gretl_matrix *tmp;

	if (l->t == NUM) {
	    tmp = gretl_matrix_copy(m);
	} else {
	    tmp = gretl_matrix_copy_transpose(m);
	}
	if (tmp == NULL) {
	    p->err = E_ALLOC;
	} else {
	    gretl_matrix_multiply_by_scalar(tmp, x);
	    ret->v.m = tmp;
	}
	return ret;
    } else {
	int i, n = m->rows * m->cols;
	double y;

	if (comp) {
	    ret->v.xval = 1;
	    if (l->t == NUM) {
		for (i=0; i<n; i++) {
		    if (xy_calc(x, m->val[i], op, MAT, p) == 0) {
			ret->v.xval = 0;
			break;
		    }
		}
	    } else {
		for (i=0; i<n; i++) {
		    if (xy_calc(m->val[i], x, op, MAT, p) == 0) {
			ret->v.xval = 0;
			break;
		    }
		}
	    }
	} else {
	    if (node_allocate_matrix(ret, m->rows, m->cols, p)) {
		return NULL;
	    }

	    if (gretl_matrix_is_dated(m)) {
		gretl_matrix_set_t1(ret->v.m, gretl_matrix_get_t1(m));
		gretl_matrix_set_t2(ret->v.m, gretl_matrix_get_t2(m));
	    }

	    if (l->t == NUM) {
		for (i=0; i<n; i++) {
		    y = xy_calc(x, m->val[i], op, MAT, p);
		    ret->v.m->val[i] = y;
		}
	    } else {
		for (i=0; i<n; i++) {
		    y = xy_calc(m->val[i], x, op, MAT, p);
		    ret->v.m->val[i] = y;
		}
	    }
	}
    }

    return ret;
}

static NODE *matrix_transpose_node (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	if (is_tmp_node(n)) {
	    /* transpose temp matrix in place */
	    p->err = gretl_matrix_transpose_in_place(n->v.m);
	    ret = n;
	} else {
	    /* create transpose as new matrix */
	    ret = aux_matrix_node(p);
	    if (!p->err) {
		ret->v.m = gretl_matrix_copy_transpose(n->v.m);
		if (ret->v.m == NULL) {
		    p->err = E_ALLOC;
		}
	    }
	}
    } else {
	ret = is_tmp_node(n) ? n : aux_matrix_node(p);
    }

    return ret;
}

/* We're looking at a string argument that is supposed to represent
   a function call: we'll do a rudimentary heuristic check here.
   FIXME this should be more rigorous?
*/

static int is_function_call (const char *s)
{
    if (!strchr(s, '(') || !strchr(s, ')')) {
	return 0;
    } else {
	return 1;
    }
}

static NODE *numeric_jacobian_or_hessian (NODE *l, NODE *m, NODE *r,
					  int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const char *s = m->v.str;
	double eps = 0.0;

	if (!is_function_call(s)) {
	    p->err = E_TYPES;
	    return NULL;
	}

	ret = aux_matrix_node(p);
	if (ret == NULL) {
	    return NULL;
	}

	if (!null_or_empty(r)) {
	    eps = node_get_scalar(r, p);
	}

	if (!p->err) {
	    if (f == F_FDJAC) {
		ret->v.m = user_fdjac(l->v.m, s, eps, p->dset, &p->err);
	    } else {
		ret->v.m = user_numhess(l->v.m, s, eps, p->dset, &p->err);
	    }
	}
    } else {
	ret = aux_matrix_node(p);
    }

    return ret;
}

static NODE *cmatrix_printf (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	const char *s = NULL;

	if (r != NULL && r->t == STR) {
	    s = r->v.str;
	}
	ret->v.xval = complex_matrix_printf(l->v.m, s, p->prn);
    }

    return ret;
}

static gretl_matrix *node_get_matrix (NODE *n, parser *p)
{
    if (n->t == U_ADDR) {
	n = n->v.b1.b;
	if (n->t != MAT) {
	    p->err = E_TYPES;
	    return NULL;
	}
    }

    return n->v.m;
}

static const char *node_get_fncall (NODE *n, parser *p)
{
    const char *ret = NULL;

    if (n->t != STR) {
	p->err = E_TYPES;
    } else {
	ret = n->v.str;
	if (!is_function_call(ret)) {
	    p->err = E_TYPES;
	}
    }

    return ret;
}

static void n_args_error (int k, int n, int f, parser *p)
{
    gretl_errmsg_sprintf( _("Number of arguments (%d) does not "
			    "match the number of\nparameters for "
			    "function %s (%d)"), k, getsymb(f), n);
    p->err = 1;
}

static NODE *BFGS_constrained_max (NODE *t, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *n = t->v.b1.b;
    NODE *ret = NULL;
    NODE *e = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *bounds = NULL;
    const char *sf = NULL;
    const char *sg = NULL;
    int i, k = n->v.bn.n_nodes;

    if (k < 3 || k > 4) {
	n_args_error(k, 3, F_BFGSCMAX, p);
    }

    for (i=0; i<k && !p->err; i++) {
	e = n->v.bn.n[i];
	if (i == 0) {
	    b = node_get_matrix(e, p);
	} else if (i == 1) {
	    e = eval(n->v.bn.n[i], p);
	    if (!p->err) {
		bounds = node_get_matrix(e, p);
	    }
	} else if (i == 2) {
	    sf = node_get_fncall(e, p);
	} else if (i == 3 && !null_or_empty(e)) {
	    sg = node_get_fncall(e, p);
	}
    }

    if (!p->err) {
	reset_p_aux(p, save_aux);
	ret = aux_scalar_node(p);
    }

    if (!p->err) {
	int minimize = (t->flags & ALS_NODE)? 1 : 0;

	ret->v.xval = user_BFGS(b, sf, sg, p->dset, bounds,
				minimize, p->prn, &p->err);
    }

    return ret;
}

static NODE *BFGS_maximize (NODE *l, NODE *m, NODE *r,
			    parser *p, NODE *t)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *b;
	const char *sf = m->v.str;
	const char *sg = NULL;

	b = node_get_matrix(l, p);

	if (!p->err) {
	    if (r->t == STR) {
		sg = r->v.str;
	    } else if (r->t != EMPTY) {
		p->err = E_TYPES;
	    }
	}

	if (!p->err && !is_function_call(sf)) {
	    p->err = E_TYPES;
	}

	if (!p->err && sg != NULL && !is_function_call(sg)) {
	    p->err = E_TYPES;
	}

	if (!p->err && gretl_is_null_matrix(b)) {
	    p->err = E_DATA;
	}

	if (p->err) {
	    return NULL;
	}

	ret = aux_scalar_node(p);

	if (ret != NULL) {
	    int minimize = (t->flags & ALS_NODE)? 1 : 0;

	    ret->v.xval = user_BFGS(b, sf, sg, p->dset, NULL,
				    minimize, p->prn, &p->err);
	}
    } else {
	ret = aux_scalar_node(p);
    }

    return ret;
}

static NODE *deriv_free_node (NODE *l, NODE *m, NODE *r,
			      parser *p, NODE *t)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *b;
	const char *sf = m->v.str;
	double tol = NADBL;
	int maxit = 0;

	b = node_get_matrix(l, p);

	if (!p->err) {
	    if (gretl_is_null_matrix(b)) {
		p->err = E_DATA;
	    } else if (!is_function_call(sf)) {
		p->err = E_TYPES;
	    }
	}

	if (!p->err) {
	    if (scalar_node(r)) {
		if (t->t == F_GSSMAX) {
		    tol = r->v.xval;
		} else {
		    maxit = node_get_int(r, p);
		}
	    } else if (!null_or_empty(r)) {
		p->err = E_TYPES;
	    }
	}

	if (p->err) {
	    return NULL;
	}

	ret = aux_scalar_node(p);

	if (ret != NULL) {
	    int minimize = (t->flags & ALS_NODE)? 1 : 0;
	    MaxMethod method = SIMANN_MAX;

	    if (t->t == F_NMMAX) {
		method = NM_MAX;
	    } else if (t->t == F_GSSMAX) {
		method = GSS_MAX;
	    }

	    ret->v.xval =
		deriv_free_optimize(method, b, sf, maxit, tol,
				    minimize, p->dset, p->prn,
				    &p->err);
	}
    } else {
	ret = aux_scalar_node(p);
    }

    return ret;
}

static void lag_calc (double *y, const double *x,
		      int k, int t1, int t2,
		      int op, double mul,
		      parser *p)
{
    int s, t;

    for (t=t1; t<=t2; t++) {
	s = t - k;
	if (dated_daily_data(p->dset)) {
	    if (s >= 0 && s < p->dset->n) {
		while (s >= 0 && xna(x[s])) {
		    s--;
		}
	    }
	} else if (p->dset->structure == STACKED_TIME_SERIES) {
	    if (s / p->dset->pd != t / p->dset->pd) {
		/* s and t pertain to different units */
		s = -1;
	    }
	}

	if (s >= 0 && s < p->dset->n) {
	    if (op == B_ASN && mul == 1.0) {
		y[t] = x[s];
	    } else if (op == B_ASN) {
		y[t] = mul * x[s];
	    } else if (op == B_ADD) {
		y[t] += mul * x[s];
	    } else {
		p->err = E_DATA;
	    }
	}
    }
}

static NODE *matrix_file_write (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const char *s = m->v.str;
	int export = 0;

	if (!null_or_empty(r)) {
	    export = (r->v.xval != 0);
	}

	ret = aux_scalar_node(p);

	if (ret != NULL) {
	    ret->v.xval = gretl_matrix_write_to_file(l->v.m, s, export);
	}
    } else {
	ret = aux_scalar_node(p);
    }

    return ret;
}

static NODE *bundle_file_write (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const char *s = m->v.str;
	int export = 0;

	if (!null_or_empty(r)) {
	    export = (r->v.xval != 0);
	}

	ret = aux_scalar_node(p);

	if (ret != NULL) {
	    ret->v.xval = gretl_bundle_write_to_file(l->v.b, s, export);
	}
    } else {
	ret = aux_scalar_node(p);
    }

    return ret;
}

/* matrix on left, scalar on right */

static NODE *matrix_scalar_func (NODE *l, NODE *r,
				 int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *m = l->v.m;
	int k = node_get_int(r, p);

	if (gretl_is_null_matrix(m)) {
	    p->err = E_INVARG;
	}

	if (p->err) {
	    return NULL;
	}

	ret = aux_matrix_node(p);
	if (ret == NULL) {
	    return NULL;
	}

	if (f == F_MSORTBY) {
	    ret->v.m = gretl_matrix_sort_by_column(m, k-1, &p->err);
	} else if (f == HF_CXTRACT) {
	    ret->v.m = gretl_cxtract(m, k, &p->err);
	}
    } else {
	ret = aux_matrix_node(p);
    }

    return ret;
}

static gretl_matrix *make_scalar_matrix (double x)
{
    gretl_matrix *m = gretl_matrix_alloc(1, 1);

    if (m != NULL) {
	m->val[0] = x;
    }

    return m;
}

/* both operands are known to be matrices or scalars */

static NODE *matrix_matrix_calc (NODE *l, NODE *r, int op, parser *p)
{
    gretl_matrix *ml = NULL, *mr = NULL;
    NODE *ret;

    if (op == B_DOTPOW || op == B_POW) {
	ret = aux_matrix_node(p);
    } else {
	/* experiment: try reusing aux matrix */
	p->flags |= P_MSAVE;
	ret = get_aux_node(p, MAT, 0, TMP_NODE);
	p->flags ^= P_MSAVE;
    }

#if EDEBUG
    fprintf(stderr, "matrix_matrix_calc: l=%p, r=%p, ret=%p\n",
	    (void *) l, (void *) r, (void *) ret);
#endif

    if (op == B_POW) {
	if (scalar_node(l) && scalar_node(r)) {
	    op = B_DOTPOW;
	} else if (!scalar_node(r)) {
	    p->err = E_TYPES;
	    return NULL;
	}
    }

    if (l->t == NUM) {
	ml = make_scalar_matrix(l->v.xval);
    } else {
	ml = l->v.m;
    }

    if (op != B_POW) {
	if (r->t == NUM) {
	    mr = make_scalar_matrix(r->v.xval);
	} else {
	    mr = r->v.m;
	}
    }

    if (ret != NULL && starting(p)) {
	if (op == B_DOTPOW) {
	    ret->v.m = gretl_matrix_dot_op(ml, mr, '^', &p->err);
	} else if (op == B_POW) {
	    int s = node_get_int(r, p);

	    if (!p->err) {
		ret->v.m = gretl_matrix_pow(ml, s, &p->err);
	    }
	} else {
	    p->err = real_matrix_calc(ml, mr, op, &ret->v.m);
	}
    }

    if (l->t == NUM) gretl_matrix_free(ml);
    if (r->t == NUM) gretl_matrix_free(mr);

    return ret;
}

static NODE *matrix_and_or (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const gretl_matrix *a = l->v.m;
	const gretl_matrix *b = r->v.m;
	int i, n = a->rows * a->cols;

	if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	    p->err = E_NONCONF;
	} else if (a->rows != b->rows || a->cols != b->cols) {
	    p->err = E_NONCONF;
	} else {
	    ret->v.m = gretl_unit_matrix_new(a->rows, a->cols);
	    if (ret->v.m == NULL) {
		p->err = E_ALLOC;
		return NULL;
	    }
	    for (i=0; i<n; i++) {
		if (op == B_AND) {
		    if (a->val[i] == 0.0 || b->val[i] == 0.0) {
			ret->v.m->val[i] = 0.0;
		    }
		} else if (op == B_OR) {
		    if (a->val[i] == 0.0 && b->val[i] == 0.0) {
			ret->v.m->val[i] = 0.0;
		    }
		}
	    }
	}
    }

    return ret;
}

/* both operands are matrices */

static NODE *matrix_bool (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret;

    if (op == B_OR || op == B_AND) {
	return matrix_and_or(l, r, op, p);
    }

    ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	const gretl_matrix *a = l->v.m;
	const gretl_matrix *b = r->v.m;
	int i, n = a->rows * a->cols;

	if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	    ret->v.xval = NADBL;
	} else if (a->rows != b->rows || a->cols != b->cols) {
	    ret->v.xval = NADBL;
	} else {
	    ret->v.xval = op == B_NEQ ? 0 : 1;
	    for (i=0; i<n; i++) {
		if (op == B_EQ && a->val[i] != b->val[i]) {
		    ret->v.xval = 0;
		    break;
		} else if (op == B_LT && a->val[i] >= b->val[i]) {
		    ret->v.xval = 0;
		    break;
		} else if (op == B_GT && a->val[i] <= b->val[i]) {
		    ret->v.xval = 0;
		    break;
		} else if (op == B_LTE && a->val[i] > b->val[i]) {
		    ret->v.xval = 0;
		    break;
		} else if (op == B_GTE && a->val[i] < b->val[i]) {
		    ret->v.xval = 0;
		    break;
		} else if (op == B_NEQ && a->val[i] != b->val[i]) {
		    ret->v.xval = 1;
		    break;
		}
	    }
	}
    }

    return ret;
}

static void matrix_error (parser *p)
{
    if (p->err == 0) {
	p->err = 1;
    }

    if (gretl_errmsg_is_set()) {
	errmsg(p->err, p->prn);
    }
}

/* functions taking a matrix argument and returning a
   scalar result */

static NODE *matrix_to_scalar_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *m;
	int tmpmat = 0;

	if (n->t == NUM) {
	    m = gretl_matrix_from_scalar(node_get_scalar(n, p));
	    tmpmat = 1;
	} else {
	    tmpmat = is_tmp_node(n);
	    m = n->v.m;
	}

	switch (f) {
	case F_ROWS:
	    ret->v.xval = m->rows;
	    break;
	case F_COLS:
	    ret->v.xval = m->cols;
	    break;
	case F_DET:
	case F_LDET:
	    ret->v.xval = user_matrix_get_determinant(m, tmpmat, f, &p->err);
	    break;
	case F_TRACE:
	    ret->v.xval = gretl_matrix_trace(m);
	    break;
	case F_NORM1:
	    ret->v.xval = gretl_matrix_one_norm(m);
	    break;
	case F_INFNORM:
	    ret->v.xval = gretl_matrix_infinity_norm(m);
	    break;
	case F_RCOND:
	    ret->v.xval = gretl_matrix_rcond(m, &p->err);
	    break;
	case F_CNUMBER:
	    ret->v.xval = gretl_matrix_cond_index(m, &p->err);
	    break;
	case F_RANK:
	    ret->v.xval = gretl_matrix_rank(m, &p->err);
	    break;
	default:
	    p->err = E_PARSE;
	    break;
	}

	if (n->t == NUM) {
	    gretl_matrix_free(m);
	}

	if (p->err) {
	    matrix_error(p);
	}
    }

    return ret;
}

static NODE *matrix_add_names (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *m = l->v.m;
	int byrow = (f == F_ROWNAMES);

	if (r->t == STR) {
	    ret->v.xval = umatrix_set_names_from_string(m, r->v.str, byrow);
	} else if (r->t == ARRAY) {
	    if (gretl_array_get_type(r->v.a) != GRETL_TYPE_STRINGS) {
		p->err = E_TYPES;
	    } else {
		ret->v.xval = umatrix_set_names_from_array(m, r->v.a, byrow);
	    }
	} else {
	    /* some sort of list-bearing node */
	    int *list = node_get_list(r, p);

	    if (p->err) {
		ret->v.xval = 1;
	    } else {
		ret->v.xval = umatrix_set_names_from_list(m, list, p->dset,
							  byrow);
	    }
	    free(list);
	}
    }

    return ret;
}

static NODE *matrix_get_colname (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	int c = node_get_int(r, p);

	ret->v.str = user_matrix_get_column_name(l->v.m, c, &p->err);
    }

    return ret;
}

static NODE *matrix_imhof (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	const gretl_matrix *m = l->v.m;
	double arg = node_get_scalar(r, p);

	ret->v.xval = imhof(m, arg, &p->err);
    }

    return ret;
}

static void matrix_minmax_indices (int f, int *mm, int *rc, int *idx)
{
    *mm = (f == F_MAXR || f == F_MAXC || f == F_IMAXR || f == F_IMAXC);
    *rc = (f == F_MINC || f == F_MAXC || f == F_IMINC || f == F_IMAXC);
    *idx = (f == F_IMINR || f == F_IMINC || f == F_IMAXR || f == F_IMAXC);
}

#define cmplx_func(f) (f == HF_CMATRIX || f == HF_CMMULT || \
		       f == HF_CINV || f == HF_CFFT || \
		       f == HF_CTRAN || f == HF_CHPROD)

static NODE *matrix_to_matrix_func (NODE *n, NODE *r, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *m = NULL;
	int tmpmat = 0;
	int a = 0, b, c;

	if (n->t == NUM) {
	    if (cmplx_func(f)) {
		p->err = E_INVARG;
	    } else {
		m = gretl_matrix_from_scalar(node_get_scalar(n, p));
		tmpmat = 1;
	    }
	} else {
	    /* the following may be dangerous? */
	    tmpmat = is_tmp_node(n);
	    m = n->v.m;
	}

	if (p->err) {
	    goto finalize;
	}

	if (gretl_is_null_matrix(m) && !nullmat_ok(f)) {
	    p->err = E_DATA;
	} else if (f == F_RESAMPLE || f == F_MREVERSE || f == F_SDC) {
	    /* the r node may be absent, but if present it should
	       hold a scalar */
	    if (!empty_or_num(r)) {
		node_type_error(f, 2, NUM, r, p);
	    }
	} else if (f == HF_CFFT) {
	    /* optional scalar (boolean) on the right */
	    a = node_get_bool(r, p, 0);
	} else if (f == F_RANKING) {
	    if (gretl_vector_get_length(m) == 0) {
		/* m must be a vector */
		p->err = E_TYPES;
	    }
	}

	if (p->err) {
	    goto finalize;
	}

	gretl_error_clear();

	switch (f) {
	case F_SUMC:
	    ret->v.m = gretl_matrix_column_sum(m, &p->err);
	    break;
	case F_SUMR:
	    ret->v.m = gretl_matrix_row_sum(m, &p->err);
	    break;
	case F_PRODC:
	    ret->v.m = gretl_matrix_column_prod(m, &p->err);
	    break;
	case F_PRODR:
	    ret->v.m = gretl_matrix_row_prod(m, &p->err);
	    break;
	case F_MEANC:
	    ret->v.m = gretl_matrix_column_mean(m, &p->err);
	    break;
	case F_MEANR:
	    ret->v.m = gretl_matrix_row_mean(m, &p->err);
	    break;
	case F_SD:
	    ret->v.m = gretl_matrix_column_sd(m, &p->err);
	    break;
	case F_SDC:
	    if (r != NULL && r->t == NUM) {
		ret->v.m = gretl_matrix_column_sd2(m, r->v.xval, &p->err);
	    } else {
		ret->v.m = gretl_matrix_column_sd(m, &p->err);
	    }
	    break;
	case F_MCOV:
	    ret->v.m = gretl_covariance_matrix(m, 0, &p->err);
	    break;
	case F_MCORR:
	    ret->v.m = gretl_covariance_matrix(m, 1, &p->err);
	    break;
	case F_CUM:
	    ret->v.m = gretl_matrix_cumcol(m, &p->err);
	    break;
	case F_DIFF:
	    ret->v.m = gretl_matrix_diffcol(m, 0, &p->err);
	    break;
	case F_DATAOK:
	    ret->v.m = gretl_matrix_isfinite(m, &p->err);
	    break;
	case F_RESAMPLE:
	    if (r != NULL && r->t == NUM) {
		ret->v.m = gretl_matrix_block_resample(m, r->v.xval, &p->err);
	    } else {
		ret->v.m = gretl_matrix_resample(m, &p->err);
	    }
	    break;
	case F_CDEMEAN:
	case F_CHOL:
	case F_PSDROOT:
	case F_INV:
	case F_INVPD:
	case F_GINV:
	case F_UPPER:
	case F_LOWER:
	    ret->v.m = user_matrix_matrix_func(m, tmpmat, f, &p->err);
	    break;
	case F_DIAG:
	    ret->v.m = gretl_matrix_get_diagonal(m, &p->err);
	    break;
	case F_TRANSP:
	    ret->v.m = gretl_matrix_copy_transpose(m);
	    break;
	case F_VEC:
	    ret->v.m = user_matrix_vec(m, &p->err);
	    break;
	case F_VECH:
	    ret->v.m = user_matrix_vech(m, &p->err);
	    break;
	case F_UNVECH:
	    ret->v.m = user_matrix_unvech(m, &p->err);
	    break;
	case F_MREVERSE:
	    if (r != NULL && r->t == NUM && r->v.xval != 0) {
		ret->v.m = gretl_matrix_reverse_cols(m);
	    } else {
		ret->v.m = gretl_matrix_reverse_rows(m);
	    }
	    break;
	case F_NULLSPC:
	    ret->v.m = gretl_matrix_right_nullspace(m, &p->err);
	    break;
	case F_MEXP:
	    ret->v.m = gretl_matrix_exp(m, &p->err);
	    break;
	case F_FFT:
	    ret->v.m = gretl_matrix_fft(m, &p->err);
	    break;
	case F_FFTI:
	    ret->v.m = gretl_matrix_ffti(m, &p->err);
	    break;
	case F_POLROOTS:
	    ret->v.m = gretl_matrix_polroots(m, &p->err);
	    break;
	case F_RANKING:
	    ret->v.m = rank_vector(m, F_SORT, &p->err);
	    break;
	case F_MINC:
	case F_MAXC:
	case F_MINR:
	case F_MAXR:
	case F_IMINC:
	case F_IMAXC:
	case F_IMINR:
	case F_IMAXR:
	    matrix_minmax_indices(f, &a, &b, &c);
	    ret->v.m = gretl_matrix_minmax(m, a, b, c, &p->err);
	    break;
	case HF_CMATRIX:
	    ret->v.m = gretl_cmatrix(m, r->v.m, &p->err);
	    break;
	case HF_CMMULT:
	    ret->v.m = gretl_zgemm(m, r->v.m, &p->err);
	    break;
	case HF_CHPROD:
	    ret->v.m = gretl_complex_hprod(m, r->v.m, &p->err);
	    break;
	case HF_CINV:
	    ret->v.m = gretl_zgetri(m, &p->err);
	    break;
	case HF_CFFT:
	    ret->v.m = gretl_complex_fft(m, a, &p->err);
	    break;
	case HF_CTRAN:
	    ret->v.m = gretl_ctran(m, &p->err);
	    break;
	default:
	    break;
	}

	if (ret->v.m == m) {
	    /* input matrix was recycled: avoid double-freeing */
	    if (n->t == MAT) {
		n->v.m = NULL;
	    } else {
		m = NULL;
	    }
	}

    finalize:

	if (n->t == NUM) {
	    /* trash the on-the-fly scalar matrix */
	    gretl_matrix_free(m);
	}

	if (ret->v.m == NULL) {
	    matrix_error(p);
	}
    }

    return ret;
}

static NODE *read_object_func (NODE *n, NODE *r, int f, parser *p)
{
    NODE *ret;

    if (f == F_MREAD) {
	ret = aux_matrix_node(p);
    } else {
	ret = aux_bundle_node(p);
    }

    if (ret != NULL && starting(p)) {
	int import = 0;

	if (!null_or_empty(r)) {
	    import = (r->v.xval != 0);
	}

	gretl_error_clear();

	switch (f) {
	case F_MREAD:
	    ret->v.m = gretl_matrix_read_from_file(n->v.str, import,
						   &p->err);
	    break;
	case F_BREAD:
	    ret->v.b = gretl_bundle_read_from_file(n->v.str, import,
						   &p->err);
	    break;
	default:
	    break;
	}

#if 0
	/* Allin, 2016-08-03: I really don't like this
	   behaviour! */
	if (f == F_MREAD && p->err == E_FOPEN) {
	    ret->v.m = gretl_null_matrix_new();
	    p->err = 0;
	}
#endif

	if (f == F_MREAD && ret->v.m == NULL) {
	    matrix_error(p);
	}
    }

    return ret;
}

static NODE *
matrix_to_matrix2_func (NODE *n, NODE *r, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const gretl_matrix *m = n->v.m;
	const char *rname;

	if (gretl_is_null_matrix(m)) {
	    p->err = E_DATA;
	    goto finalize;
	}

	gretl_error_clear();

	/* on the right: address of matrix or null */
	if (r->t == EMPTY) {
	    rname = "null";
	} else {
	    /* note: switch to the 'content' sub-node */
	    r = r->v.b1.b;
	    if (umatrix_node(r)) {
		rname = r->vname;
	    } else {
		p->err = E_PARSE;
		gretl_errmsg_set("Expected the address of a matrix");
		return ret;
	    }
	}

	switch (f) {
	case F_QR:
	    ret->v.m = user_matrix_QR_decomp(m, rname, &p->err);
	    break;
	case F_EIGSYM:
	    ret->v.m = user_matrix_eigen_analysis(m, rname, 1, &p->err);
	    break;
	case F_EIGGEN:
	    ret->v.m = user_matrix_eigen_analysis(m, rname, 0, &p->err);
	    break;
	case HF_CEIGH:
	    {
		gretl_matrix *V = NULL;

		if (r->t != EMPTY) {
		    V = r->v.m;
		}
		ret->v.m = gretl_zheev(m, V, &p->err);
	    }
	    break;
	}

    finalize:

	if (ret->v.m == NULL) {
	    matrix_error(p);
	}
    }

    return ret;
}

static int ok_matrix_dim (int r, int c, int f)
{
    if (f == F_IMAT || f == F_ZEROS || f == F_ONES || f == F_MUNIF || \
	f == F_MNORM) {
	/* zero is OK for matrix creation functions, which then
	   return an empty matrix
	*/
	return (r >= 0 && c >= 0);
    } else {
	double xm = (double) r * (double) c;

	return (r > 0 && c > 0 && xm < INT_MAX);
    }
}

static NODE *matrix_fill_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	int cols = 0, rows = node_get_int(l, p);

	if (!p->err) {
	    if (f == F_IMAT) {
		cols = rows;
	    } else {
		cols = node_get_int(r, p);
	    }
	}

	if (!p->err && !ok_matrix_dim(rows, cols, f)) {
	    p->err = E_INVARG;
	    matrix_error(p);
	}

	if (p->err) {
	    return ret;
	}

	switch (f) {
	case F_IMAT:
	    ret->v.m = gretl_identity_matrix_new(rows);
	    break;
	case F_ZEROS:
	    ret->v.m = gretl_zero_matrix_new(rows, cols);
	    break;
	case F_ONES:
	    ret->v.m = gretl_unit_matrix_new(rows, cols);
	    break;
	case F_MUNIF:
	    ret->v.m = gretl_random_matrix_new(rows, cols,
					       D_UNIFORM);
	    break;
	case F_MNORM:
	    ret->v.m = gretl_random_matrix_new(rows, cols,
					       D_NORMAL);
	    break;
	default:
	    break;
	}
    }

    return ret;
}

#if EDEBUG > 1

static void print_mspec (matrix_subspec *mspec)
{
    const char *mstypes[] = {
	"SEL_RANGE",
	"SEL_ELEMENT",
	"SEL_MATRIX",
	"SEL_DIAG",
	"SEL_ALL",
	"SEL_NULL"
    };

    fprintf(stderr, "mspec at %p:\n", (void *) mspec);

    if (mspec != NULL) {
	int i;

	for (i=0; i<2; i++) {
	    fprintf(stderr, "type[%d] = %s\n", i, mstypes[mspec->type[i]]);
	    if (mspec->type[i] == SEL_RANGE) {
		fprintf(stderr, "sel[%d].range[0] = %d\n",
			i, mspec->sel[i].range[0]);
		fprintf(stderr, "sel[%d].range[1] = %d\n",
			i, mspec->sel[i].range[1]);
	    } else if (mspec->type[i] == SEL_MATRIX) {
		gretl_matrix_print(mspec->sel[i].m, "sel matrix");
	    }
	}
    }
}

#endif

/* compose a sub-matrix specification, from scalars and/or
   index matrices */

static void build_mspec (NODE *targ, NODE *l, NODE *r, parser *p)
{
    matrix_subspec *spec = targ->v.mspec;

    if (spec == NULL) {
	spec = matrix_subspec_new();
	if (spec == NULL) {
	    p->err = E_ALLOC;
	    return;
	}
    }

#if EDEBUG > 1
    fprintf(stderr, "build_mspec: l->t=%d (%s)\n", l->t, getsymb(l->t));
    if (r == NULL) {
	fprintf(stderr, " r = NULL\n");
    } else {
	fprintf(stderr, " r->t=%d (%s)\n", r->t, getsymb(r->t));
    }
#endif

    if (l->t == DUM) {
	if (l->v.idnum == DUM_DIAG) {
	    spec->type[0] = SEL_DIAG;
	    spec->type[1] = SEL_ALL;
	} else {
	    p->err = E_TYPES;
	}
	goto finished;
    } else if (l->t == NUM && r != NULL && r->t == NUM) {
	spec->type[0] = spec->type[1] = SEL_ELEMENT;
	mspec_set_row_index(spec, l->v.xval);
	mspec_set_col_index(spec, r->v.xval);
	goto finished;
    }

    if (l->t == NUM) {
	spec->type[0] = SEL_RANGE;
	mspec_set_row_index(spec, l->v.xval);
    } else if (l->t == IVEC) {
	spec->type[0] = SEL_RANGE;
	spec->sel[0].range[0] = l->v.ivec[0];
	spec->sel[0].range[1] = l->v.ivec[1];
    } else if (l->t == MAT) {
	if (gretl_matrix_is_scalar(l->v.m)) {
	    spec->type[0] = SEL_RANGE;
	    mspec_set_row_index(spec, l->v.m->val[0]);
	} else {
	    spec->type[0] = SEL_MATRIX;
	    spec->sel[0].m = l->v.m;
	}
    } else if (l->t == EMPTY) {
	spec->type[0] = SEL_ALL;
    } else {
	p->err = E_TYPES;
	goto finished;
    }

    if (r == NULL) {
	spec->type[1] = SEL_NULL;
    } else if (r->t == NUM) {
	spec->type[1] = SEL_RANGE;
	mspec_set_col_index(spec, r->v.xval);
    } else if (r->t == IVEC) {
	spec->type[1] = SEL_RANGE;
	spec->sel[1].range[0] = r->v.ivec[0];
	spec->sel[1].range[1] = r->v.ivec[1];
    } else if (r->t == MAT) {
	if (gretl_matrix_is_scalar(r->v.m)) {
	    spec->type[1] = SEL_RANGE;
	    mspec_set_col_index(spec, r->v.m->val[0]);
	} else {
	    spec->type[1] = SEL_MATRIX;
	    spec->sel[1].m = r->v.m;
	}
    } else if (r->t == EMPTY) {
	spec->type[1] = SEL_ALL;
    } else {
	p->err = E_TYPES;
    }

 finished:

#if EDEBUG > 1
    print_mspec(spec);
#endif

    if (p->err && spec != NULL) {
	free(spec);
	spec = NULL;
    }

    targ->v.mspec = spec;
}

/* node holding evaluated result of matrix specification */

static NODE *mspec_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_mspec_node(p);

    if (ret != NULL && starting(p)) {
	build_mspec(ret, l, r, p);
    }

    return ret;
}

static NODE *submatrix_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *a = NULL;
	matrix_subspec *spec = NULL;

	if (r->t != MSPEC) {
	    fprintf(stderr, "submatrix_node: couldn't find mspec\n");
	    p->err = E_TYPES;
	    return NULL;
	}

	spec = r->v.mspec;

	if (l->t == MAT) {
	    a = matrix_get_submatrix(l->v.m, spec, 0, &p->err);
	} else if (l->t == STR) {
	    /* FIXME never reached? */
	    a = user_matrix_get_submatrix(l->v.str, spec, &p->err);
	} else {
	    p->err = E_TYPES;
	}

	if (a != NULL) {
	    ret = aux_matrix_node(p);
	    if (ret == NULL) {
		gretl_matrix_free(a);
	    } else {
		ret->v.m = a;
	    }
	}
    } else {
	ret = aux_matrix_node(p);
    }

    return ret;
}

static int stored_list_check (const int *list, const DATASET *dset)
{
    int badv = 0;
    int err = 0;

    if (dset == NULL || dset->n == 0) {
	err = E_NODATA;
    } else {
	int i;

	for (i=1; i<=list[0] && !err; i++) {
	    if (list[i] >= dset->v ||
		(list[i] < 0 && list[i] != LISTSEP)) {
		badv = list[i];
		err = E_DATA;
	    }
	}
    }

    if (badv != 0) {
	gretl_errmsg_sprintf("list check: series ID %d "
			     "is out of bounds", badv);
    }

    return err;
}

static NODE *array_element_node (gretl_array *a, int i,
				 parser *p)
{
    NODE *ret = NULL;
    GretlType type = 0;
    void *data;

    data = gretl_array_get_element(a, i-1, &type, &p->err);

    if (p->err == E_INVARG) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i);
    }

    if (!p->err) {
	if (type == GRETL_TYPE_STRING) {
	    /* revised 2017-05-21 */
	    ret = string_pointer_node(p);
	    if (ret != NULL) {
		ret->v.str = data;
	    }
	} else if (type == GRETL_TYPE_MATRIX) {
	    ret = matrix_pointer_node(p);
	    if (ret != NULL) {
		ret->v.m = data;
	    }
	} else if (type == GRETL_TYPE_BUNDLE) {
	    ret = bundle_pointer_node(p);
	    if (ret != NULL) {
		ret->v.b = data;
	    }
	} else if (type == GRETL_TYPE_LIST) {
	    /* revised 2017-05-21 */
	    p->err = stored_list_check((const int *) data, p->dset);
	    if (!p->err) {
		ret = list_pointer_node(p);
		if (ret != NULL) {
		    ret->v.ivec = data;
		}
	    }
	}
    }

    return ret;
}

static int check_range_spec (int range[], int n)
{
    int err = 0;

    if (range[0] < 1 || range[0] > n) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), range[0]);
	err = E_DATA;
    } else if (range[1] == -999) {
	range[1] = n;
    }

    if (!err && (range[1] < 1 || range[1] > n || range[1] < range[0])) {
	err = E_DATA;
    }

    return err;
}

static NODE *array_range_node (gretl_array *a, int range[],
			       parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	p->err = check_range_spec(range, gretl_array_get_length(a));
	if (!p->err) {
	    ret = aux_array_node(p);
	    if (ret != NULL) {
		ret->v.a = gretl_array_copy_range(a, range[0] - 1,
						  range[1] - 1,
						  &p->err);
	    }
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *list_range_node (int *list, int range[], parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	p->err = check_range_spec(range, list[0]);
	if (!p->err) {
	    int n = range[1] - range[0] + 1;

	    ret = aux_list_node(p);
	    if (ret != NULL) {
		ret->v.ivec = gretl_list_new(n);
		if (ret->v.ivec != NULL) {
		    int i, j = 1;

		    for (i=range[0]; i<=range[1]; i++) {
			ret->v.ivec[j++] = list[i];
		    }
		}
	    }
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *string_range_node (const char *s, int r1, int r2, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	int range[2] = {r1, r2};

	p->err = check_range_spec(range, strlen(s));
	if (!p->err) {
	    int n = r2 - r1 + 1;

	    ret = aux_string_node(p);
	    if (ret != NULL) {
		ret->v.str = gretl_strndup(s + r1 - 1, n);
		if (ret->v.str == NULL) {
		    p->err = E_ALLOC;
		}
	    }
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *real_list_series_node (int *list, int i, parser *p)
{
    NODE *ret = NULL;
    int v = 0;

    if (i < 1 || i > list[0]) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i);
	p->err = E_INVARG;
    } else {
	v = list[i];
	if (v < 0 || v >= p->dset->v) {
	    gretl_errmsg_sprintf(_("Variable number %d is out of bounds"), v);
	    p->err = E_DATA;
	}
    }

    if (!p->err) {
	ret = aux_empty_series_node(p);
	if (!p->err) {
	    /* scrub TMP_NODE, because using dset->Z member! */
	    ret->flags = AUX_NODE;
	    ret->vnum = v;
	    ret->v.xvec = p->dset->Z[v];
	}
    }

    return ret;
}

/* coming from a context where we have @list and @i */

static NODE *list_member_node (int *list, int i, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	ret = real_list_series_node(list, i, p);
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static int mspec_get_series_index (matrix_subspec *s,
				   parser *p)
{
    int t = -1;

    if (s->type[0] == SEL_RANGE &&
	s->type[1] == SEL_NULL) {
	if (s->sel[0].range[0] == s->sel[0].range[1]) {
	    t = s->sel[0].range[0];
	} else {
	    /* allow for dates such as "2008:4" */
	    gchar *tmp;

	    tmp = g_strdup_printf("%d:%d", s->sel[0].range[0],
				  s->sel[0].range[1]);
	    t = get_observation_number(tmp, p->dset);
	    g_free(tmp);
	}
    }

    if (t < 1 || t > p->dset->n) {
	p->err = E_DATA;
    }

    return t;
}

/* Here we're checking @spec for suitability in specifying
   a "slice" of an array, list or string: unlike the matrix
   case, @spec must be one-dimensional, and moreover (at
   present) it must come down to a single element or a
   straightfoward range (so not a vector).
*/

static int test_for_single_range (matrix_subspec *spec,
				  parser *p)
{
    int ret = 0;

    if (spec->type[0] == SEL_RANGE &&
	spec->type[1] == SEL_NULL) {
	if (spec->sel[0].range[0] == spec->sel[0].range[1]) {
	    ret = spec->sel[0].range[0];
	}
    } else {
	p->err = E_TYPES;
	ret = -1;
    }

    return ret;
}

static NODE *subobject_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	if (l->t == MAT && r->t == MSPEC) {
	    return submatrix_node(l, r, p);
	} else if ((l->t == ARRAY || l->t == LIST || l->t == STR) &&
		   r->t == MSPEC) {
	    int i = test_for_single_range(r->v.mspec, p);

	    if (!p->err && i > 0) {
		if (l->t == ARRAY) {
		    ret = array_element_node(l->v.a, i, p);
		} else if (l->t == LIST) {
		    ret = list_member_node(l->v.ivec, i, p);
		} else {
		    ret = string_range_node(l->v.str, i, i, p);
		}
	    } else if (!p->err) {
		int range[2];

		range[0] = r->v.mspec->sel[0].range[0];
		range[1] = r->v.mspec->sel[0].range[1];
		if (l->t == ARRAY) {
		    ret = array_range_node(l->v.a, range, p);
		} else if (l->t == LIST) {
		    ret = list_range_node(l->v.ivec, range, p);
		} else {
		    ret = string_range_node(l->v.str, range[0], range[1], p);
		}
	    }
	} else if (l->t == SERIES && r->t == MSPEC) {
	    int t = mspec_get_series_index(r->v.mspec, p);

	    if (!p->err) {
		ret = aux_scalar_node(p);
		if (!p->err) {
		    ret->v.xval = l->v.xvec[t-1];
		}
	    }
	} else {
	    fprintf(stderr, "subobject_node: l='%s', r='%s'\n",
		    getsymb(l->t), getsymb(r->t));
	    p->err = E_TYPES;
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *process_subslice (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	if (scalar_node(l) && (scalar_node(r) || r->t == EMPTY)) {
	    ret = aux_ivec_node(p, 2);
	    if (ret != NULL) {
		ret->v.ivec[0] = node_get_int(l, p);
		ret->v.ivec[1] = (r->t == EMPTY)?
		    MSEL_MAX : node_get_int(r, p);
	    }
	} else {
	    p->err = E_TYPES;
	}
    } else {
	ret = aux_ivec_node(p, 2);
    }

    return ret;
}

static double real_apply_func (double x, int f, parser *p)
{
    double y;

    errno = 0;

    if (na(x)) {
	switch (f) {
	case F_MISSING:
	    return 1.0;
	case F_DATAOK:
	case F_MISSZERO:
	    return 0.0;
	default:
	    if (na(x)) {
		return NADBL;
	    }
	}
    }

    switch (f) {
    case U_NEG:
	return -x;
    case U_POS:
	return x;
    case U_NOT:
	return x == 0;
    case F_ABS:
	return fabs(x);
    case F_TOINT:
	return (double) (int) x;
    case F_CEIL:
	return ceil(x);
    case F_FLOOR:
	return floor(x);
    case F_ROUND:
	return gretl_round(x);
    case F_SIN:
	return sin(x);
    case F_COS:
	return cos(x);
    case F_TAN:
	return tan(x);
    case F_ASIN:
	return asin(x);
    case F_ACOS:
	return acos(x);
    case F_ATAN:
	return atan(x);
    case F_SINH:
	return sinh(x);
    case F_COSH:
	return cosh(x);
    case F_TANH:
	return tanh(x);
    case F_ASINH:
	return asinh(x);
    case F_ACOSH:
	return acosh(x);
    case F_ATANH:
	return atanh(x);
    case F_CNORM:
	return normal_cdf(x);
    case F_DNORM:
	return normal_pdf(x);
    case F_QNORM:
	return normal_cdf_inverse(x);
    case F_LOGISTIC:
	return logistic_cdf(x);
    case F_GAMMA:
	y = gamma_function(x);
	if (na(y)) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_LNGAMMA:
	y = ln_gamma(x);
	if (na(y)) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_DIGAMMA:
	y = digamma(x);
	if (na(y)) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_MISSING:
	return 0.0;
    case F_DATAOK:
	return 1.0;
    case F_MISSZERO:
	return x;
    case F_ZEROMISS:
	return (x == 0.0)? NADBL : x;
    case F_SQRT:
	y = sqrt(x);
	if (errno) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_LOG:
	y = log(x);
	if (errno) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_LOG10:
	y = log10(x);
	if (errno) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_LOG2:
	y = log2(x);
	if (errno) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_EXP:
	y = exp(x);
	if (errno) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_INVMILLS:
	y = invmills(x);
	if (na(y)) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_EASTER:
	y = easterdate(x);
	return y;
    default:
	return 0.0;
    }
}

static NODE *apply_scalar_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	ret->v.xval = real_apply_func(n->v.xval, f, p);
    }

    return ret;
}

static NODE *sleep_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	int s = node_get_int(n, p);

	g_usleep(G_USEC_PER_SEC * s);
	ret->v.xval = 0;
    }

    return ret;
}

static NODE *scalar_isnan_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	double x = node_get_scalar(n, p);

	if (!p->err) {
	    ret->v.xval = isnan(x);
	}
    }

    return ret;
}

static NODE *matrix_isnan_node (NODE *n, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const gretl_matrix *m = n->v.m;

	if (m->rows == 0 || m->cols == 0) {
	    p->err = E_DATA;
	} else {
	    ret->v.m = gretl_matrix_alloc(m->rows, m->cols);
	    if (ret->v.m == NULL) {
		p->err = E_ALLOC;
	    } else {
		int i, n = m->rows * m->cols;

		for (i=0; i<n; i++) {
		    ret->v.m->val[i] = isnan(m->val[i]);
		}
	    }
	}
    }

    return ret;
}

static NODE *apply_series_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_series_node(p);
    int t;

#if 0 /* initial parallelization experiment */
    if (ret != NULL) {
	const double *x;

	if (n->t == SERIES) {
	    x = n->v.xvec;
	} else {
	    x = get_colvec_as_series(n, f, p);
	}

	if (!p->err) {
	    if (autoreg(p)) {
		ret->v.xvec[p->obs] = real_apply_func(x[p->obs], f, p);
	    } else {
#if defined(_OPENMP)
#pragma omp parallel for private(t)
		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    ret->v.xvec[t] = real_apply_func(x[t], f, p);
		}
#else
		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    ret->v.xvec[t] = real_apply_func(x[t], f, p);
		}
#endif
	    }
	}
    }
#else
    if (ret != NULL) {
	const double *x;

	if (n->t == SERIES) {
	    x = n->v.xvec;
	} else {
	    x = get_colvec_as_series(n, f, p);
	}

	if (!p->err) {
	    if (autoreg(p)) {
		ret->v.xvec[p->obs] = real_apply_func(x[p->obs], f, p);
	    } else {
		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    ret->v.xvec[t] = real_apply_func(x[t], f, p);
		}
	    }
	}
    }
#endif

    return ret;
}

/* argument is series or list; value returned is list */

static NODE *dummify_func (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	int *list = NULL;
	double oddval = NADBL;

	if (!null_or_empty(r)) {
	    if (r->t != NUM) {
		p->err = E_TYPES;
		return ret;
	    } else {
		oddval = r->v.xval;
	    }
	}

	if (l->t == LIST) {
	    list = gretl_list_copy(l->v.ivec);
	} else if (useries_node(l)) {
	    list = gretl_list_new(1);
	    list[1] = l->vnum;
	} else {
	    p->err = E_TYPES;
	}

	if (p->err) {
	    ; /* don't do anything more */
	} else if (list == NULL) {
	    p->err = E_ALLOC;
	} else if (r->t == EMPTY) {
	    /* got just one argument */
	    p->err = list_dumgenr(&list, p->dset, OPT_F);
	    ret->v.ivec = list;
	} else if (list[0] > 1) {
	    gretl_errmsg_set("dummify(x, y): first argument should be a single variable");
	    free(list);
	    p->err = E_DATA;
	} else {
	    p->err = dumgenr_with_oddval(&list, p->dset, oddval);
	    ret->v.ivec = list;
	}
    }

    return ret;
}

/* argument is list; value returned is list */

static NODE *cdummify_func (NODE *n, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	int *list = NULL;

	if (n->t == LIST) {
	    list = gretl_list_copy(n->v.ivec);
	} else if (useries_node(n)) {
	    list = gretl_list_new(1);
	    list[1] = n->vnum;
	} else {
	    p->err = E_TYPES;
	}

	if (p->err) {
	    ; /* don't do anything more */
	} else if (list == NULL) {
	    p->err = E_ALLOC;
	} else {
	    p->err = auto_dummify_list(&list, p->dset);
	    ret->v.ivec = list;
	}
    }

    return ret;
}

static NODE *get_series_info (NODE *n, parser *p)
{
    NODE *ret = aux_bundle_node(p);

    if (ret != NULL && starting(p)) {
	int v = 0;

	if (useries_node(n)) {
	    v = n->vnum;
	} else {
	    v = node_get_int(n, p);
	}

	if (!p->err) {
	    ret->v.b = series_info_bundle(p->dset, v, &p->err);
	}
    }

    return ret;
}

/* middle argument is series or list; value returned is list in
   either case */

static NODE *list_make_lags (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	gretlopt opt = OPT_NONE;
	gretl_matrix *v = NULL;
	int *list = NULL;
	int k = 0;

	/* ordering of the results? */
	if (node_get_bool(r, p, 0) > 0 && !p->err) {
	    opt = OPT_L; /* by lags */
	}

	if (!p->err) {
	    /* scalar lag order or vector */
	    if (l->t == NUM) {
		k = node_get_int(l, p);
	    } else {
		v = l->v.m;
	    }
	}

	if (!p->err) {
	    list = node_get_list(m, p);
	}

	if (!p->err) {
	    if (list[0] > 0) {
		p->err = list_laggenr(&list, 1, k, v,
				      p->dset, 0, opt);
	    }
	    ret->v.ivec = list;
	}
    }

    return ret;
}

static NODE *matrix_make_lags (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	gretlopt opt = OPT_NONE;
	gretl_matrix *kvec = NULL;
	gretl_matrix *src = NULL;

	if (node_get_bool(r, p, 0) > 0 && !p->err) {
	    opt = OPT_L; /* by lags */
	}

	if (!p->err) {
	    /* scalar max lag order or vector */
	    if (l->t == NUM) {
		int i, k = node_get_int(l, p);

		if (!p->err && k <= 0) {
		    /* scalar k must be positive */
		    p->err = E_INVARG;
		} else if (!p->err) {
		    kvec = gretl_vector_alloc(k);
		    if (kvec == NULL) {
			p->err = E_ALLOC;
		    } else {
			for (i=0; i<k; i++) {
			    kvec->val[i] = i+1;
			}
		    }
		}
	    } else {
		kvec = l->v.m;
	    }
	}

	if (!p->err) {
	    src = m->v.m;
	}

	if (!p->err) {
	    ret->v.m = gretl_matrix_lag(src, kvec, opt, 0.0);
	}

	if (kvec != l->v.m) {
	    gretl_matrix_free(kvec);
	}
    }

    return ret;
}

/* args are minlag, maxlag, MIDAS-list-to-be-lagged */

static NODE *hf_list_make_lags (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	int *list = NULL;
	int lmin, lmax = 0;
	int cfac = 0;

	lmin = node_get_int(l, p);
	if (!p->err) {
	    lmax = node_get_int(m, p);
	}
	if (!p->err) {
	    list = node_get_list(r, p);
	}

	if (!p->err) {
	    /* compaction factor for high-frequency data */
	    cfac = list[0];
	    if (cfac < 2) {
		fprintf(stderr, "hflags: not a MIDAS list\n");
		p->err = E_INVARG;
		free(list);
	    }
	}

	if (!p->err) {
	    if (list[0] > 0) {
		p->err = list_laggenr(&list, lmin, lmax, NULL,
				      p->dset, cfac, OPT_L);
	    }
	    ret->v.ivec = list;
	}
    }

    return ret;
}

#define ok_list_func(f) (f == F_LOG || f == F_DIFF || \
			 f == F_LDIFF || f == F_SDIFF || \
			 f == F_SQUARE || f == F_ODEV || \
			 f == F_RESAMPLE || f == F_DROPCOLL || \
			 f == F_HFDIFF || f == F_HFLDIFF)

/* The following handles functions that are "basically" for series,
   but which can also be applied to lists -- except for F_DROPCOLL,
   F_HFDIFF and F_HDLDIFF, which require a list argument.
*/

static NODE *apply_list_func (NODE *n, NODE *r, int f, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (!ok_list_func(f)) {
	p->err = E_TYPES;
	return ret;
    }

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(n, p);
	gretlopt opt = OPT_NONE;
	double parm = NADBL;
	int t = 0;

	if (f == F_SQUARE) {
	    if (r != NULL && node_is_true(r, p)) {
		opt = OPT_O;
	    }
	} else if (f == F_DROPCOLL || f == F_HFDIFF ||
		   f == F_HFLDIFF) {
	    /* handle optional parameter */
	    if (!null_or_empty(r)) {
		parm = node_get_scalar(r, p);
		if (p->err) {
		    return ret;
		}
	    }
	}

	/* note: @list is modified by the library functions
	   called below */

	if (list != NULL) {
	    /* note: an empty list argument produces an
	       empty list return
	    */
	    if (list[0] > 0) {
		switch (f) {
		case F_LOG:
		    p->err = list_loggenr(list, p->dset);
		    break;
		case F_DIFF:
		case F_LDIFF:
		case F_SDIFF:
		    if (f == F_DIFF) t = DIFF;
		    else if (f == F_LDIFF) t = LDIFF;
		    else if (f == F_SDIFF) t = SDIFF;
		    p->err = list_diffgenr(list, t, p->dset);
		    break;
		case F_SQUARE:
		    p->err = list_xpxgenr(&list, p->dset, opt);
		    break;
		case F_ODEV:
		    p->err = list_orthdev(list, p->dset);
		    break;
		case F_RESAMPLE:
		    p->err = list_resample(list, p->dset);
		    break;
		case F_DROPCOLL:
		    p->err = list_dropcoll(list, parm, p->dset);
		    break;
		case F_HFDIFF:
		case F_HFLDIFF:
		    t = (f == F_HFDIFF)? DIFF : LDIFF;
		    p->err = hf_list_diffgenr(list, t, parm, p->dset);
		    break;
		default:
		    break;
		}
	    }
	    ret->v.ivec = list;
	}
    }

    return ret;
}

static NODE *hf_list_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    gretl_matrix *v = l->v.m;
    int f_ratio = node_get_int(m, p);
    char *pfx = r->v.str;
    NODE *ret = NULL;

    if (!p->err) {
	int n = gretl_vector_get_length(v);

	if (n == 0) {
	    p->err = E_NONCONF;
	} else if (f_ratio < 3) {
	    p->err = E_INVARG;
	} else if (*pfx == '\0' || !gretl_is_ascii(pfx) ||
		   strlen(pfx) > 24) {
	    p->err = E_INVARG;
	} else {
	    int T = sample_size(p->dset);

	    if (n != f_ratio * T) {
		p->err = E_INVARG;
	    }
	}
    }

    if (!p->err) {
	ret = aux_list_node(p);
	if (ret != NULL) {
	    ret->v.ivec = vector_to_midas_list(v, f_ratio, pfx,
					       p->dset, &p->err);
	}
    }

    return ret;
}

static NODE *dataset_list_node (parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	int *list = full_var_list(p->dset, NULL);

	if (list == NULL) {
	    list = gretl_null_list();
	}
	if (list == NULL) {
	    p->err = E_DATA;
	}
	ret->v.ivec = list;
    }

    return ret;
}

static NODE *trend_node (parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	ret = aux_empty_series_node(p);
	if (!p->err) {
	    p->err = gen_time(p->dset, 1, &ret->vnum);
	    if (!p->err) {
		ret->v.xvec = p->dset->Z[ret->vnum];
		/* not TMP_NODE because we're borrowing a Z column */
		ret->flags &= ~TMP_NODE;
	    }
	}
    }

    return ret;
}

static NODE *seasonals_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (!dataset_is_seasonal(p->dset) &&
	!dataset_is_seasonal_panel(p->dset)) {
	p->err = E_PDWRONG;
    } else {
	int ref = 0, center = 0;

	if (!null_or_empty(l)) {
	    ref = node_get_int(l, p);
	}
	if (!null_or_empty(r)) {
	    center = node_is_true(r, p);
	}

	if (!p->err) {
	    ret = aux_list_node(p);
	}

	if (ret != NULL) {
	    ret->v.ivec = seasonals_list(p->dset, ref, center, &p->err);
	}
    }

    return ret;
}

static NODE *get_lag_list (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	int *list = NULL, *srclist = NULL;
	int i, imin = 1, imax = 1;
	int lv = 0;

	if (!useries_node(l) && !ulist_node(l)) {
	    /* we need a named series or list on the left */
	    p->err = E_TYPES;
	} else if (r->t != IVEC && r->t != NUM) {
	    /* we need one or more integers on the right */
	    p->err = E_TYPES;
	}

	if (p->err) {
	    return NULL;
	}

	if (ulist_node(l)) {
	    srclist = l->v.ivec;
	    imax = srclist[0];
	} else {
	    lv = l->vnum;
	}

	if (imax == 0) {
	    /* empty list on input -> empty on output */
	    list = gretl_list_new(0);
	    if (list == NULL) {
		p->err = E_ALLOC;
	    }
	    goto loopdone;
	}

	for (i=imin; i<=imax && !p->err; i++) {
	    if (srclist != NULL) {
		lv = srclist[i];
	    }
	    if (r->t == IVEC) {
		int fromlag = r->v.ivec[0];
		int tolag = r->v.ivec[1];

		if (list == NULL) {
		    list = laggenr_from_to(lv, fromlag, tolag,
					   p->dset, &p->err);
		} else {
		    int *tmp;

		    tmp = laggenr_from_to(lv, fromlag, tolag,
					   p->dset, &p->err);
		    if (!p->err) {
			p->err = gretl_list_add_list(&list, tmp);
			free(tmp);
		    }
		}
	    } else {
		int lag = -r->v.xval;

		lv = laggenr(lv, lag, p->dset);
		if (lv > 0) {
		    list = gretl_list_append_term(&list, lv);
		    if (list == NULL) {
			p->err = E_ALLOC;
		    }
		}
	    }
	}

    loopdone:

	if (list != NULL) {
	    ret = aux_list_node(p);
	    if (ret != NULL) {
		ret->v.ivec = list;
	    } else {
		free(list);
	    }
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

int *list_from_strings_array (gretl_array *a, parser *p)
{
    GretlType type = gretl_array_get_type(a);
    int *list = NULL;

    if (type != GRETL_TYPE_STRINGS) {
	p->err = E_TYPES;
    } else {
	int i, vi, n = 0;
	char **S = gretl_array_get_strings(a, &n);

	for (i=0; i<n && !p->err; i++) {
	    vi = current_series_index(p->dset, S[i]);
	    if (vi < 0) {
		gretl_errmsg_sprintf("'%s' is not a known series", S[i]);
		p->err = E_UNKVAR;
	    }
	}

	if (!p->err) {
	    list = gretl_list_new(n);
	    if (list == NULL) {
		p->err = E_ALLOC;
	    } else {
		for (i=0; i<n; i++) {
		    list[i+1] = current_series_index(p->dset, S[i]);
		}
	    }
	}
    }

    return list;
}

/* get an *int list from node @n: note that the list is always
   newly allocated, and so should be freed by the caller if
   it's just for temporary use
*/

int *node_get_list (NODE *n, parser *p)
{
    int *list = NULL;
    int v = 0;

    if (n->t == LIST) {
	list = gretl_list_copy(n->v.ivec);
    } else if (n->t == SERIES || n->t == NUM) {
	v = (n->t == SERIES)? n->vnum : node_get_int(n, p);
	if (!p->err) {
	    if (v < 0 || v >= p->dset->v) {
		p->err = E_UNKVAR;
	    } else {
		list = gretl_list_new(1);
		if (list == NULL) {
		    p->err = E_ALLOC;
		} else {
		    list[1] = v;
		}
	    }
	}
    } else if (n->t == EMPTY) {
	list = gretl_null_list();
    } else if (dataset_dum(n)) {
	list = full_var_list(p->dset, NULL);
    } else if (n->t == MAT) {
	list = gretl_list_from_vector(n->v.m, p->dset, &p->err);
    } else {
	p->err = E_TYPES;
    }

    if (!p->err && list == NULL) {
	p->err = E_ALLOC;
    } else if (p->err == E_UNKVAR && v != 0) {
	gretl_errmsg_sprintf(_("Variable number %d is out of bounds"), v);
    }

    return list;
}

static NODE *eval_lcat (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	int *list1, *list2 = NULL;

	list1 = node_get_list(l, p); /* note, copied */
	if (list1 != NULL) {
	    list2 = node_get_list(r, p); /* copied */
	}
	if (list2 != NULL) {
	    p->err = gretl_list_add_list(&list1, list2);
	}
	ret->v.ivec = list1;
	free(list2);
    }

    return ret;
}

static NODE *list_list_op (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	int *llist, *rlist = NULL;
	int *list = NULL;

	llist = node_get_list(l, p);
	if (llist != NULL) {
	    rlist = node_get_list(r, p);
	}
	if (rlist != NULL) {
	    if (f == B_AND) {
		list = gretl_list_intersection(llist, rlist, &p->err);
	    } else if (f == B_OR) {
		list = gretl_list_union(llist, rlist, &p->err);
	    } else if (f == B_SUB) {
		list = gretl_list_drop(llist, rlist, &p->err);
	    } else if (f == B_POW) {
		list = gretl_list_product(llist, rlist, p->dset, &p->err);
	    } else if (f == B_ADD) {
		list = gretl_list_plus(llist, rlist, &p->err);
	    }
	}
	ret->v.ivec = list;
	free(llist);
	free(rlist);
    }

    return ret;
}

/* Binary operator applied to two bundles: at present only '+'
   (for union) is supported.
*/

static NODE *bundle_op (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_bundle_node(p);

    if (ret != NULL && starting(p)) {
	gretl_bundle *bl = l->v.b;
	gretl_bundle *br = r->v.b;

	if (!p->err) {
	    if (f == B_ADD) {
		ret->v.b = gretl_bundle_union(bl, br, &p->err);
	    } else {
		p->err = E_TYPES;
	    }
	}
    }

    return ret;
}

/* Binary operator applied to two arrays: at present only '+'
   (for appending) is supported.
*/

static NODE *array_op (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_array_node(p);

    if (ret != NULL && starting(p)) {
	gretl_array *al = l->v.a;
	gretl_array *ar = r->v.a;

	if (!p->err) {
	    if (f == B_ADD) {
		ret->v.a = gretl_arrays_join(al, ar, &p->err);
	    } else {
		p->err = E_TYPES;
	    }
	}
    }

    return ret;
}

/* in case we switched the LHS and RHS in a boolean comparison */

static int reversed_comp (int f)
{
    if (f == B_GT) {
	return B_LT;
    } else if (f == B_LT) {
	return B_GT;
    } else if (f == B_GTE) {
	return B_LTE;
    } else if (f == B_LTE) {
	return B_GTE;
    } else {
	return f;
    }
}

/* Boolean test of all vars in list against a scalar or series, for
   each observation in the sample, hence generating a series.
   The list will always be on the left-hand node; the 'reversed'
   flag is set if the list was originally on the right.
*/

static NODE *list_bool_comp (NODE *l, NODE *r, int f, int reversed,
			     parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(l, p);
	double *x = ret->v.xvec;
	double xit, targ = NADBL;
	double *tvec = NULL;
	int i, t;

	if (r->t == NUM) {
	    targ = r->v.xval;
	} else {
	    tvec = r->v.xvec;
	}

	if (reversed) {
	    f = reversed_comp(f);
	}

	if (list != NULL) {
	    for (t=p->dset->t1; t<=p->dset->t2; t++) {
		if (tvec != NULL) {
		    targ = tvec[t];
		}
		if (na(targ)) {
		    x[t] = NADBL;
		    continue;
		}
		x[t] = 1.0; /* assume 'true' */
		for (i=1; i<=list[0]; i++) {
		    xit = p->dset->Z[list[i]][t];
		    if (na(xit)) {
			x[t] = NADBL;
			break;
		    } else if (f == B_EQ && xit != targ) {
			x[t] = 0.0;
		    } else if (f == B_NEQ && xit == targ) {
			x[t] = 0.0;
		    } else if (f == B_LT && xit >= targ) {
			x[t] = 0.0;
		    } else if (f == B_GT && xit <= targ) {
			x[t] = 0.0;
		    } else if (f == B_LTE && xit > targ) {
			x[t] = 0.0;
		    } else if (f == B_GTE && xit < targ) {
			x[t] = 0.0;
		    }
		}
	    }
	    free(list);
	}
    }

    return ret;
}

/* Test for whether or not two lists are identical.  Note that
   using gretl_list_cmp() the order of the members matters.
   Perhaps the order shouldn't matter?
*/

static NODE *list_list_comp (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	int *llist = node_get_list(l, p);
	int *rlist = node_get_list(r, p);

	if (llist != NULL && rlist != NULL) {
	    int d = gretl_list_cmp(llist, rlist);

	    if (f == B_NEQ) {
		ret->v.xval = d;
	    } else if (f == B_EQ) {
		ret->v.xval = !d;
	    } else {
		p->err = E_TYPES;
	    }
	}
	free(llist);
	free(rlist);
    }

    return ret;
}

static NODE *num_string_comp (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	NODE *xnode = l->t == NUM ? l : r;
	NODE *snode = l->t == STR ? l : r;
	const char *s = snode->v.str;
	int v = xnode->vnum;

	if (v <= 0 || !is_string_valued(p->dset, v)) {
	    p->err = E_TYPES;
	} else {
	    double sx = series_decode_string(p->dset, v, s);

	    if (na(sx)) {
		ret->v.xval = NADBL;
	    } else if (f == B_EQ) {
		ret->v.xval = (sx == xnode->v.xval);
	    } else {
		ret->v.xval = (sx != xnode->v.xval);
	    }
	}
    }

    return ret;
}

/* argument is list; value returned is series */

static NODE *list_to_series_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(n, p);

	if (list != NULL) {
	    p->err = cross_sectional_stat(ret->v.xvec, list,
					  p->dset, f);
	    free(list);
	}
    }

    return ret;
}

/* arguments are series on left, list on right: we add all members
   of list to series, or subtract all members */

static NODE *series_list_calc (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(r, p);

	if (list != NULL) {
	    double xt, xi;
	    int i, t;

	    for (t=p->dset->t1; t<=p->dset->t2; t++) {
		xt = l->v.xvec[t];
		if (!na(xt)) {
		    for (i=1; i<=list[0]; i++) {
			xi = p->dset->Z[list[i]][t];
			if (na(xi)) {
			    xt = NADBL;
			    break;
			} else if (f == B_ADD) {
			    xt += xi;
			} else {
			    xt -= xi;
			}
		    }
		}
		ret->v.xvec[t] = xt;
	    }
	    free(list);
	}
    }

    return ret;
}

static int node_get_midas_method (NODE *n, parser *p)
{
    int ret = -1;

    if (n->t == STR) {
	if (!strcmp(n->v.str, "nealmon")) {
	    ret = 1;
	} else if (!strcmp(n->v.str, "beta0")) {
	    ret = 2;
	} else if (!strcmp(n->v.str, "betan")) {
	    ret = 3;
	} else if (!strcmp(n->v.str, "almonp")) {
	    ret = 4;
	} else {
	    p->err = E_INVARG;
	}
    } else {
	ret = node_get_int(n, p);
    }

    return ret;
}

static NODE *lincomb_func (NODE *l, NODE *m, NODE *r, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(l, p);
	const gretl_matrix *b = m->v.m;
	int k = 0;

	if (!p->err && (list == NULL || gretl_is_null_matrix(b))) {
	    p->err = E_DATA;
	}

	if (!p->err && f == F_MLINCOMB) {
	    k = node_get_midas_method(r, p);
	}

	if (!p->err) {
	    if (f == F_MLINCOMB) {
		p->err = midas_linear_combo(ret->v.xvec, list, b, k, p->dset);
	    } else {
		p->err = list_linear_combo(ret->v.xvec, list, b, p->dset);
	    }
	}

	free(list);
    }

    return ret;
}

static NODE *list_list_series_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
	int *llist = node_get_list(l, p);
	int *rlist = node_get_list(r, p);

	if (llist != NULL && rlist != NULL) {
	    p->err = x_sectional_weighted_stat(ret->v.xvec,
					       llist, rlist,
					       p->dset, f);
	}
	free(llist);
	free(rlist);
    }

    return ret;
}

/* check for missing obs in a list of variables */

static NODE *list_ok_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(n, p);
	int i, v, t;
	double x;

	if (list == NULL || list[0] == 0) {
	    free(list);
	    return ret;
	}

	for (t=p->dset->t1; t<=p->dset->t2; t++) {
	    x = (f == F_DATAOK)? 1 : 0;
	    for (i=1; i<=list[0]; i++) {
		v = list[i];
		if (na(p->dset->Z[v][t])) {
		    x = (f == F_DATAOK)? 0 : 1;
		    break;
		}
	    }
	    ret->v.xvec[t] = x;
	}

	free(list);
    }

    return ret;
}

/* functions taking (up to) two scalars as arguments and
   returning a series result */

static NODE *
series_fill_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
	double x, y;

	x = (l->t == EMPTY)? NADBL : node_get_scalar(l, p);
	y = (r->t == EMPTY)? NADBL : node_get_scalar(r, p);

	switch (f) {
	case F_RUNIFORM:
	    p->err = gretl_rand_uniform_minmax(ret->v.xvec,
					       p->dset->t1,
					       p->dset->t2,
					       x, y);
	    break;
	case F_RNORMAL:
	    p->err = gretl_rand_normal_full(ret->v.xvec,
					    p->dset->t1,
					    p->dset->t2,
					    x, y);
	    break;
	default:
	    break;
	}
    }

    return ret;
}

static gretl_matrix *fc_matrix_from_list (NODE *n, int n1,
					  parser *p)
{
    gretl_matrix *ret = NULL;

    if (sample_size(p->dset) != n1) {
	p->err = E_NONCONF;
    } else {
	ret = gretl_matrix_data_subset(n->v.ivec,
				       p->dset,
				       p->dset->t1,
				       p->dset->t2,
				       M_MISSING_OK,
				       &p->err);
    }

    return ret;
}

/* Functions taking two series as arguments and returning a scalar
   or matrix result. We also accept as arguments two matrices if
   they are vectors of the same length. In the case of F_NAALEN
   and F_KMEIER we can accept input with no right-hand argument
   (meaning no censoring).
*/

static NODE *series_2_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *Fmat = NULL;
	const double *x = NULL, *y = NULL;
	int free_Fmat = 0;
	int n = 0, n2 = 0;

	if (l->t == SERIES) {
	    /* series on left */
	    n = sample_size(p->dset);
	    x = l->v.xvec + p->dset->t1;
	} else {
	    /* must be matrix on left */
	    n = gretl_vector_get_length(l->v.m);
	    if (n == 0) {
		p->err = E_TYPES;
	    } else {
		x = l->v.m->val;
	    }
	}

	if (!p->err) {
	    if (null_or_empty(r)) {
		; /* OK for duration funcs */
	    } else if (r->t == SERIES) {
		/* series on right */
		n2 = sample_size(p->dset);
		if (n2 != n) {
		    p->err = E_NONCONF;
		} else {
		    y = r->v.xvec + p->dset->t1;
		}
	    } else if (r->t == MAT) {
		/* matrix on right */
		n2 = gretl_vector_get_length(r->v.m);
		if (n2 != n) {
		    if (f == F_FCSTATS && r->v.m->rows == n) {
			Fmat = r->v.m;
		    } else {
			p->err = E_NONCONF;
		    }
		} else {
		    y = r->v.m->val;
		}
	    } else if (r->t == LIST && f == F_FCSTATS) {
		/* a list may work for fcstats */
		Fmat = fc_matrix_from_list(r, n, p);
		if (!p->err) {
		    free_Fmat = 1;
		}
	    } else {
		p->err = E_TYPES;
	    }
	}

	if (p->err) {
	    return ret;
	}

	if (f == F_FCSTATS || f == F_NAALEN || f == F_KMEIER) {
	    ret = aux_matrix_node(p);
	} else {
	    ret = aux_scalar_node(p);
	}

	if (ret == NULL) {
	    return NULL;
	}

	/* n is taken as inclusive below */
	n--;

	switch (f) {
	case F_COR:
	    ret->v.xval = gretl_corr(0, n, x, y, NULL);
	    break;
	case F_COV:
	    ret->v.xval = gretl_covar(0, n, x, y, NULL);
	    break;
	case F_FCSTATS:
	    if (Fmat != NULL) {
		ret->v.m = matrix_fc_stats(x, Fmat, OPT_D, &p->err);
	    } else {
		ret->v.m = forecast_stats(x, y, 0, n, OPT_D, &p->err);
	    }
	    break;
	case F_NAALEN:
	    ret->v.m = duration_func(x, y, 0, n, OPT_NONE, &p->err);
	    break;
	case F_KMEIER:
	    ret->v.m = duration_func(x, y, 0, n, OPT_K, &p->err);
	    break;
	default:
	    break;
	}

	if (free_Fmat) {
	    gretl_matrix_free(Fmat);
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static gretlopt get_npcorr_option (NODE *n, parser *p)
{
    gretlopt opt = OPT_NONE;

    if (null_or_empty(n)) {
	; /* OK */
    } else {
	/* screened already: must be string */
	const char *s = n->v.str;

	if (!strcmp(s, "kendall")) {
	    opt = OPT_K;
	} else if (!strcmp(s, "spearman")) {
	    opt = OPT_S;
	} else {
	    p->err = E_INVARG;
	}
    }

    return opt;
}

static NODE *npcorr_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const double *x = NULL, *y = NULL;
	gretlopt opt = OPT_NONE;
	int n1 = 0, n2 = 0;

	if (l->t == SERIES) {
	    x = l->v.xvec + p->dset->t1;
	    n1 = sample_size(p->dset);
	} else {
	    n1 = gretl_vector_get_length(l->v.m);
	    if (n1 == 0) {
		p->err = E_INVARG;
	    } else {
		x = l->v.m->val;
	    }
	}

	if (!p->err && m->t == SERIES) {
	    y = m->v.xvec + p->dset->t1;
	    n2 = sample_size(p->dset);
	} else if (!p->err) {
	    n2 = gretl_vector_get_length(m->v.m);
	    if (n2 == 0) {
		p->err = E_INVARG;
	    } else {
		y = m->v.m->val;
	    }
	}

	if (!p->err && n1 != n2) {
	    p->err = E_NONCONF;
	} else if (!p->err) {
	    opt = get_npcorr_option(r, p);
	}

	if (!p->err) {
	    if (opt & OPT_S) {
		ret->v.m = spearman_rho_func(x, y, n1, &p->err);
	    } else {
		ret->v.m = kendall_tau_func(x, y, n1, &p->err);
	    }
	}
    }

    return ret;
}

/* takes two series or two matrices as arguments */

static NODE *mxtab_func (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	if (l->t == MAT && r->t == MAT) {
	    ret->v.m = matrix_matrix_xtab(l->v.m, r->v.m, &p->err);
	} else if (l->t == SERIES && r->t == SERIES) {
	    const double *x = l->v.xvec;
	    const double *y = r->v.xvec;

	    ret->v.m = gretl_matrix_xtab(p->dset->t1, p->dset->t2,
					 x, y, &p->err);
	} else {
	    p->err = E_TYPES;
	}
    }

    return ret;
}

static int type_translate_to_int (GretlType type)
{
    if (type == GRETL_TYPE_DOUBLE) {
	return 1;
    } else if (type == GRETL_TYPE_SERIES) {
	return 2;
    } else if (type == GRETL_TYPE_MATRIX) {
	return 3;
    } else if (type == GRETL_TYPE_STRING) {
	return 4;
    } else if (type == GRETL_TYPE_BUNDLE) {
	return 5;
    } else if (type == GRETL_TYPE_ARRAY) {
	return 6;
    } else if (type == GRETL_TYPE_LIST) {
	return 7;
    } else {
	return 0;
    }
}

static NODE *object_status (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	const char *s = n->v.str;

	ret->v.xval = NADBL;

	if (f == F_ISNULL) {
	    ret->v.xval = 1;
	    if (gretl_is_series(s, p->dset)) {
		ret->v.xval = 0.0;
	    } else if (gretl_is_user_var(s)) {
		ret->v.xval = 0.0;
	    }
	} else if (f == F_TYPEOF || f == F_EXISTS) {
	    GretlType type = user_var_get_type_by_name(s);

	    if (type == 0 && gretl_is_series(s, p->dset)) {
		type = GRETL_TYPE_SERIES;
	    }
	    ret->v.xval = type_translate_to_int(type);
	} else if (f == F_ISSTRING) {
	    /* temporarily reinstated */
	    ret->v.xval = gretl_is_string(s);
	} else if (f == F_ISDISCR) {
	    int v = current_series_index(p->dset, s);

	    if (v >= 0) {
		ret->v.xval = series_is_discrete(p->dset, v);
	    }
	} else if (f == F_OBSNUM) {
	    int t = get_observation_number(s, p->dset);

	    if (t > 0) {
		ret->v.xval = t;
	    }
	} else if (f == F_STRLEN) {
	    /* ret->v.xval = strlen(s); */
	    ret->v.xval = g_utf8_strlen(s, -1);
	} else if (f == F_NLINES) {
	    ret->v.xval = count_lines(s);
	} else if (f == F_REMOVE) {
	    gretl_maybe_switch_dir(s);
	    ret->v.xval = gretl_remove(s);
	}
    }

    return ret;
}

static NODE *generic_typeof_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);
    GretlType t = GRETL_TYPE_NONE;

    switch (n->t) {
    case NUM:
	t = GRETL_TYPE_DOUBLE;
	break;
    case SERIES:
	t = GRETL_TYPE_SERIES;
	break;
    case MAT:
	t = GRETL_TYPE_MATRIX;
	break;
    case STR:
	t = GRETL_TYPE_STRING;
	break;
    case BUNDLE:
	t = GRETL_TYPE_BUNDLE;
	break;
    case ARRAY:
	t = GRETL_TYPE_ARRAY;
	break;
    case LIST:
	t = GRETL_TYPE_LIST;
	break;
    default:
	break;
    }
	
    ret->v.xval = type_translate_to_int(t);

    return ret;
}

/* return scalar node holding the number of elements in
   the object associated with node @n
*/

static NODE *n_elements_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	if (n->t == MAT) {
	    gretl_matrix *m = n->v.m;

	    ret->v.xval = m->rows * m->cols;
	} else if (n->t == ARRAY) {
	    gretl_array *a = n->v.a;

	    ret->v.xval = gretl_array_get_length(a);
	} else if (n->t == BUNDLE) {
	    gretl_bundle *b = n->v.b;

	    ret->v.xval = gretl_bundle_get_n_members(b);
	} else if (ok_list_node(n)) {
	    int *list = node_get_list(n, p);

	    if (list != NULL) {
		ret->v.xval = list[0];
		free(list);
	    }
	} else if (n->t == STR) {
	    /* backward compatibility: _name_ of list */
	    int *list = get_list_by_name(n->v.str);

	    if (list != NULL) {
		ret->v.xval = list[0];
	    } else {
		p->err = E_TYPES;
	    }
	} else {
	    p->err = E_TYPES;
	}
    }

    return ret;
}

static int look_up_vname (const char *s, const DATASET *dset)
{
    int i;

    for (i=0; i<dset->v; i++) {
	if (!strcmp(s, dset->varname[i])) {
	    return i;
	}
    }

    return -1;
}

/* return scalar node holding the position of the series
   associated with node @r in the list associated with node
   @l, or zero if the series is not present in the list
*/

static NODE *in_list_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (p->err == 0 && (p->dset == NULL || p->dset->v == 0)) {
	p->err = E_NODATA;
    }

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(l, p);

	if (list != NULL) {
	    int k = -1;

	    if (useries_node(r)) {
		k = r->vnum;
	    } else if (r->t == NUM) {
		if (r->v.xval >= 0 && r->v.xval < p->dset->v) {
		    k = (int) r->v.xval;
		}
	    } else if (r->t == STR) {
		k = look_up_vname(r->v.str, p->dset);
		if (k < 0) {
		    ret->v.xval = 0;
		}
	    } else {
		node_type_error(F_INLIST, 2, SERIES, r, p);
	    }
	    if (k >= 0) {
		ret->v.xval = in_gretl_list(list, k);
	    }
	    free(list);
	}
    }

    return ret;
}

static NODE *argname_from_uvar (NODE *n, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	const char *vname = NULL;

	if (n->t == SERIES) {
	    vname = p->dset->varname[n->vnum];
	} else if (n->t == NUM) {
	    vname = n->vname;
	}

	if (vname == NULL) {
	    p->err = E_DATA;
	} else {
	    ret->v.str = gretl_func_get_arg_name(vname, &p->err);
	}
    }

    return ret;
}

static NODE *varnum_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	if (n->t == STR) {
	    int v = current_series_index(p->dset, n->v.str);

	    ret->v.xval = (v >= 0)? v : NADBL;
	} else {
	    p->err = E_DATA;
	}
    }

    return ret;
}

static NODE *int_to_string_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	int i;

	if (scalar_node(n)) {
	    i = node_get_int(n, p);
	} else if (n->t == SERIES && f == F_VARNAME) {
	    i = n->vnum;
	} else {
	    node_type_error(f, 0, NUM, n, p);
	    return NULL;
	}

	if (f == F_OBSLABEL) {
	    ret->v.str = retrieve_date_string(i, p->dset, &p->err);
	} else if (f == F_VARNAME) {
	    if (i >= 0 && i < p->dset->v) {
		ret->v.str = gretl_strdup(p->dset->varname[i]);
	    } else {
		p->err = E_DATA;
	    }
	} else {
	    p->err = E_DATA;
	}

	if (!p->err && ret->v.str == NULL) {
	    p->err = E_ALLOC;
	}
    }

    return ret;
}

static NODE *list_to_string_func (NODE *n, int f, parser *p)
{
    NODE *ret;

    if (f == F_VARNAMES) {
	ret = aux_array_node(p);
    } else {
	ret = aux_string_node(p);
    }

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(n, p);

	if (p->err) {
	    return ret;
	}

	if (f == F_VARNAME) {
	    ret->v.str = gretl_list_get_names(list, p->dset,
					      &p->err);
	} else if (f == F_VARNAMES) {
	    char **S = gretl_list_get_names_array(list, p->dset,
						  &p->err);
	    int ns = list[0];

	    if (!p->err) {
		ret->v.a = gretl_array_from_strings(S, ns, 0, &p->err);
	    }
	} else {
	    p->err = E_DATA;
	}

	free(list);
    }

    return ret;
}

/* handles both getenv (string value of variable) and
   ngetenv (numerical value of variable)
*/

static NODE *do_getenv (NODE *l, int f, parser *p)
{
    NODE *ret = (f == F_GETENV)? aux_string_node(p) :
	aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	int defined = 0;
	char *estr;

	estr = gretl_getenv(l->v.str, &defined, &p->err);

	if (f == F_GETENV) {
	    ret->v.str = estr;
	} else {
	    /* ngetenv */
	    if (defined) {
		char *test = NULL;
		double x;

		errno = 0;
		x = strtod(estr, &test);
		if (*test == '\0' && errno == 0) {
		    ret->v.xval = x;
		}
	    }
	    free(estr);
	}
    }

    return ret;
}

static NODE *single_string_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	const char *s = n->v.str;

	if (f == F_ARGNAME) {
	    ret->v.str = gretl_func_get_arg_name(s, &p->err);
	} else if (f == F_BACKTICK) {
	    ret->v.str = gretl_backtick(s, &p->err);
	} else if (f == F_STRSTRIP) {
	    ret->v.str = gretl_strstrip_copy(s, &p->err);
	} else if (f == F_FIXNAME) {
	    ret->v.str = calloc(VNAMELEN, 1);
	    normalize_join_colname(ret->v.str, s, 0);
	} else {
	    p->err = E_DATA;
	}
    }

    return ret;
}

static NODE *readfile_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	const char *codeset = NULL;

	if (!null_or_empty(r)) {
	    if (r->t == STR) {
		codeset = r->v.str;
	    } else {
		node_type_error(F_READFILE, 2, STR, r, p);
	    }
	}

	if (!p->err) {
	    ret->v.str = retrieve_file_content(l->v.str, codeset, &p->err);
	}
    }

    return ret;
}

static void strstr_escape (char *s)
{
    int i, n = strlen(s);

    for (i=0; i<n; i++) {
	if (s[i] == '\\' && (i == 0 || s[i-1] != '\\')) {
	    if (s[i+1] == 'n') {
		s[i] = '\n';
		shift_string_left(s + i + 1, 1);
		i++;
	    } else if (s[i+1] == 't') {
		s[i] = '\t';
		shift_string_left(s + i + 1, 1);
		i++;
	    }
	}
    }
}

static NODE *two_string_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret;

    if (starting(p)) {
	const char *sl = l->v.str;
	const char *sr = r->v.str;

	ret = aux_string_node(p);

	if (ret == NULL) {
	    return NULL;
	}

	if (f == F_STRSTR) {
	    char *sret, *tmp = gretl_strdup(sr);

	    if (tmp != NULL) {
		strstr_escape(tmp);
		sret = strstr(sl, tmp);
		if (sret != NULL) {
		    ret->v.str = gretl_strdup(sret);
		} else {
		    ret->v.str = gretl_strdup("");
		}
		free(tmp);
	    }
	} else if (f == B_HCAT) {
	    int n1 = strlen(l->v.str);
	    int n2 = strlen(r->v.str);

	    ret->v.str = malloc(n1 + n2 + 1);
	    if (ret->v.str != NULL) {
		*ret->v.str = '\0';
		strcat(ret->v.str, l->v.str);
		strcat(ret->v.str, r->v.str);
	    }
	} else if (f == F_JSONGET) {
	    char *(*jfunc) (const char *, const char *,
			    int *, int *);

	    jfunc = get_plugin_function("json_get");
	    if (jfunc == NULL) {
		p->err = E_FOPEN;
	    } else {
		ret->v.str = jfunc(l->v.str, r->v.str, NULL, &p->err);
	    }
	} else if (f == F_XMLGET) {
	    char *(*xfunc) (const char *, void *, GretlType,
			    int *, int *);

	    xfunc = get_plugin_function("xml_get");
	    if (xfunc == NULL) {
		p->err = E_FOPEN;
	    } else {
		if (r->t == ARRAY) {
		    ret->v.str = xfunc(l->v.str, r->v.a, GRETL_TYPE_ARRAY,
				       NULL, &p->err);
		} else {
		    ret->v.str = xfunc(l->v.str, r->v.str, GRETL_TYPE_STRING,
				       NULL, &p->err);
		}
	    }
	} else {
	    p->err = E_DATA;
	}

	if (!p->err && ret->v.str == NULL) {
	    p->err = E_ALLOC;
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *one_string_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	char *s;

	if (f == F_TOLOWER) {
	    s = ret->v.str = gretl_strdup(n->v.str);
	    while (s && *s) {
		*s = tolower(*s);
		s++;
	    }
	} else if (f == F_TOUPPER) {
	    s = ret->v.str = gretl_strdup(n->v.str);
	    while (s && *s) {
		*s = toupper(*s);
		s++;
	    }
	} else {
	    p->err = E_DATA;
	}

	if (!p->err && ret->v.str == NULL) {
	    p->err = E_ALLOC;
	}
    }

    return ret;
}

static const char *advance_to_split (const char *s)
{
    int n = strcspn(s, " \t\r\n");
    int m = strlen(s);

    return (n == m)? NULL : s + n;
}

static NODE *strsplit_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const char *s = l->v.str;
	int k = 0;

	if (r->t != EMPTY) {
	    k = node_get_int(r, p);
	    if (k < 1) {
		p->err = E_DATA;
	    } else {
		ret = aux_string_node(p);
	    }
	} else {
	    ret = aux_array_node(p);
	}

	if (!p->err && k > 0) {
	    /* returning a single string */
	    const char *q;
	    int i;

	    for (i=1; i<k; i++) {
		q = advance_to_split(s);
		if (q != NULL) {
		    q += strspn(q, " \t\r\n");
		    s = q;
		} else {
		    s = "";
		    break;
		}
	    }
	    ret->v.str = gretl_strndup(s, strcspn(s, " \t\r\n"));
	    if (!p->err && ret->v.str == NULL) {
		p->err = E_ALLOC;
	    }
	} else if (!p->err) {
	    /* returning an array of strings */
	    char **S;
	    int ns = 0;

	    S = gretl_string_split_quoted(s, &ns, NULL, &p->err);
	    if (!p->err) {
		ret->v.a = gretl_array_from_strings(S, ns, 0, &p->err);
	    }
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *errmsg_node (NODE *l, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	const char *src = NULL;

	if (l->t == NUM) {
	    int errval = l->v.xval;

	    if (errval < 0 || errval >= E_MAX) {
		p->err = E_DATA;
	    } else {
		src = errmsg_get_with_default(errval);
	    }
	} else {
	    /* empty input node */
	    src = gretl_errmsg_get();
	}

	if (src != NULL) {
	    ret->v.str = gretl_strdup(src);
	    if (ret->v.str == NULL) {
		p->err = E_ALLOC;
	    }
	}
    }

    return ret;
}

static NODE *isodate_node (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (!scalar_node(l) && l->t != SERIES) {
	node_type_error(f, 1, NUM, l, p);
    } else if (!scalar_node(r) && r->t != EMPTY) {
	node_type_error(f, 2, NUM, r, p);
    }

    if (!p->err) {
	int julian = (f == F_JULDATE);
	
	if (scalar_node(l)) {
	    /* epoch day node is scalar */
	    int as_string = scalar_node(r)? node_get_int(r, p) : 0;

	    if (!p->err) {
		ret = as_string ? aux_string_node(p) : aux_scalar_node(p);
	    }
	    if (ret != NULL) {
		double x = node_get_scalar(l, p);

		if (!as_string && xna(x)) {
		    ret->v.xval = NADBL;
		} else if (x >= 1 && x <= UINT_MAX) {
		    if (as_string) {
			ret->v.str = ymd_extended_from_epoch_day((guint32) x,
								 julian, &p->err);
		    } else {
			ret->v.xval = ymd_basic_from_epoch_day((guint32) x,
							       julian, &p->err);
		    }
		} else {
		    p->err = E_INVARG;
		}
	    }
	} else {
	    /* epoch day node is series */
	    ret = aux_series_node(p);
	    if (ret != NULL) {
		double xt;
		int t;

		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    xt = l->v.xvec[t];
		    if (na(xt)) {
			ret->v.xvec[t] = NADBL;
		    } else if (xt >= 1 && xt <= UINT_MAX) {
			ret->v.xvec[t] = ymd_basic_from_epoch_day((guint32) xt,
								  julian, &p->err);
		    } else {
			p->err = E_INVARG;
			break;
		    }
		}
	    }
	}
    }

    return ret;
}

static NODE *atof_node (NODE *l, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	char *endptr = NULL;

	errno = 0;
	gretl_push_c_numeric_locale();
	ret->v.xval = strtod(l->v.str, &endptr);
	if (errno || endptr == l->v.str) {
	    errno = 0;
	    ret->v.xval = NADBL;
	}
	gretl_pop_c_numeric_locale();
    }

    return ret;
}

static void strip_newline (char *s)
{
    if (s != NULL && *s != '\0') {
	int i, len = strlen(s);

	for (i=len-1; i>=0; i--) {
	    if (s[i] == '\n' || s[i] == '\r') {
		s[i] = '\0';
	    } else {
		break;
	    }
	}
    }
}

static NODE *getline_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	const char *buf = l->v.str;

	if (null_or_empty(r)) {
	    bufgets_finalize(buf);
	} else if (l->vname == NULL) {
	    gretl_errmsg_set("getline: the source must be a named string variable");
	    p->err = E_INVARG;
	} else if (r->vname == NULL) {
	    gretl_errmsg_set("getline: the target must be a named string variable");
	    p->err = E_INVARG;
	} else {
	    p->err = query_bufgets_init(buf);
	    if (!p->err) {
		size_t len = bufgets_peek_line_length(buf);

		if (len == 0) {
		    bufgets_finalize(buf);
		    r->v.str = user_string_reset(r->vname, &p->err);
		    ret->v.xval = 0;
		} else {
		    r->v.str = user_string_resize(r->vname, len, &p->err);
		    if (!p->err) {
			bufgets(r->v.str, len, buf);
			strip_newline(r->v.str);
			ret->v.xval = 1;
		    }
		}
	    }
	}
    }

    return ret;
}

static int series_get_start (int t1, int t2, const double *x)
{
    int t;

    for (t=t1; t<=t2; t++) {
	if (!xna(x[t])) {
	    break;
	}
    }

    return t + 1;
}

static int series_get_end (int t1, int t2, const double *x)
{
    int t;

    for (t=t2; t>=t1; t--) {
	if (!xna(x[t])) {
	    break;
	}
    }

    return t + 1;
}

static void cast_to_series (NODE *n, int f, gretl_matrix **tmp,
			    int *t1, int *t2, parser *p)
{
    gretl_matrix *m = n->v.m;
    int len = gretl_vector_get_length(m);

    if (gretl_is_null_matrix(m)) {
	p->err = E_DATA;
    } else if (len > 0 && len == p->dset->n) {
	*tmp = m;
	n->v.xvec = m->val;
    } else if (len > 0 && t1 != NULL && t2 != NULL) {
	*tmp = m;
	n->v.xvec = m->val;
	*t1 = 0;
	*t2 = len - 1;
    } else {
	node_type_error(f, 1, SERIES, n, p);
    }
}

/* Functions taking a series or vector as argument and returning
   a scalar.
*/

static NODE *
series_scalar_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *tmp = NULL;
	int t1 = p->dset->t1;
	int t2 = p->dset->t2;
	const double *x;

	if (n->t == MAT) {
	    if (f == F_SUM) {
		/* we'll sum all elements of a matrix */
		tmp = n->v.m;
		t1 = 0;
		t2 = tmp->rows * tmp->cols - 1;
		n->v.xvec = tmp->val;
	    } else if (f == F_T1 || f == F_T2) {
		cast_to_series(n, f, &tmp, NULL, NULL, p);
	    } else {
		cast_to_series(n, f, &tmp, &t1, &t2, p);
	    }
	    if (p->err) {
		return NULL;
	    }
	}

	x = n->v.xvec;

	switch (f) {
	case F_SUM:
	    ret->v.xval = gretl_sum(t1, t2, x);
	    break;
	case F_SUMALL:
	    ret->v.xval = series_sum_all(t1, t2, x);
	    break;
	case F_MEAN:
	    ret->v.xval = gretl_mean(t1, t2, x);
	    break;
	case F_SD:
	    ret->v.xval = gretl_stddev(t1, t2, x);
	    break;
	case F_VCE:
	    ret->v.xval = gretl_variance(t1, t2, x);
	    break;
	case F_SST:
	    ret->v.xval = gretl_sst(t1, t2, x);
	    break;
	case F_SKEWNESS:
	    ret->v.xval = gretl_skewness(t1, t2, x);
	    break;
	case F_KURTOSIS:
	    ret->v.xval = gretl_kurtosis(t1, t2, x);
	    break;
	case F_MIN:
	    ret->v.xval = gretl_min(t1, t2, x);
	    break;
	case F_MAX:
	    ret->v.xval = gretl_max(t1, t2, x);
	    break;
	case F_MEDIAN:
	    ret->v.xval = gretl_median(t1, t2, x);
	    break;
	case F_GINI:
	    ret->v.xval = gretl_gini(t1, t2, x);
	    break;
	case F_NOBS:
	    ret->v.xval = series_get_nobs(t1, t2, x);
	    break;
	case F_ISCONST:
	    ret->v.xval = gretl_isconst(t1, t2, x);
	    break;
	case F_ISDUMMY:
	    ret->v.xval = gretl_isdummy(t1, t2, x);
	    break;
	case F_T1:
	    ret->v.xval = series_get_start(0, p->dset->n - 1, x);
	    break;
	case F_T2:
	    ret->v.xval = series_get_end(0, p->dset->n - 1, x);
	    break;
	default:
	    break;
	}

	if (n->t == MAT) {
	    n->v.m = tmp;
	}
    }

    return ret;
}

static gretlopt get_normtest_option (NODE *n, parser *p)
{
    gretlopt opt = OPT_NONE;

    if (null_or_empty(n)) {
	; /* OK */
    } else if (n->t == STR) {
	const char *s = n->v.str;

	if (!strcmp(s, "swilk")) {
	    opt = OPT_W;
	} else if (!strcmp(s, "lillie")) {
	    opt = OPT_L;
	} else if (!strcmp(s, "jbera")) {
	    opt = OPT_J;
	} else if (!strcmp(s, "dhansen")) {
	    opt = OPT_D;
	} else {
	    p->err = E_INVARG;
	}
    } else {
	p->err = E_TYPES;
    }

    return opt;
}

static NODE *series_matrix_node (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const double *x = NULL;
	int n = 0, t1 = 0, t2 = 0;
	gretlopt opt = OPT_NONE;

	if (f == F_NORMTEST) {
	    opt = get_normtest_option(r, p);
	    if (p->err) {
		return NULL;
	    }
	}

	if (l->t == MAT) {
	    n = gretl_vector_get_length(l->v.m);
	    if (n == 0) {
		p->err = E_TYPES;
	    } else {
		x = l->v.m->val;
		t2 = n - 1;
	    }
	} else if (f == F_NORMTEST) {
	    x = l->v.xvec;
	    t1 = p->dset->t1;
	    t2 = p->dset->t2;
	} else {
	    x = l->v.xvec + p->dset->t1;
	    n = sample_size(p->dset);
	}

	if (!p->err) {
	    ret = aux_matrix_node(p);
	}

	if (!p->err) {
	    if (f == F_NORMTEST) {
		ret->v.m = gretl_normtest_matrix(x, t1, t2, opt, &p->err);
	    } else {
		ret->v.m = empirical_cdf(x, n, &p->err);
	    }
	}
    }

    return ret;
}

/* There must be a matrix in @l; @r may hold a vector or
   a scalar value */

static NODE *matrix_quantiles_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *pmat = NULL;
	int free_pmat = 0;

	if (r->t == MAT) {
	    pmat = r->v.m;
	} else {
	    double x = node_get_scalar(r, p);

	    pmat = gretl_matrix_from_scalar(x);
	    free_pmat = 1;
	}

	ret = aux_matrix_node(p);
	if (ret != NULL) {
	    ret->v.m = gretl_matrix_quantiles(l->v.m, pmat, &p->err);
	}

	if (free_pmat) {
	    gretl_matrix_free(pmat);
	}
    }

    return ret;
}

/* functions taking a series and a scalar as arguments and returning
   a scalar
*/

static NODE *
series_scalar_scalar_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	double rval = node_get_scalar(r, p);
	const double *xvec;
	int t1 = p->dset->t1;
	int t2 = p->dset->t2;
	int pd = 1;

	if (l->t == MAT) {
	    int n = gretl_vector_get_length(l->v.m);

	    if (n == 0) {
		p->err = E_TYPES;
		return NULL;
	    }
	    t1 = 0;
	    t2 = n - 1;
	    xvec = l->v.m->val;
	} else {
	    /* got a series on the left */
	    pd = p->dset->pd;
	    xvec = l->v.xvec;
	}

	ret = aux_scalar_node(p);
	if (p->err) {
	    return ret;
	}

	switch (f) {
	case F_LRVAR:
	    ret->v.xval = gretl_long_run_variance(t1, t2, xvec, (int) rval);
	    break;
	case F_QUANTILE:
	    ret->v.xval = gretl_quantile(t1, t2, xvec, rval, OPT_NONE, &p->err);
	    break;
	case F_NPV:
	    ret->v.xval = gretl_npv(t1, t2, xvec, rval, pd, &p->err);
	    break;
	case F_ISCONST:
	    ret->v.xval = panel_isconst(t1, t2, pd, xvec, (int) rval);
	    break;
	default:
	    break;
	}

    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *isconst_or_dum_node (NODE *l, NODE *r, parser *p, int f)
{
    if (f == F_ISDUMMY || null_or_empty(r)) {
	return series_scalar_func(l, f, p);
    } else if (l->t == MAT) {
	node_type_error(f, 1, SERIES, l, p);
	return NULL;
    } else if (!dataset_is_panel(p->dset)) {
	p->err = E_PDWRONG;
	return NULL;
    } else {
	return series_scalar_scalar_func(l, r, f, p);
    }
}

/* Series on left, scalar or string on right, as in
   x[23] or somevar["CA"]. We return the selected
   scalar value from the series -- unless the series
   is string-valued, in which case we return the
   string value for the given observation.
*/

static NODE *series_obs (NODE *l, NODE *r, parser *p)
{
    int strval = stringvec_node(l);
    NODE *ret;

    ret = strval ? aux_string_node(p) : aux_scalar_node(p);

    if (ret != NULL) {
	int t = -1; /* invalid */

	if (r->t == STR) {
	    t = dateton(r->v.str, p->dset);
	} else {
	    /* plain integer */
	    t = node_get_int(r, p);
	    if (!p->err && t > 0 && t <= p->dset->n) {
		t--; /* convert to zero based */
	    }
	}

	if (!p->err) {
	    if (t >= 0 && t < p->dset->n) {
		if (strval) {
		    const char *s =
			series_get_string_for_obs(p->dset, l->vnum, t);

		    ret->v.str = gretl_strdup(s);
		} else {
		    ret->v.xval = l->v.xvec[t];
		}
		ret->vnum = l->vnum; /* added 2013-09-14 */
	    } else if (strval) {
		ret->v.str = gretl_strdup("");
	    } else {
		ret->v.xval = NADBL;
	    }
	}
    }

    return ret;
}

static NODE *series_ljung_box (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	const double *x = l->v.xvec;
	int k = node_get_int(r, p);
	int t1 = p->dset->t1;
	int t2 = p->dset->t2;

	if (!p->err && k <= 0) {
	    gretl_errmsg_sprintf(_("Invalid lag order %d"), k);
	    p->err = E_DATA;
	}

	if (!p->err) {
	    p->err = series_adjust_sample(x, &t1, &t2);
	}

	if (!p->err) {
	    ret->v.xval = ljung_box(k, t1, t2, x, &p->err);
	}
    }

    return ret;
}

static NODE *series_polyfit (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
	const double *x = l->v.xvec;
	int order = node_get_int(r, p);

	if (!p->err) {
	    p->err = poly_trend(x, ret->v.xvec, p->dset, order);
	}
    }

    return ret;
}

static NODE *series_lag (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;
    const double *x = l->v.xvec;
    int k = -(node_get_int(r, p));
    int t1, t2;

    if (!p->err) {
	ret = aux_series_node(p);
    }

    if (ret == NULL) {
	return NULL;
    }

    t1 = (autoreg(p))? p->obs : p->dset->t1;
    t2 = (autoreg(p))? p->obs : p->dset->t2;

    lag_calc(ret->v.xvec, x, k, t1, t2, B_ASN, 1.0, p);

    return ret;
}

static NODE *series_sort_by (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
	if (l->t == SERIES && r->t == SERIES) {
	    p->err = gretl_sort_by(l->v.xvec, r->v.xvec, ret->v.xvec, p->dset);
	} else {
	    p->err = E_TYPES;
	}
    }

    return ret;
}

static NODE *vector_sort (NODE *l, int f, parser *p)
{
    NODE *ret = (l->t == SERIES)? aux_series_node(p) :
	aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	if (l->t == SERIES) {
	    p->err = sort_series(l->v.xvec, ret->v.xvec, f, p->dset);
	} else if (gretl_is_null_matrix(l->v.m)) {
	    p->err = E_DATA;
	} else {
	    int descending = (f == F_DSORT);

	    ret->v.m = gretl_vector_sort(l->v.m, descending, &p->err);
	}
    }

    return ret;
}

static NODE *vector_values (NODE *l, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const double *x = NULL;
	int n = 0;

	if (l->t == SERIES) {
	    n = sample_size(p->dset);
	    x = l->v.xvec + p->dset->t1;
	} else if (!gretl_is_null_matrix(l->v.m)) {
	    n = gretl_vector_get_length(l->v.m);
	    x = l->v.m->val;
	}

	if (n > 0 && x != NULL) {
	    gretlopt opt = (f == F_VALUES)? OPT_S : OPT_NONE;

	    ret->v.m = gretl_matrix_values(x, n, opt, &p->err);
	} else {
	    p->err = E_DATA;
	}
    }

    return ret;
}

static NODE *do_irr (NODE *l, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	const double *x = NULL;
	int pd = 1, n = 0;

	if (l->t == SERIES) {
	    n = sample_size(p->dset);
	    x = l->v.xvec + p->dset->t1;
	    pd = p->dset->pd;
	} else if (!gretl_is_null_matrix(l->v.m)) {
	    n = gretl_vector_get_length(l->v.m);
	    x = l->v.m->val;
	}

	if (n > 0 && x != NULL) {
	    ret->v.xval = gretl_irr(x, n, pd, &p->err);
	} else {
	    p->err = E_DATA;
	}
    }

    return ret;
}

/* Takes a series as argument and returns a matrix:
   right now only F_FREQ does this
*/

static NODE *series_matrix_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL) {
	gretl_matrix *tmp = NULL;
	int t1 = p->dset->t1;
	int t2 = p->dset->t2;

	if (n->t == MAT) {
	    cast_to_series(n, f, &tmp, &t1, &t2, p);
	}

	if (!p->err) {
	    ret->v.m = freqdist_matrix(n->v.xvec, t1, t2, &p->err);
	    if (n->t == MAT) {
		/* restore matrix on @n after "cast" above */
		n->v.m = tmp;
	    }
	}
    }

    return ret;
}

static int get_logtrans (const char *s)
{
    if (s != NULL) {
	if (*s != 'T' && strchr(s, 'l')) {
	    return 1;
	}
    }

    return 0;
}

#define use_tramo(s) (s != NULL && (s[0] == 't' || s[0] == 'T'))

#define is_panel_stat(f) (f == F_PNOBS || \
			  f == F_PMIN ||  \
			  f == F_PMAX ||  \
			  f == F_PSUM || \
			  f == F_PMEAN || \
			  f == F_PXSUM ||  \
			  f == F_PXNOBS ||  \
			  f == F_PSD)

/* Functions taking a series as argument and returning a series.
   Note that the 'r' node may contain an auxiliary parameter;
   in that case the aux value should be a scalar, unless
   we're doing F_DESEAS, in which case it should be a string,
   or one of the panel stats functions, in which case it should
   be a series.
*/

static NODE *series_series_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;
    int rtype = NUM; /* the optional right-node type */

    if (f == F_SDIFF && !dataset_is_seasonal(p->dset)) {
	p->err = E_PDWRONG;
	return NULL;
    }

    if (f == F_DESEAS) {
	rtype = STR;
    } else if (is_panel_stat(f)) {
	rtype = SERIES;
    }

    if (null_or_empty(r)) {
	rtype = 0; /* not present, OK */
    } else if (rtype == NUM) {
	if (!scalar_node(r)) {
	    node_type_error(f, 2, rtype, r, p);
	    return NULL;
	}
    } else if (r->t != rtype) {
	node_type_error(f, 2, rtype, r, p);
	return NULL;
    }

    ret = aux_series_node(p);

    if (ret != NULL) {
	gretl_matrix *tmp = NULL;
	double parm = NADBL;
	const double *z = NULL;
	const double *x;
	double *y;

	if (l->t == MAT) {
	    cast_to_series(l, f, &tmp, NULL, NULL, p);
	}

	if (rtype == SERIES) {
	    z = r->v.xvec;
	} else if (rtype == NUM) {
	    parm = node_get_scalar(r, p);
	}

	if (p->err) {
	    return NULL;
	}

	x = l->v.xvec;
	y = ret->v.xvec;

	switch (f) {
	case F_HPFILT:
	    p->err = hp_filter(x, y, p->dset, parm, OPT_NONE);
	    break;
	case F_FRACDIFF:
	    p->err = fracdiff_series(x, y, parm, 1, (autoreg(p))? p->obs : -1, p->dset);
	    break;
	case F_FRACLAG:
	    p->err = fracdiff_series(x, y, parm, 0, (autoreg(p))? p->obs : -1, p->dset);
	    break;
	case F_BOXCOX:
	    p->err = boxcox_series(x, y, parm, p->dset);
	    break;
	case F_DIFF:
	case F_LDIFF:
	case F_SDIFF:
	    p->err = diff_series(x, y, f, p->dset);
	    break;
	case F_ODEV:
	    p->err = orthdev_series(x, y, p->dset);
	    break;
	case F_CUM:
	    p->err = cum_series(x, y, p->dset);
	    break;
	case F_DESEAS:
	    if (rtype == STR) {
		int tramo = use_tramo(r->v.str);
		int logt = get_logtrans(r->v.str);

		p->err = seasonally_adjust_series(x, y, p->dset, tramo, logt);
	    } else {
		p->err = seasonally_adjust_series(x, y, p->dset, 0, 0);
	    }
	    break;
	case F_TRAMOLIN:
	    p->err = tramo_linearize_series(x, y, p->dset);
	    break;
	case F_RESAMPLE:
	    if (rtype == NUM) {
		p->err = block_resample_series(x, y, parm, p->dset);
	    } else {
		p->err = resample_series(x, y, p->dset);
	    }
	    break;
	case F_PNOBS:
	case F_PMIN:
	case F_PMAX:
	case F_PSUM:
	case F_PMEAN:
	case F_PXSUM:
	case F_PXNOBS:
	case F_PSD:
	    p->err = panel_statistic(x, y, p->dset, f, z);
	    break;
	case F_RANKING:
	    p->err = rank_series(x, y, F_SORT, p->dset);
	    break;
	default:
	    break;
	}

	if (l->t == MAT) {
	    l->v.m = tmp;
	}
    }

    return ret;
}

static NODE *do_panel_shrink (NODE *l, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	ret->v.m = panel_shrink(l->v.xvec, p->dset, &p->err);
    }

    return ret;
}

static NODE *do_panel_expand (NODE *l, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
	p->err = panel_expand(l->v.m, ret->v.xvec, OPT_NONE, p->dset);
    }

    return ret;
}

/* pergm function takes series or column vector arg, returns matrix:
   if we come up with more functions on that pattern, the following
   could be extended
*/

static NODE *pergm_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (!empty_or_num(r)) {
	/* optional 'r' node must be scalar */
	node_type_error(F_PERGM, 2, NUM, r, p);
    } else if (l->t == MAT && gretl_vector_get_length(l->v.m) == 0) {
	/* if 'l' node is not a series, must be a vector */
	node_type_error(F_PERGM, 1, SERIES, l, p);
    } else {
	ret = aux_matrix_node(p);
    }

    if (!p->err) {
	const double *x = NULL;
	int t1 = 0, t2 = 0;
	int width = -1;

	if (l->t == SERIES) {
	    x = l->v.xvec;
	    t1 = p->dset->t1;
	    t2 = p->dset->t2;
	} else if (l->t == MAT) {
	    x = l->v.m->val;
	    t1 = 0;
	    t2 = gretl_vector_get_length(l->v.m) - 1;
	}

	if (r != NULL && r->t == NUM) {
	    width = r->v.xval;
	}

	ret->v.m = periodogram_matrix(x, t1, t2, width, &p->err);
    }

    return ret;
}

/* application of scalar function to each element of matrix */

static NODE *apply_matrix_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const gretl_matrix *m = n->v.m;
	int i, n = m->rows * m->cols;
	double x;

	if (node_allocate_matrix(ret, m->rows, m->cols, p)) {
	    return NULL;
	}

	for (i=0; i<n && !p->err; i++) {
	    /* FIXME error handling? */
	    x = real_apply_func(m->val[i], f, p);
	    ret->v.m->val[i] = x;
	}
    }

    return ret;
}

static gretl_matrix *matrix_from_scalars (GPtrArray *a, int m,
					  int nsep, int seppos,
					  parser *p)
{
    gretl_matrix *M;
    NODE *n;
    int r = nsep + 1;
    int c = (seppos > 0)? seppos : m;
    int nelem = m - nsep;
    double x;
    int i, j, k;

    /* check that all rows are the same length */

    if (nelem != r * c) {
	p->err = E_PARSE;
    } else if (nsep > 0) {
	k = 0;
	for (i=0; i<m; i++) {
	    n = g_ptr_array_index(a, i);
	    if (n->t == EMPTY) {
		if (i - k != seppos) {
		    p->err = E_PARSE;
		    break;
		}
		k = i + 1;
	    }
	}
    }

    if (p->err) {
	pprintf(p->prn, _("Matrix specification is not coherent"));
	pputc(p->prn, '\n');
	return NULL;
    }

#if EDEBUG
    fprintf(stderr, "matrix_from_scalars: m=%d, nsep=%d, seppos=%d, nelem=%d\n",
	    m, nsep, seppos, nelem);
#endif

    M = gretl_matrix_alloc(r, c);
    if (M == NULL) {
	p->err = E_ALLOC;
    } else {
	k = 0;
	for (i=0; i<r && !p->err; i++) {
	    for (j=0; j<c; j++) {
		n = g_ptr_array_index(a, k++);
		if (n->t == EMPTY) {
		    n = g_ptr_array_index(a, k++);
		}
		x = node_get_scalar(n, p);
		if (na(x)) {
		    x = M_NA;
		}
		gretl_matrix_set(M, i, j, x);
	    }
	}
    }

    return M;
}

static int *full_series_list (const DATASET *dset, int *err)
{
    int *list = NULL;

    if (dset->v < 2) {
	*err = E_DATA;
	return NULL;
    }

    list = gretl_consecutive_list_new(1, dset->v - 1);
    if (list == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    return list;
}

static gretl_matrix *real_matrix_from_list (const int *list,
					    const DATASET *dset,
					    parser *p)
{
    gretl_matrix *M;

    if (list != NULL && list[0] == 0) {
	M = gretl_null_matrix_new();
    } else {
	const gretl_matrix *mmask = get_matrix_mask();

	if (mmask != NULL) {
	    M = gretl_matrix_data_subset_special(list, dset,
						 mmask, &p->err);
	} else {
	    int missop;

	    if (libset_get_bool(SKIP_MISSING)) {
		missop = M_MISSING_SKIP;
	    } else {
		missop = M_MISSING_OK;
	    }

	    M = gretl_matrix_data_subset(list, dset, dset->t1, dset->t2,
					 missop, &p->err);
	}
    }

    return M;
}

static gretl_matrix *matrix_from_list (NODE *n, parser *p)
{
    gretl_matrix *M = NULL;
    int *list = NULL;
    int freelist = 0;

    if (n != NULL) {
	if (n->t == LIST) {
	    list = n->v.ivec;
	} else {
	    p->err = E_DATA;
	}
    } else {
	list = full_series_list(p->dset, &p->err);
	freelist = 1;
    }

    if (!p->err) {
	M = real_matrix_from_list(list, p->dset, p);
    }

    if (freelist) {
	free(list);
    }

    return M;
}

static void *arg_get_data (NODE *n, int ref, GretlType *type)
{
    void *data = NULL;

    if (n->t == SERIES) {
	if (ref) {
	    *type = GRETL_TYPE_SERIES_REF;
	    data = &n->vnum;
	} else if (n->vname != NULL) {
	    /* FIXME conditionality here? */
	    *type = GRETL_TYPE_USERIES;
	    data = &n->vnum;
	} else {
	    *type = GRETL_TYPE_SERIES;
	    data = n->v.xvec;
	}
    } else if (n->t == NUM) {
	*type = ref ? GRETL_TYPE_SCALAR_REF : GRETL_TYPE_DOUBLE;
	data = &n->v.xval;
    } else if (n->t == MAT) {
	*type = ref ? GRETL_TYPE_MATRIX_REF : GRETL_TYPE_MATRIX;
	data = n->v.m;
    } else if (n->t == BUNDLE) {
	*type = ref ? GRETL_TYPE_BUNDLE_REF : GRETL_TYPE_BUNDLE;
	data = n->v.b;
    } else if (n->t == ARRAY) {
	*type = ref ? GRETL_TYPE_ARRAY_REF : GRETL_TYPE_ARRAY;
	data = n->v.a;
    } else if (n->t == STR) {
	*type = ref ? GRETL_TYPE_STRING_REF : GRETL_TYPE_STRING;
	data = n->v.str;
    } else if (n->t == LIST) {
	*type = GRETL_TYPE_LIST;
	data = n->v.ivec;
    } else {
	*type = GRETL_TYPE_NONE;
    }

    return data;
}

static NODE *suitable_ufunc_ret_node (parser *p,
				      GretlType t)
{
    if (t == GRETL_TYPE_DOUBLE) {
	return aux_scalar_node(p);
    } else if (t == GRETL_TYPE_SERIES) {
	return aux_empty_series_node(p);
    } else if (t == GRETL_TYPE_MATRIX) {
	return aux_matrix_node(p);
    } else if (t == GRETL_TYPE_LIST) {
	return aux_list_node(p);
    } else if (t == GRETL_TYPE_STRING) {
	return aux_string_node(p);
    } else if (t == GRETL_TYPE_BUNDLE) {
	return aux_bundle_node(p);
    } else if (gretl_array_type(t)) {
	return aux_array_node(p);
    } else {
	p->err = E_TYPES;
	return NULL;
    }
}

#define ok_ufunc_sym(s) (s == NUM || s == SERIES || s == MAT || \
                         s == LIST || s == U_ADDR || s == DUM || \
                         s == STR || s == EMPTY || s == BUNDLE || \
			 s == ARRAY)

/* evaluate a user-defined function */

static NODE *eval_ufunc (NODE *t, parser *p)
{
    NODE *l = t->v.b2.l;
    NODE *r = t->v.b2.r;
    NODE *save_aux = p->aux;
    NODE *ret = NULL;
    const char *funname = l->vname;
    ufunc *uf = l->v.ptr;
    fncall *fc = NULL;
    GretlType rtype = 0;
    int exec_done = 0;
    int i, nparam, argc = 0;

    rtype = user_func_get_return_type(uf);

    if (!p->err && rtype == GRETL_TYPE_VOID) {
	if (p->targ == UNK && p->lh.name[0] == '\0') {
	    /* can we ever get here? */
	    p->targ = EMPTY;
	} else if (p->targ != EMPTY) {
	    gretl_errmsg_sprintf(_("The function %s does not return any value"),
				 funname);
	    p->err = E_TYPES;
	}
    }

    if (!p->err) {
	/* get the argument and param counts */
	argc = r->v.bn.n_nodes;
	nparam = fn_n_params(uf);
	if (argc > nparam) {
	    gretl_errmsg_sprintf(_("Number of arguments (%d) does not "
				   "match the number of\nparameters for "
				   "function %s (%d)"),
				 argc, funname, nparam);
	    p->err = E_DATA;
	}
    }

    if (p->err) {
	/* no sense in continuing */
	return NULL;
    }

#if 1 /* for now, just warn */
    if (t == p->tree && (p->flags & P_CATCH)) {
	gretl_warnmsg_set(_("\"catch\" should not be used on calls to "
			    "user-defined functions"));
    }
#else
    if (t == p->tree && (p->flags & P_CATCH)) {
	p->err = E_BADCATCH;
	return NULL;
    }
#endif

    fc = fncall_new(uf);
    if (fc == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

    /* evaluate the function argument nodes */

    for (i=0; i<argc && !p->err; i++) {
	NODE *arg, *ni = r->v.bn.n[i];
	GretlType argt = 0;
	int reftype = 0;
	void *data;

	if (starting(p)) {
	    /* evaluate all nodes */
	    arg = eval(ni, p);
	} else if (ni->vname != NULL) {
	    /* otherwise let named variables through "as is" */
	    arg = ni;
	} else {
	    arg = eval(ni, p);
	}

	if (p->err || arg == NULL) {
	    fprintf(stderr, "%s: failed to evaluate arg %d\n", funname, i);
	    p->err = (p->err == 0)? E_DATA : p->err;
	} else if (!ok_ufunc_sym(arg->t)) {
	    gretl_errmsg_sprintf("%s: invalid argument type %s", funname,
				 typestr(arg->t));
	    p->err = E_TYPES;
	}

#if EDEBUG
	fprintf(stderr, "%s: arg %d is of type %s (err=%d)\n", funname, i,
		arg == NULL ? "?" : getsymb(arg->t), p->err);
#endif

	if (!p->err && arg->t == U_ADDR) {
	    /* address node: switch to the 'content' sub-node */
	    reftype = 1;
	    arg = arg->v.b1.b;
	}

	if (!p->err && arg->t == DUM && arg->v.idnum != DUM_NULL) {
	    p->err = E_TYPES;
	}

	if (!p->err) {
	    /* assemble info and push argument */
	    data = arg_get_data(arg, reftype, &argt);
	    p->err = push_function_arg(fc, arg->vname, argt, data);
	}

	if (p->err) {
	    fprintf(stderr, "%s: error evaluating arg %d\n", funname, i);
	}
    }

    /* try sending args to function */

    if (!p->err) {
	char **pdescrip = NULL;
	double xret = NADBL;
	double *Xret = NULL;
	gretl_matrix *mret = NULL;
	gretl_bundle *bret = NULL;
	gretl_array *aret = NULL;
	char *sret = NULL;
	int *iret = NULL;
	void *retp = NULL;

	if (rtype == GRETL_TYPE_DOUBLE) {
	    retp = &xret;
	} else if (rtype == GRETL_TYPE_SERIES) {
	    retp = &Xret;
	} else if (rtype == GRETL_TYPE_MATRIX) {
	    if (p->targ == UNK && p->tree == t) {
		/* target type not specified, and function returns
		   a matrix -> set target type to matrix
		*/
		p->targ = MAT;
	    }
	    retp = &mret;
	} else if (rtype == GRETL_TYPE_LIST) {
	    if (p->targ == EMPTY && p->tree == t) {
		/* this function offers a list return, but the
		   caller hasn't assigned it and it's not being
		   used as an argument to a further function, so
		   ignore the return value
		*/
		;
	    } else {
		retp = &iret;
	    }
	} else if (rtype == GRETL_TYPE_STRING) {
	    retp = &sret;
	} else if (rtype == GRETL_TYPE_BUNDLE) {
	    retp = &bret;
	} else if (gretl_array_type(rtype)) {
	    retp = &aret;
	}

	if ((p->flags & P_UFRET) && rtype == GRETL_TYPE_SERIES) {
	    /* pick up description of generated series, if any */
	    pdescrip = &p->lh.label;
	}

	p->err = gretl_function_exec(fc, rtype, p->dset, retp,
				     pdescrip, p->prn);
	exec_done = 1;

	if (!p->err && retp != NULL) {
	    reset_p_aux(p, save_aux);
	    ret = suitable_ufunc_ret_node(p, rtype);
	}

	if (!p->err && ret != NULL) {
	    if (rtype == GRETL_TYPE_DOUBLE) {
		ret->v.xval = xret;
	    } else if (rtype == GRETL_TYPE_SERIES) {
		if (ret->v.xvec != NULL) {
		    free(ret->v.xvec);
		}
		ret->v.xvec = Xret;
	    } else if (rtype == GRETL_TYPE_MATRIX) {
		if (is_tmp_node(ret)) {
		    gretl_matrix_free(ret->v.m);
		}
		ret->v.m = mret;
	    } else if (rtype == GRETL_TYPE_LIST) {
		if (is_tmp_node(ret)) {
		    free(ret->v.ivec);
		}
		if (iret != NULL) {
		    ret->v.ivec = iret;
		} else {
		    ret->v.ivec = gretl_list_new(0);
		}
	    } else if (rtype == GRETL_TYPE_STRING) {
		if (is_tmp_node(ret)) {
		    free(ret->v.str);
		}
		ret->v.str = sret;
	    } else if (rtype == GRETL_TYPE_BUNDLE) {
		if (is_tmp_node(ret)) {
		    gretl_bundle_destroy(ret->v.b);
		}
		ret->t = BUNDLE;
		ret->v.b = bret;
	    } else if (gretl_array_type(rtype)) {
		if (is_tmp_node(ret)) {
		    gretl_array_destroy(ret->v.a);
		}
		ret->t = ARRAY;
		ret->v.a = aret;
	    }
	}
    }

    if (!exec_done) {
	/* avoid leaking memory */
	fncall_destroy(fc);
    }

#if EDEBUG
    fprintf(stderr, "eval_ufunc: p->err = %d, ret = %p\n",
	    p->err, (void *) ret);
#endif

    return ret;
}

#ifdef USE_RLIB

/* evaluate an R function */

static NODE *eval_Rfunc (NODE *t, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *l = t->v.b2.l;
    NODE *r = t->v.b2.r;
    int i, argc = r->v.bn.n_nodes;
    const char *funname = l->v.str;
    int rtype = GRETL_TYPE_NONE;
    NODE *ret = NULL;

    /* first find the function */
    p->err = gretl_R_get_call(funname, argc);
    if (p->err) {
	return NULL;
    }

    /* evaluate the function arguments */
    for (i=0; i<argc && !p->err; i++) {
	NODE *arg = eval(r->v.bn.n[i], p);

	if (arg == NULL) {
	    fprintf(stderr, "%s: failed to evaluate arg %d\n", funname, i);
	} else if (!ok_ufunc_sym(arg->t)) {
	    fprintf(stderr, "%s: node type %d: not OK\n", funname, arg->t);
	    p->err = E_TYPES;
	}

	if (p->err) {
	    break;
	}

#if EDEBUG
	fprintf(stderr, "%s: arg[%d] is of type %d\n", funname, i, arg->t);
#endif

	if (arg->t == NUM) {
	    p->err = gretl_R_function_add_scalar(arg->v.xval);
	} else if (arg->t == SERIES) {
	    gretl_matrix *m = tmp_matrix_from_series(arg, p);

	    if (m != NULL) {
		p->err = gretl_R_function_add_matrix(m);
		gretl_matrix_free(m);
	    }
	} else if (arg->t == MAT) {
	    p->err = gretl_R_function_add_matrix(arg->v.m);
	} else {
	    fprintf(stderr, "eval_Rfunc: argument not supported\n");
	    p->err = E_TYPES;
	    return NULL;
	}
	if (p->err) {
	    fprintf(stderr, "eval_Rfunc: error evaluating arg %d\n", i);
	}
    }

    /* try sending args to function */

    if (!p->err) {
	double xret = NADBL;
	void *retp = &xret;

	p->err = gretl_R_function_exec(funname, &rtype, &retp);

	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    if (gretl_scalar_type(rtype)) {
		ret = aux_scalar_node(p);
		if (ret != NULL) {
		    ret->v.xval = xret;
		}
	    } else if (rtype == GRETL_TYPE_MATRIX) {
		ret = aux_matrix_node(p);
		if (ret != NULL) {
		    if (is_tmp_node(ret)) {
			gretl_matrix_free(ret->v.m);
		    }
		    ret->v.m = (gretl_matrix *) retp;
		}
	    }
	}
    }

#if EDEBUG
    fprintf(stderr, "eval_Rfunc: p->err = %d, ret = %p\n",
	    p->err, (void *) ret);
#endif

    return ret;
}

#endif

static NODE *model_var_via_accessor (const char *key, parser *p)
{
    NODE *ret = NULL;
    gchar *dkey;
    int mv;

    dkey = g_strdup_printf("$%s", key);

    if ((mv = mvar_lookup(dkey)) > 0) {
	/* OK: @key (with dollar-sign prepended) is a
	   regular model-related accessor
	*/
	NODE tmp;

	memset(&tmp, 0, sizeof tmp);
	tmp.t = MVAR;
	tmp.v.idnum = mv;
	ret = object_var_node(&tmp, p);
    }

    g_free(dkey);

    return ret;
}

/* Getting an object from within a bundle: on the left is the
   bundle reference, on the right should be a string -- the
   key to look up to get content.
*/

static NODE *get_bundle_member (NODE *l, NODE *r, parser *p)
{
    char *key = r->v.str;
    GretlType type;
    int copied = 0;
    int size = 0;
    void *val = NULL;
    NODE *ret = NULL;

#if EDEBUG
    if (l->t == DBUNDLE) {
	fprintf(stderr, "get_bundle_member: %s[\"%s\"]\n",
		bvarname(l->v.idnum), key);
    } else {
	fprintf(stderr, "get_bundle_member: %s[\"%s\"]\n",
		l->vname, key);
    }
#endif

    if (l->t == DBUNDLE) {
	if (l->v.idnum == B_MODEL) {
	    /* Treating the 'last model' as a bundle: we'll first
	       try for a regular model-related accessor then,
	       unless we already hit an error, for a "special"
	       (a named item of model data that's not represented
	       by a regular accessor).
	    */
	    ret = model_var_via_accessor(key, p);
	    if (ret != NULL || p->err) {
		return ret;
	    } else {
		val = last_model_get_data(key, &type, &size, &copied, &p->err);
	    }
	} else if (l->v.idnum == B_SYSINFO) {
	    val = sysinfo_bundle_get_data(key, &type, &p->err);
	}
    } else {
	val = gretl_bundle_get_data(l->v.b, key, &type, &size, &p->err);
    }

    if (p->err) {
	return ret;
    }

    if (type == GRETL_TYPE_INT) {
	ret = aux_scalar_node(p);
	if (ret != NULL) {
	    int *ip = val;

	    ret->v.xval = *ip;
	}
    } else if (type == GRETL_TYPE_DOUBLE) {
	ret = aux_scalar_node(p);
	if (ret != NULL) {
	    double *dp = val;

	    ret->v.xval = *dp;
	}
    } else if (type == GRETL_TYPE_STRING) {
	ret = string_pointer_node(p);
	if (ret != NULL) {
	    ret->v.str = (char *) val;
	}
    } else if (type == GRETL_TYPE_MATRIX) {
	if (copied) {
	    ret = aux_matrix_node(p);
	} else {
	    ret = matrix_pointer_node(p);
	}
	if (ret != NULL) {
	    ret->v.m = (gretl_matrix *) val;
	}
    } else if (type == GRETL_TYPE_BUNDLE) {
	if (copied) {
	    ret = aux_bundle_node(p);
	} else {
	    ret = bundle_pointer_node(p);
	}
	if (ret != NULL) {
	    ret->v.b = (gretl_bundle *) val;
	}
    } else if (type == GRETL_TYPE_ARRAY) {
	ret = array_pointer_node(p);
	if (ret != NULL) {
	    ret->v.a = (gretl_array *) val;
	}
    } else if (type == GRETL_TYPE_SERIES) {
	const double *x = val;

	if (size <= p->dset->n) {
	    ret = aux_series_node(p);
	    if (ret != NULL) {
		int t;

		for (t=p->dset->t1; t<=p->dset->t2 && t<size; t++) {
		    ret->v.xvec[t] = x[t];
		}
	    }
	} else if (size > 0) {
	    ret = aux_matrix_node(p);
	    if (ret != NULL) {
		ret->v.m = gretl_vector_from_array(x, size,
						   GRETL_MOD_NONE);
		if (ret->v.m == NULL) {
		    p->err = E_ALLOC;
		}
	    }
	} else {
	    p->err = E_DATA;
	}
    } else if (type == GRETL_TYPE_LIST) {
	p->err = stored_list_check((const int *) val, p->dset);
	if (!p->err) {
	    if (copied) {
		ret = aux_list_node(p);
	    } else {
		ret = list_pointer_node(p);
	    }
	}
	if (!p->err) {
	    ret->v.ivec = (int *) val;
	}
    } else {
	p->err = E_DATA;
    }

    return ret;
}

static NODE *test_bundle_key (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	gretl_bundle *bundle = l->v.b;
	const char *key = r->v.str;
	GretlType type = 0;
	int err = 0;

	gretl_bundle_get_data(bundle, key, &type, NULL, &err);
	ret->v.xval = type_translate_to_int(type);
	if (err) {
	    gretl_error_clear();
	}
    }

    return ret;
}

static const char *optional_bundle_get (gretl_bundle *b,
					const char *key,
					double *px,
					int *err)
{
    const char *s = NULL;

    if (!*err) {
	/* proceed only if we haven't already hit an error */
	if (px != NULL) {
	    *px = gretl_bundle_get_scalar(b, key, err);
	} else {
	    s = gretl_bundle_get_string(b, key, err);
	}
	if (*err == E_DATA) {
	    /* non-existence of item (E_DATA) is OK, but
	       wrong type (E_TYPES) is not
	    */
	    gretl_error_clear();
	    *err = 0;
	}
    }

    return s;
}

static NODE *curl_bundle_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

#ifndef USE_CURL
    gretl_errmsg_set(_("Internet access not supported"));
    p->err = E_DATA;
#else
    if (ret != NULL) {
	gretl_bundle *b = NULL;
	int curl_err = 0;

	if (n->t != U_ADDR) {
	    p->err = E_TYPES;
	} else {
	    /* switch to 'content' sub-node */
	    n = n->v.b1.b;
	    if (n->t != BUNDLE) {
		p->err = E_TYPES;
	    } else {
		b = n->v.b;
	    }
	}

	if (!p->err) {
	    const char *url = NULL;
	    const char *header = NULL;
	    const char *postdata = NULL;
	    char *output = NULL;
	    char *errmsg = NULL;
	    double xinclude = 0;

	    url = gretl_bundle_get_string(b, "URL", &p->err);
	    header = optional_bundle_get(b, "header", NULL, &p->err);
	    postdata = optional_bundle_get(b, "postdata", NULL, &p->err);
	    optional_bundle_get(b, "include", &xinclude, &p->err);

	    if (!p->err) {
		int include = (xinclude == 1.0);

		curl_err = gretl_curl(url, header, postdata, include,
				      &output, &errmsg);
	    }
	    if (output != NULL) {
		p->err = gretl_bundle_set_string(b, "output", output);
		free(output);
	    } else if (errmsg != NULL) {
		p->err = gretl_bundle_set_string(b, "errmsg", errmsg);
		free(errmsg);
	    }
	}

	if (!p->err) {
	    ret->v.xval = curl_err;
	}
    }
#endif /* curl supported in libgretl */

    return ret;
}

static NODE *svm_predict_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_series_node(p);

#ifndef HAVE_LIBSVM
    gretl_errmsg_set(_("Libsvm is not supported"));
    p->err = E_DATA;
#else
    if (ret != NULL) {
	gretl_bundle *rbundle = NULL;

	if (!null_or_empty(r)) {
	    /* we should have a bundle pointer on the right */
	    if (r->t == U_ADDR) {
		r = r->v.b1.b;
		if (r->t == BUNDLE) {
		    rbundle = r->v.b;
		} else {
		    p->err = E_INVARG;
		}
	    } else {
		p->err = E_INVARG;
	    }
	}

	if (!p->err) {
	    int (*pfunc) (const int *, gretl_bundle *,
			  gretl_bundle *, double *,
			  int *, DATASET *, PRN *);
	    int got_yhat = 0;

	    pfunc = get_plugin_function("gretl_svm_predict");
	    if (pfunc == NULL) {
		p->err = E_FOPEN;
	    } else {
		p->err = pfunc(l->v.ivec, m->v.b, rbundle, ret->v.xvec,
			       &got_yhat, p->dset, p->prn);
		if (!p->err && !got_yhat) {
		    /* change the return type to scalar NA */
		    free(ret->v.xvec);
		    ret->t = NUM;
		    ret->v.xval = NADBL;
		}
	    }
	}
    }
#endif    

    return ret;
}

static gretl_bundle *bvar_get_bundle (NODE *n, parser *p)
{
    gretl_bundle *b = NULL;

    if (n->v.idnum == B_MODEL) {
	b = bundle_from_model(NULL, p->dset, &p->err);
    } else if (n->v.idnum == B_SYSINFO) {
	gretl_bundle *tmp = get_sysinfo_bundle(&p->err);

	if (!p->err) {
	    b = gretl_bundle_copy(tmp, &p->err);
	}
    } else {
	p->err = E_DATA;
    }

    return b;
}

static NODE *dollar_bundle_node (NODE *n, parser *p)
{
    NODE *ret = aux_bundle_node(p);

    if (ret != NULL) {
	ret->v.b = bvar_get_bundle(n, p);
    }

    return ret;
}

static NODE *type_string_node (NODE *n, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL) {
	int t = node_get_int(n, p);

	if (!p->err) {
	    const char *s = "null";

	    if (t == 1) {
		s = "scalar";
	    } else if (t == 2) {
		s = "series";
	    } else if (t == 3) {
		s = "matrix";
	    } else if (t == 4) {
		s = "string";
	    } else if (t == 5) {
		s = "bundle";
	    } else if (t == 6) {
		s = "array";
	    } else if (t == 7) {
		s = "list";
	    }

	    ret->v.str = gretl_strdup(s);
	    if (ret->v.str == NULL) {
		p->err = E_ALLOC;
	    }
	}
    }

    return ret;
}

static double *scalar_to_series (NODE *n, parser *p)
{
    double *ret = NULL;
    int t;

    if (p->dset == NULL || p->dset->n == 0) {
	p->err = E_NODATA;
    } else {
	ret = malloc(p->dset->n * sizeof *ret);
	if (ret == NULL) {
	    p->err = E_ALLOC;
	} else {
	    for (t=0; t<p->dset->n; t++) {
		if (t >= p->dset->t1 && t <= p->dset->t2) {
		    ret[t] = n->v.xval;
		} else {
		    ret[t] = NADBL;
		}
	    }
	}
    }

    return ret;
}

static GretlType gretl_type_of (int t)
{
    if (t == NUM) {
	return GRETL_TYPE_DOUBLE;
    } else if (t == SERIES) {
	return GRETL_TYPE_SERIES;
    } else if (t == MAT) {
	return GRETL_TYPE_MATRIX;
    } else if (t == STR) {
	return GRETL_TYPE_STRING;
    } else if (t == BUNDLE) {
	return GRETL_TYPE_BUNDLE;
    } else if (t == LIST) {
	return GRETL_TYPE_LIST;
    } else if (t == ARRAY) {
	return GRETL_TYPE_ARRAY;
    } else {
	return 0;
    }
}

/* note: we come here only if we setting a bundle-member or
   an element of an array
*/

static int lhs_type_check (GretlType spec, GretlType rhs,
			   int t)
{
    int err = 0;

    if (spec != 0 && spec != rhs) {
	if (t == BUNDLE) {
	    gretl_errmsg_sprintf(_("Expected %s but got %s"),
				 gretl_type_get_name(spec),
				 gretl_type_get_name(rhs));
	} else {
	    gretl_errmsg_sprintf(_("Specified type %s does not match array type %s"),
				 gretl_type_get_name(spec),
				 gretl_type_get_name(rhs));
	}
	err = E_TYPES;
    }

    return err;
}

static void *get_mod_assign_result (void *lp, GretlType ltype,
				    NODE *r, int *size, parser *p)
{
    void *ret = NULL;
    NODE *l, *op;

    if (p->op == INC || p->op == DEC) {
	/* handle increment/decrement postfix operator */
	if (ltype == GRETL_TYPE_DOUBLE) {
	    double x = *(double *) lp;

	    if (!na(x)) {
		x += (p->op == INC)? 1 : -1;
		*(double *) lp = x;
	    }
	    ret = lp;
	} else {
	    p->err = E_TYPES;
	}
	return ret; /* handled */
    }

    /* create binary tree to hold left and right */
    l = newempty();
    op = newb2(p->op, l, r);

    /* if that went OK, put the relevant type specifier
       and pointer onto the LHS sub-node
    */
    if (op == NULL || l == NULL) {
	p->err = E_ALLOC;
    } else if (ltype == GRETL_TYPE_MATRIX) {
	l->t = MAT;
	l->v.m = lp;
    } else if (ltype == GRETL_TYPE_DOUBLE) {
	l->t = NUM;
	l->v.xval = *(double *) lp;
    } else if (ltype == GRETL_TYPE_STRING) {
	l->t = STR;
	l->v.str = lp;
	if (p->op == B_ADD) {
	    /* reinterpret '+' when appending to a string */
	    op->t = B_HCAT;
	}
    } else if (ltype == GRETL_TYPE_BUNDLE) {
	l->t = BUNDLE;
	l->v.b = lp;
    } else if (ltype == GRETL_TYPE_ARRAY) {
	l->t = ARRAY;
	l->v.a = lp;
    } else if (ltype == GRETL_TYPE_LIST) {
	l->t = LIST;
	l->v.ivec = lp;
	if (p->op == B_ADD) {
	    /* reinterpret '+' when appending to a list */
	    op->t = B_LCAT;
	}
    } else if (ltype == GRETL_TYPE_SERIES) {
	l->t = SERIES;
	l->v.xvec = lp;
    } else {
	p->err = E_TYPES;
    }

    if (!p->err) {
	/* FIXME parser state variables? */
	int saveflags = p->flags;
	int savetarg = p->targ;
	NODE *ev;

#if LHDEBUG
	fputs("*** op tree, before ***\n", stderr);
	print_tree(op, p, 0);
#endif
	p->targ = l->t;
	p->flags = P_START;
	ev = eval(op, p);
#if LHDEBUG
	fputs("*** ev tree, after ***\n", stderr);
	print_tree(ev, p, 0);
#endif

	if (!p->err) {
	    /* get @ret off node @ev and clean up */
	    if (ev->t == MAT) {
		ret = ev->v.m;
	    } else if (ev->t == BUNDLE) {
		ret = ev->v.b;
	    } else if (ev->t == STR) {
		ret = ev->v.str;
	    } else if (ev->t == ARRAY) {
		ret = ev->v.a;
	    } else if (ev->t == LIST) {
		ret = ev->v.ivec;
	    } else if (ev->t == SERIES) {
		/* FIXME sample range? */
		if (size != NULL) {
		    *size = p->dset->n;
		}
		ret = ev->v.xvec;
	    } else if (ev->t == NUM) {
		ret = lp;
		*(double *) lp = ev->v.xval;
	    } else {
		p->err = E_TYPES;
	    }
	}

	p->targ = savetarg;
	p->flags = saveflags;
	free(ev);
    }

    /* thought: if @p is reusable, should we try preserving the
       nodes allocated here?
    */

    /* trash temporary nodes */
    free(op);
    free(l);

    if (ret == NULL && !p->err) {
	p->err = E_DATA;
    }

    return ret;
}

/* ".=" : we need a scalar on the RHS */

static int dot_assign_to_matrix (gretl_matrix *m, parser *p,
				 int *prechecked)
{
    int err = 0;
    
    if (p->ret->t == NUM) {
	double x = p->ret->v.xval;

	if (na(x)) {
	    x = M_NA;
	    set_gretl_warning(W_GENNAN);
	}
	gretl_matrix_fill(m, x);
	if (prechecked != NULL) {
	    *prechecked = 1;
	}
    } else {
	err = E_TYPES;
    }

    return err;
}

/* Setting an object in a bundle under a given key string. We get here
   only if p->lh.expr is non-NULL.
*/

static int set_bundle_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->v.b2.l;
    NODE *lh2 = lhs->v.b2.r;
    GretlType targ = 0;
    GretlType type = 0;
    gretl_bundle *bundle;
    void *ptr = NULL;
    char *key = NULL;
    int size = 0;
    int donate = 0;
    int err = 0;

    if (lh1->t != BUNDLE) {
	return E_DATA;
    }

    bundle = lh1->v.b;
    key = lh2->v.str;

    if (bundle == NULL || key == NULL) {
	return E_DATA;
    }

#if LHDEBUG
    fprintf(stderr, "set_bundle_value: bundle = %p, key = '%s'\n",
	    (void *) bundle, key);
#endif

    if (p->op != B_ASN) {
	/* We must have an existing bundle member under @key, and
	   its type will determine the type of the result of
	   inflected assignment.
	*/
	GretlType ltype = 0;
	void *lp;

	lp = gretl_bundle_get_data(bundle, key, &ltype, &size, &err);
	if (!err) {
	    targ = gretl_type_of(p->targ);
	    err = lhs_type_check(targ, ltype, BUNDLE);
	}
	if (p->op == B_DOTASN) {
	    /* accepted only for matrices */
	    if (!err) {
		if (ltype == GRETL_TYPE_MATRIX) {
		    err = dot_assign_to_matrix(lp, p, NULL);
		} else {
		    err = E_TYPES;
		}
	    }
	    return err; /* handled */
	}
	if (!err) {
	    ptr = get_mod_assign_result(lp, ltype, rhs, &size, p);
	    err = p->err;
	    if (p->op == INC || p->op == DEC) {
		return err; /* handled */
	    }
	}
	if (!err) {
	    type = ltype;
	    if (ptr != lp) {
		donate = 1; /* donate: is this always right? */
	    }
	}
	goto push_data;
    }    

    /* Note: @targ is the gretl type specified by the caller for
       the bundle member (if any, this need not be supplied), and
       @type is the gretl type of the object arising on the RHS.
       It's an error if @targ is non-zero and @type does not
       agree with it -- except for the case where @targ is given
       as "series" and we get a suitable matrix on the right. As
       of 2015-10-03, when we get a request to put a series into
       a bundle we actually put in a matrix, which in fact makes it
       easier to get a series back out again.
    */

    if (p->targ == ARRAY) {
	targ = p->lh.gtype;
    } else {
	targ = gretl_type_of(p->targ);
    }

    if (!err) {
	switch (rhs->t) {
	case NUM:
	    if (targ == GRETL_TYPE_SERIES) {
		ptr = scalar_to_series(rhs, p);
		if (p->err) {
		    err = p->err;
		} else {
		    type = GRETL_TYPE_SERIES;
		    size = p->dset->n;
		    donate = 1;
		}
	    } else if (targ == GRETL_TYPE_MATRIX) {
		ptr = gretl_matrix_from_scalar(rhs->v.xval);
		type = GRETL_TYPE_MATRIX;
		donate = 1;
	    } else {
		ptr = &rhs->v.xval;
		type = GRETL_TYPE_DOUBLE;
	    }
	    break;
	case STR:
	    ptr = rhs->v.str;
	    type = GRETL_TYPE_STRING;
	    donate = !reusable(p) && is_tmp_node(rhs);
	    break;
	case MAT:
	    if (targ == GRETL_TYPE_DOUBLE && scalar_matrix_node(rhs)) {
		ptr = &rhs->v.m->val[0];
		type = GRETL_TYPE_DOUBLE;
	    } else if (targ == GRETL_TYPE_SERIES) {
		ptr = (double *) get_colvec_as_series(rhs, 0, p);
		if (!p->err) {
		    type = GRETL_TYPE_SERIES;
		    size = p->dset->n;
		}
	    } else if (targ == GRETL_TYPE_LIST) {
		ptr = gretl_list_from_vector(rhs->v.m, p->dset, &p->err);
		if (!p->err) {
		    type = GRETL_TYPE_LIST;
		    donate = 1;
		}
	    } else {
		ptr = rhs->v.m;
		type = GRETL_TYPE_MATRIX;
		donate = is_tmp_node(rhs);
	    }
	    break;
	case SERIES:
	    ptr = rhs->v.xvec;
	    type = GRETL_TYPE_SERIES;
	    size = p->dset->n;
	    donate = !reusable(p) && is_tmp_node(rhs);
	    break;
	case BUNDLE:
	    ptr = rhs->v.b;
	    type = GRETL_TYPE_BUNDLE;
	    break;
	case ARRAY:
	    ptr = rhs->v.a;
	    /* get more specific type for comparison with
	       what the user specified (if anything)
	    */
	    type = gretl_array_get_type(rhs->v.a);
	    donate = !reusable(p) && is_tmp_node(rhs);
	    break;
	case LIST:
	    if (targ == GRETL_TYPE_MATRIX) {
		ptr = gretl_list_to_vector(rhs->v.ivec, &p->err);
		if (!p->err) {
		    type = GRETL_TYPE_MATRIX;
		    donate = 1;
		}
	    } else {
		ptr = rhs->v.ivec;
		type = GRETL_TYPE_LIST;
		donate = is_tmp_node(rhs);
	    }
	    break;
	default:
	    err = E_TYPES;
	    break;
	}
    }

    if (!err) {
	/* check for result type-incompatible with user's spec */
	err = lhs_type_check(targ, type, BUNDLE);
    }

 push_data:

    if (!err) {
	if (gretl_is_array_type(type)) {
	    /* revert to generic array type for the functions below */
	    type = GRETL_TYPE_ARRAY;
	}
	if (donate) {
	    /* it's OK to hand over the data pointer */
	    err = gretl_bundle_donate_data(bundle, key, ptr, type, size);
	    if (ptr == rhs->v.ptr) {
		rhs->v.ptr = NULL; /* avoid freeing! */
	    }
	} else {
	    /* the data must be copied into the bundle */
	    err = gretl_bundle_set_data(bundle, key, ptr, type, size);
	}
	if (!err && type == GRETL_TYPE_MATRIX) {
	    /* for use by genr_get_output_matrix() */
	    p->lh.mret = ptr;
	}
    }

    return err;
}

static int array_index_from_mspec (matrix_subspec *spec, int *err)
{
    int idx = 0;

    if (spec->type[0] == SEL_RANGE &&
	spec->type[1] == SEL_NULL &&
	spec->sel[0].range[0] == spec->sel[0].range[1]) {
	idx = spec->sel[0].range[0];
    } else {
	gretl_errmsg_set("Invalid left-hand side index value");
	*err = E_TYPES;
    }

    return idx;
}

static int set_array_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->v.b2.l;
    NODE *lh2 = lhs->v.b2.r;
    GretlType atype = 0;
    GretlType type = 0;
    GretlType targ = 0;
    gretl_array *array = NULL;
    void *ptr = NULL;
    int idx = 0;
    int donate = 0;
    int err = 0;

    if (lh1->t != ARRAY) {
	return E_TYPES;
    }

    if (lh2->t == MSPEC) {
	/* FIXME allow a range here? */
	idx = array_index_from_mspec(lh2->v.mspec, &err);
	if (err) {
	    return err;
	}
    } else {
	idx = lh2->v.xval;
    }

    array = lh1->v.a;
    if (array == NULL) {
	return E_DATA;
    } else if (idx <= 0 || idx > gretl_array_get_length(array)) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), idx);
	return E_DATA;
    }

    atype = gretl_array_get_content_type(array);
    targ = gretl_type_of(p->targ);
    err = lhs_type_check(targ, atype, ARRAY);

#if LHDEBUG
    fprintf(stderr, "set_array_value: array %p, idx=%d, atype=%s, targ=%s, err=%d\n",
	    (void *) array, idx, gretl_type_get_name(atype),
	    gretl_type_get_name(targ), err);
#endif

    idx--; /* convert index to 0-based */

    if (!err && p->op != B_ASN) {
	GretlType ltype = 0;
	void *lp;

	lp = gretl_array_get_element(array, idx, &ltype, &err);	
	if (p->op == B_DOTASN) {
	    if (!err) {
		if (ltype == GRETL_TYPE_MATRIX) {
		    err = dot_assign_to_matrix(lp, p, NULL);
		} else {
		    err = E_TYPES;
		}
	    }
	    return err; /* handled */
	}
	if (!err) {
	    ptr = get_mod_assign_result(lp, ltype, rhs, NULL, p);
	    err = p->err;
	    if (p->op == INC || p->op == DEC) {
		return err; /* handled */
	    }
	}
	if (!err) {
	    type = ltype;
	    if (ptr != lp) {
		donate = 1; /* donate: always right? */
	    }
	}
	goto push_data;
    }

    if (!err && atype == GRETL_TYPE_LIST) {
	ptr = node_get_list(rhs, p);
	err = p->err;
	if (!err) {
	    type = GRETL_TYPE_LIST;
	    donate = 1;
	}
    } else if (!err) {
	switch (rhs->t) {
	case NUM:
	    if (atype == GRETL_TYPE_MATRIX) {
		ptr = gretl_matrix_from_scalar(rhs->v.xval);
		type = GRETL_TYPE_MATRIX;
		donate = 1;
	    }
	    break;
	case STR:
	    ptr = rhs->v.str;
	    type = GRETL_TYPE_STRING;
	    donate = !reusable(p) && is_tmp_node(rhs);
	    break;
	case MAT:
	    ptr = rhs->v.m;
	    type = GRETL_TYPE_MATRIX;
	    donate = is_tmp_node(rhs);
	    break;
	case BUNDLE:
	    ptr = rhs->v.b;
	    type = GRETL_TYPE_BUNDLE;
	    break;
	case LIST:
	    ptr = rhs->v.ivec;
	    type = GRETL_TYPE_LIST;
	    donate = is_tmp_node(rhs);
	    break;
	default:
	    err = E_TYPES;
	    break;
	}
    }

    if (!err && type != atype) {
	err = E_TYPES;
    }

 push_data:

    if (!err) {
	if (donate) {
	    /* it's OK to hand over the data pointer */
	    err = gretl_array_set_element(array, idx, ptr, type, 0);
	    if (ptr == rhs->v.ptr) {
		rhs->v.ptr = NULL; /* gone! */
	    }
	} else {
	    /* the data must be copied into the array */
	    err = gretl_array_set_element(array, idx, ptr, type, 1);
	}
	if (!err && type == GRETL_TYPE_MATRIX) {
	    /* for use by genr_get_output_matrix() */
	    p->lh.mret = ptr;
	}
    }

    return err;
}

/* setting member of list: only straight assignment is accepted */

static int set_list_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->v.b2.l;
    NODE *lh2 = lhs->v.b2.r;
    int *list = NULL;
    int idx = 0, v = -1;
    int err = 0;

    if (p->op != B_ASN) {
	gretl_errmsg_sprintf(_("'%s' : not implemented for this type"),
			     get_opstr(p->op));
	return E_TYPES;
    }

    if (lh2->t == MSPEC) {
	idx = array_index_from_mspec(lh2->v.mspec, &err);
	if (err) {
	    return err;
	}
    } else {
	idx = lh2->v.xval;
    }

    list = lh1->v.ivec;
    if (list == NULL) {
	return E_DATA;
    } else if (idx < 1 || idx > list[0]) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), idx);
	return E_DATA;
    }

#if LHDEBUG
    fprintf(stderr, "set_list_value: list = %p, idx = %d\n",
	    (void *) list, idx);
#endif

    if (rhs->t == NUM) {
	v = node_get_int(rhs, p);
    } else if (rhs->t == SERIES) {
	v = rhs->vnum;
    } else {
	p->err = E_TYPES;
    }

    if (!p->err && (v < 0 || v >= p->dset->v)) {
	gretl_errmsg_set(_("Invalid list element"));
	p->err = E_DATA;
    }

    if (!p->err) {
	list[idx] = v;
    }

    return p->err;
}

/* setting element of string: only straight assignment accepted */

static int set_string_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->v.b2.l;
    NODE *lh2 = lhs->v.b2.r;
    char *s1 = NULL;
    char *s2 = NULL;
    int idx = 0, err = 0;

    if (p->op != B_ASN) {
	gretl_errmsg_sprintf(_("'%s' : not implemented for this type"),
			     get_opstr(p->op));
	return E_TYPES;
    } else if (rhs->t != STR) {
	return E_TYPES;
    }

    if (lh2->t == MSPEC) {
	idx = array_index_from_mspec(lh2->v.mspec, &err);
	if (err) {
	    return err;
	}
    } else {
	idx = lh2->v.xval;
    }

    s1 = lh1->v.str;
    if (s1 == NULL) {
	return E_DATA;
    } else if (idx < 1 || idx > strlen(s1)) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), idx);
	return E_DATA;
    }

#if LHDEBUG
    fprintf(stderr, "set_string_value: s1 = %s, idx = %d\n",
	    s1, idx);
#endif

    s2 = rhs->v.str;
    if (strlen(s2) != 1) {
	return E_INVARG; /* FIXME */
    } else {
	s1[idx - 1] = s2[0];
    }

    return p->err;
}

static int set_series_obs_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->v.b2.l;
    NODE *lh2 = lhs->v.b2.r;
    double **Z = p->dset->Z;
    char *label = NULL;
    double x = NADBL;
    int op = p->op;
    int v, t;

    if (lh1->t != SERIES || lh2->t != NUM) {
	return E_TYPES;
    } else {
	v = lh1->vnum;
	if (v <= 0 || v >= p->dset->v) {
	    return E_DATA;
	}
	t = node_get_int(lh2, p);
	if (t < 1 || t > p->dset->n) {
	    return E_DATA;
	}
	/* convert to 0-based */
	t--;
    }

    if (rhs == NULL) {
	if (p->op == INC) {
	    x = 1;
	    op = B_ADD;
	} else if (p->op == DEC) {
	    x = 1;
	    op = B_SUB;
	} else {
	    return E_TYPES; /* ? */
	}
    } else if (rhs->t == STR) {
	if (is_string_valued(p->dset, v)) {
	    label = rhs->v.str;
	} else {
	    return E_TYPES;
	}
    } else {
	x = node_get_scalar(rhs, p);
	if (p->err) {
	    return p->err;
	}
    }

    if (is_string_valued(p->dset, v)) {
	if (label != NULL) {
	    if (op != B_ASN) {
		p->err = E_TYPES;
	    } else {
		p->err = series_set_string_val(p->dset, v, t, label);
	    }
	} else {
	    x = xy_calc(Z[v][t], x, op, NUM, p);
	    if (!p->err) {
		p->err = string_series_assign_value(p->dset, v, t, x);
	    }
	}
    } else {
	Z[v][t] = xy_calc(Z[v][t], x, op, NUM, p);
    }
    
    if (p->err == 0) {
	/* made a change to an element of a series */
	p->flags |= P_OBSVAL;
	set_dataset_is_changed();
    }

    return p->err;
}

#define ok_submatrix_op(o) (o == B_ASN  || o == B_DOTASN || \
			    o == B_ADD  || o == B_SUB ||    \
			    o == B_MUL  || o == B_DIV)

/* Here we're replacing a sub-matrix of the original LHS matrix, by
   either straight or inflected assignment. The value that we're
   using for replacement will be either a matrix or a scalar.
*/

static int set_matrix_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->v.b2.l;
    NODE *lh2 = lhs->v.b2.r;
    gretl_matrix *m1, *m2 = NULL;
    matrix_subspec *spec;
    double y = NADBL;
    int rhs_scalar = 0;
    int prechecked = 0;
    int free_m2 = 0;

    if (!ok_submatrix_op(p->op)) {
	gretl_errmsg_sprintf(_("The operator '%s' is not valid in this context"),
			     get_opstr(p->op));
	return E_TYPES;
    } else if (lh1->t != MAT) {
	fprintf(stderr, "set_matrix_value: got %s, not matrix!\n",
		getsymb(lh1->t));
	return E_DATA;
    }

    m1 = lh1->v.m;
    spec = lh2->v.mspec;

    if (m1 == NULL || spec == NULL) {
	return E_DATA;
    }

    /* check the validity of the subspec we got, and
       adjust it if need be in the light of the
       dimensions of @m.
    */
    p->err = check_matrix_subspec(spec, m1);
    if (p->err) {
	return p->err;
    }

#if EDEBUG > 1
    gretl_matrix_print(m1, "m1, in set_matrix_value");
    print_mspec(spec);
#endif

    if (scalar_node(rhs)) {
	/* single value on RHS */
	y = (rhs->t == NUM)? rhs->v.xval: rhs->v.m->val[0];
	rhs_scalar = 1;
    } else if (rhs->t == MAT) {
	/* not a scalar: get the RHS matrix */
	m2 = rhs->v.m;
    } else if (rhs->t == SERIES) {
	/* legacy: this has long been accepted */
	m2 = tmp_matrix_from_series(rhs, p);
	prechecked = 1;
	free_m2 = 1; /* flag temporary status */
    } else {
	p->err = E_TYPES;
    }

    if (p->err) {
	return p->err;
    }

    if (rhs_scalar && spec->type[0] == SEL_ELEMENT) {
    	/* assignment (possibly inflected) of a scalar value
	   to a single element of an existing matrix
	*/
	int i = mspec_get_row_index(spec);
	int j = mspec_get_col_index(spec);
	double x = matrix_get_element(m1, i, j, &p->err);

	if (!p->err) {
	    if (p->op == B_ASN) {
		x = y;
	    } else {
		x = xy_calc(x, y, p->op, MAT, p);
	    }
	    if (xna(x)) {
		if (na(x)) {
		    x = M_NA;
		}
		set_gretl_warning(W_GENNAN);
	    }
	    gretl_matrix_set(m1, i-1, j-1, x);
	}
	return p->err; /* note, we're done */
    }    

    if (rhs_scalar && p->op == B_ASN) {
	/* straight assignment of a scalar value to a
	   non-scalar submatrix */
	if (xna(y)) {
	    if (na(y)) {
		y = M_NA;
	    }
	    set_gretl_warning(W_GENNAN);
	}
	p->err = assign_scalar_to_submatrix(m1, y, spec);
	return p->err; /* note, we're done */
    }

    if (p->op != B_ASN) {
	/* Here we're doing '+=' or some such, in which case a new
	   submatrix must be calculated using the original
	   submatrix @a and the newly generated matrix (or
	   scalar value).
	*/
	gretl_matrix *a = matrix_get_submatrix(m1, spec, 1, &p->err);

	if (!p->err) {
	    if (rhs_scalar) {
		int i, n = a->rows * a->cols;

		for (i=0; i<n; i++) {
		    a->val[i] = xy_calc(a->val[i], y, p->op, MAT, p);
		}
		/* assign computed matrix to m2, and mark it for
		   freeing */
		m2 = a;
		free_m2 = 1;
		prechecked = 1;
	    } else {
		gretl_matrix *b = NULL;

		p->err = real_matrix_calc(a, m2, p->op, &b);
		gretl_matrix_free(a);
		/* replace RHS m2 with computed result */
		if (free_m2) {
		    /* m2 was temp result of series conversion */
		    gretl_matrix_free(m2);
		}
		free_m2 = 1;
		m2 = b;
	    }
	}
    }

    if (!p->err) {
	/* Write new submatrix @m2 into place: note that we come here
	   directly if none of the special conditions above are
	   satisfied -- for example, if the newly generated value
	   is a matrix and the task is straight assignment. Also
	   check for numerical "breakage" in the replacement
	   submatrix.
	*/
	p->err = matrix_replace_submatrix(m1, m2, spec);
#if MATRIX_NA_CHECK
	if (!p->err && !prechecked && gretl_matrix_xna_check(m1)) {
	    set_gretl_warning(W_GENNAN);
	}
#endif
    }

    if (free_m2) {
	gretl_matrix_free(m2);
    }

    return p->err;
}

static int set_bundle_note (gretl_bundle *b, const char *key,
			    const char *note, parser *p)
{
    int err = 0;

    if (b == NULL) {
	p->err = E_UNKVAR;
    } else {
	err = gretl_bundle_set_note(b, key, note);
    }

    return err;
}

static gretl_matrix *get_corrgm_matrix (NODE *l,
					NODE *m,
					NODE *r,
					parser *p)
{
    int xcf = (r->t != EMPTY);
    int *list = NULL;
    gretl_matrix *A = NULL;
    int k;

    /* ensure we've got an order */
    k = node_get_int(m, p);
    if (p->err) {
	return NULL;
    }

    /* if we're supposed to have a list, check that we
       actually have one */
    if (l->t != SERIES && l->t != MAT) {
	list = node_get_list(l, p);
	if (p->err) {
	    return NULL;
	}
    }

    /* if third node is matrix, must be col vector */
    if (r->t == MAT) {
	if (r->v.m->cols != 1) {
	    p->err = E_NONCONF;
	    return NULL;
	}
    }

    if (!xcf) {
	/* acf/pacf */
	if (l->t == SERIES) {
	    A = acf_matrix(l->v.xvec, k, p->dset, 0, &p->err);
	} else if (l->t == MAT) {
	    A = multi_acf(l->v.m, NULL, NULL, k, &p->err);
	} else {
	    /* it must be a list */
	    A = multi_acf(NULL, list, p->dset, k, &p->err);
	}
    } else {
	/* cross-correlogram */
	void *px = NULL, *py = NULL;
	int xtype = SERIES;

	if (list != NULL) {
	    px = list;
	    xtype = LIST;
	} else if (l->t == MAT) {
	    px = l->v.m;
	    xtype = MAT;
	} else {
	    px = l->v.xvec;
	}

	py = (r->t == MAT)? (void *) r->v.m : (void *) r->v.xvec;

	A = multi_xcf(px, xtype, py, r->t, p->dset, k, &p->err);
    }

    free(list);

    return A;
}

static const char *ptr_node_get_matrix_name (NODE *t, parser *p)
{
    const char *name = NULL;

    if (t->t == U_ADDR) {
	NODE *n = t->v.b1.b;

	if (umatrix_node(n)) {
	    name = n->vname;
	} else {
	    p->err = E_TYPES;
	}
    }

    return name;
}

static gretl_matrix *get_density_matrix (NODE *t, double bws,
					 int ctrl, parser *p)
{
    gretl_matrix *(*kdfunc) (const double *, int, double,
			     gretlopt, int *);
    gretlopt opt = ctrl ? OPT_O : OPT_NONE;
    gretl_matrix *m = NULL;
    const double *x = NULL;
    int n;

    if (t->t == SERIES) {
	n = sample_size(p->dset);
	x = t->v.xvec + p->dset->t1;
    } else {
	n = gretl_vector_get_length(t->v.m);
	if (n == 0) {
	    p->err = E_TYPES;
	} else {
	    x = t->v.m->val;
	}
    }

    if (!p->err) {
	kdfunc = get_plugin_function("kernel_density_matrix");
	if (kdfunc == NULL) {
	    p->err = E_FOPEN;
	} else {
	    m = (*kdfunc)(x, n, bws, opt, &p->err);
	}
    }

    return m;
}

static int aggregate_discrete_check (const int *list, const DATASET *dset)
{
    int i, vi;

    for (i=1; i<=list[0]; i++) {
	vi = list[i];
	if (!series_is_discrete(dset, vi) &&
	    !gretl_isdiscrete(dset->t1, dset->t2, dset->Z[vi])) {
	    gretl_errmsg_sprintf(_("The variable '%s' is not discrete"),
				 dset->varname[vi]);
	    return E_DATA;
	}
    }

    return 0;
}

static gretl_matrix *mshape_scalar (double x, int r, int c, int *err)
{
    gretl_matrix *m = gretl_matrix_alloc(r, c);

    if (m == NULL) {
	*err = E_ALLOC;
    } else {
	int i, n = r * c;

	for (i=0; i<n; i++) {
	    m->val[i] = x;
	}
    }

    return m;
}

/* evaluate a built-in function that has three arguments */

static NODE *eval_3args_func (NODE *l, NODE *m, NODE *r, int f, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *ret = NULL;
    gretl_matrix *A = NULL;
    int post_process = 1;

    if (f == F_MSHAPE) {
	if (l->t != MAT && l->t != NUM) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (!scalar_node(r)) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    int k1 = node_get_int(m, p);
	    int k2 = node_get_int(r, p);

	    if (!p->err) {
		if (l->t == NUM) {
		    A = mshape_scalar(l->v.xval, k1, k2, &p->err);
		} else {
		    A = gretl_matrix_shape(l->v.m, k1, k2, &p->err);
		}
	    }
	}
    } else if (f == F_TRIMR) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (!scalar_node(r)) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    int k1 = node_get_int(m, p);
	    int k2 = node_get_int(r, p);

	    if (!p->err) {
		A = gretl_matrix_trim_rows(l->v.m, k1, k2, &p->err);
	    }
	}
    } else if (f == F_SVD) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != U_ADDR && m->t != EMPTY) {
	    node_type_error(f, 2, U_ADDR, m, p);
	} else if (r->t != U_ADDR && r->t != EMPTY) {
	    node_type_error(f, 3, U_ADDR, r, p);
	} else {
	    const char *uname, *vname;

	    uname = ptr_node_get_matrix_name(m, p);
	    vname = ptr_node_get_matrix_name(r, p);
	    A = user_matrix_SVD(l->v.m, uname, vname, &p->err);
	}
    } else if (f == F_TOEPSOLV || f == F_VARSIMUL) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != MAT) {
	    node_type_error(f, 2, MAT, m, p);
	} else if (r->t != MAT) {
	    node_type_error(f, 3, MAT, r, p);
	} else {
	    if (f == F_TOEPSOLV) {
		A = gretl_toeplitz_solve(l->v.m, m->v.m, r->v.m, &p->err);
	    } else {
		A = gretl_matrix_varsimul(l->v.m, m->v.m, r->v.m, &p->err);
	    }
	}
    } else if (f == F_CORRGM) {
	if (l->t != SERIES && l->t != MAT && !ok_list_node(l)) {
	    node_type_error(f, 1, SERIES, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != EMPTY && r->t != SERIES && r->t != MAT) {
	    node_type_error(f, 3, SERIES, r, p);
	} else {
	    A = get_corrgm_matrix(l, m, r, p);
	}
    } else if (f == F_SEQ) {
	if (!scalar_node(l)) {
	    node_type_error(f, 1, NUM, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (!scalar_node(r) && r->t != EMPTY) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    double start = node_get_scalar(l, p);
	    double end = node_get_scalar(m, p);
	    double step = (r->t != EMPTY)? node_get_scalar(r, p) : 1.0;

	    if (!p->err) {
		A = gretl_matrix_seq(start, end, step, &p->err);
	    }
	}
    } else if (f == F_STRNCMP) {
	post_process = 0;
	if (l->t != STR) {
	    node_type_error(f, 1, STR, l, p);
	} else if (m->t != STR) {
	    node_type_error(f, 2, STR, m, p);
	} else if (!empty_or_num(r)) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		if (r != NULL && r->t == NUM) {
		    ret->v.xval = strncmp(l->v.str, m->v.str,
					  (int) r->v.xval);
		} else {
		    ret->v.xval = strcmp(l->v.str, m->v.str);
		}
	    }
	}
    } else if (f == F_WEEKDAY) {
	post_process = 0;
	if (l->t == NUM && m->t == NUM && r->t == NUM) {
	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		int yr = node_get_int(l, p);
		int mo = node_get_int(m, p);
		int day = node_get_int(r, p);

		if (!p->err) {
		    int julian = 0;
		    
		    if (yr < 0) {
			yr = -yr;
			julian = 1;
		    }
		    ret->v.xval = day_of_week(yr, mo, day, julian, &p->err);
		}
	    }
	} else if (l->t == SERIES && m->t == SERIES && r->t == SERIES) {
	    reset_p_aux(p, save_aux);
	    ret = aux_series_node(p);
	    if (ret != NULL) {
		p->err = fill_day_of_week_array(ret->v.xvec,
						l->v.xvec,
						m->v.xvec,
						r->v.xvec,
						p->dset);
	    }
	} else {
	    p->err = E_TYPES;
	}
    } else if (f == F_DAYSPAN) {
	post_process = 0;
	if (l->t == NUM && m->t == NUM && r->t == NUM) {
	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		int ed1 = l->v.xval;
		int ed2 = m->v.xval;
		int wkdays = r->v.xval;

		ret->v.xval = day_span(ed1, ed2, wkdays, &p->err);
	    }
	} else {
	    p->err = E_TYPES;
	}
    } else if (f == F_SMPLSPAN) {
	post_process = 0;
	if (l->t == STR && m->t == STR && r->t == NUM) {
	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		int pd = node_get_int(r, p);

		if (!p->err) {
		    ret->v.xval = sample_span(l->v.str, m->v.str,
					      pd, &p->err);
		}
	    }
	} else {
	    p->err = E_TYPES;
	}	
    } else if (f == F_KDENSITY) {
	if (l->t != SERIES && l->t != MAT) {
	    node_type_error(f, 1, SERIES, l, p);
	} else if (m->t != NUM && m->t != EMPTY) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != NUM && r->t != EMPTY) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    double bws = (m->t != EMPTY)? m->v.xval : 1.0;
	    int ctrl = (r->t != EMPTY)? (int) r->v.xval : 0;

	    A = get_density_matrix(l, bws, ctrl, p);
	}
    } else if (f == F_MONTHLEN) {
	post_process = 0;
	if (!scalar_node(l)) {
	    node_type_error(f, 1, NUM, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (!scalar_node(r)) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    int mo = node_get_int(l, p);
	    int yr = node_get_int(m, p);
	    int wk = node_get_int(r, p);

	    if (p->err) {
		; /* from node_get_int() */
	    } else if (mo < 1 || mo > 12 ||
		(wk != 5 && wk != 6 && wk != 7)) {
		p->err = E_INVARG;
	    } else {
		reset_p_aux(p, save_aux);
		ret = aux_scalar_node(p);
		if (!p->err) {
		    int julian = 0;

		    if (yr < 0) {
			yr = -yr;
			julian = 1;
		    }
		    ret->v.xval = get_days_in_month(mo, yr, wk, julian);
		}
	    }
	}
    } else if (f == F_SETNOTE) {
	post_process = 0;
	if (l->t != BUNDLE) {
	    node_type_error(f, 1, BUNDLE, l, p);
	} else if (m->t != STR) {
	    node_type_error(f, 2, STR, m, p);
	} else if (r->t != STR) {
	    node_type_error(f, 3, STR, r, p);
	} else {
	    reset_p_aux(p, save_aux);
	    ret = aux_scalar_node(p);
	    if (!p->err) {
		ret->v.xval = set_bundle_note(l->v.b, m->v.str,
					      r->v.str, p);
	    }
	}
    } else if (f == F_BWFILT) {
	gretl_matrix *tmp = NULL;

	post_process = 0;
	if (l->t != SERIES) {
	    if (l->t == MAT) {
		cast_to_series(l, f, &tmp, NULL, NULL, p);
	    } else {
		node_type_error(f, 1, SERIES, l, p);
	    }
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != NUM) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    ret = aux_series_node(p);
	    if (!p->err) {
		p->err = butterworth_filter(l->v.xvec, ret->v.xvec, p->dset,
					    m->v.xval, r->v.xval);
	    }
	}
	if (tmp != NULL) {
	    l->v.m = tmp;
	}
    } else if (f == F_CHOWLIN) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != MAT && r->t != EMPTY) {
	    node_type_error(f, 3, MAT, r, p);
	} else {
	    const gretl_matrix *X = (r->t == MAT)? r->v.m : NULL;

	    A = matrix_chowlin(l->v.m, X, m->v.xval, &p->err);
	}
    } else if (f == F_IRF) {
	if (l->t != NUM) {
	    node_type_error(f, 1, NUM, l, p);
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != NUM && r->t != EMPTY) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    double alpha = (r->t == NUM)? r->v.xval : 0.0;
	    int targ = (int) l->v.xval - 1;
	    int shock = (int) m->v.xval - 1;

	    A = last_model_get_irf_matrix(targ, shock, alpha,
					  p->dset, &p->err);
	}
    } else if (f == F_MLAG) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != MAT && m->t != NUM) {
	    /* scalar or vector */
	    node_type_error(f, 2, MAT, m, p);
	} else if (r->t != NUM && r->t != EMPTY) {
	    /* optional scalar */
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    double missval = (r->t == NUM)? r->v.xval : 0.0;
	    gretl_matrix *mm;

	    if (m->t == NUM) {
		/* promote arg2 if scalar */
		mm = make_scalar_matrix(m->v.xval);
	    } else {
		mm = m->v.m;
	    }
	    A = gretl_matrix_lag(l->v.m, mm, OPT_L, missval);
	    if (m->t == NUM && mm != NULL) {
		gretl_matrix_free(mm);
	    }
	}
    } else if (f == F_EIGSOLVE) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != MAT) {
	    node_type_error(f, 2, MAT, m, p);
	} else if (r->t != EMPTY && r->t != U_ADDR) {
	    /* optional matrix-pointer */
	    node_type_error(f, 3, U_ADDR, r, p);
	} else {
	    const char *rname;

	    rname = (r->t == U_ADDR)? ptr_node_get_matrix_name(r, p) : "null";
	    A = user_gensymm_eigenvals(l->v.m, m->v.m, rname, &p->err);
	}
    } else if (f == F_NADARWAT) {
	post_process = 0;
	if (l->t != SERIES) {
	    node_type_error(f, 1, SERIES, l, p);
	} else if (m->t != SERIES) {
	    node_type_error(f, 2, SERIES, m, p);
	} else if (r->t != NUM) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    reset_p_aux(p, save_aux);
	    ret = aux_series_node(p);
	    if (!p->err) {
		p->err = nadaraya_watson(l->v.xvec, m->v.xvec,
					 r->v.xval, p->dset,
					 ret->v.xvec);
	    }
	}
    } else if (f == F_PRINCOMP) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != EMPTY && r->t != NUM) {
	    /* optional boolean */
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    int cov = (r->t == EMPTY)? 0 : node_get_int(r, p);
	    int k = node_get_int(m, p);

	    if (!p->err) {
		A = gretl_matrix_pca(l->v.m, k,
				     cov ? OPT_C : OPT_NONE,
				     &p->err);
	    }
	}
    } else if (f == F_HALTON) {
	if (l->t != NUM) {
	    node_type_error(f, 1, NUM, l, p);
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != EMPTY && r->t != NUM) {
	    /* optional offset */
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    int offset = (r->t == EMPTY)? 10 : node_get_int(r, p);
	    int rows = node_get_int(l, p);
	    int cols = node_get_int(m, p);

	    if (!p->err) {
		A = halton_matrix(rows, cols, offset, &p->err);
	    }
	}
    } else if (f == F_IWISHART) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != EMPTY && r->t != NUM) {
	    /* optional number of replications */
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    const gretl_matrix *S = l->v.m;
	    int v = node_get_int(m, p);
	    int N = (r->t == EMPTY)? 0 : node_get_int(r, p);

	    if (!p->err) {
		if (N == 0) {
		    A = inverse_wishart_matrix(S, v, &p->err);
		} else {
		    A = inverse_wishart_sequence(S, v, N, &p->err);
		}
	    }
	}
    } else if (f == F_AGGRBY) {
	if (l->t != SERIES && l->t != LIST && !null_or_empty(l)) {
	    node_type_error(f, 1, SERIES, l, p);
	} else if (m->t != SERIES && m->t != LIST) {
	    node_type_error(f, 2, SERIES, m, p);
	} else if (r->t != STR && !null_or_empty(r)) {
	    node_type_error(f, 3, STR, r, p);
	} else {
	    const char *fncall = NULL;
	    const double *x = NULL;
	    const double *y = NULL;
	    const int *xlist = NULL;
	    const int *ylist = NULL;

	    if (r->t == STR) {
		fncall = r->v.str;
	    }
	    if (l->t == SERIES) {
		x = l->v.xvec;
	    } else if (l->t == LIST) {
		xlist = l->v.ivec;
	    }
	    if (m->t == SERIES) {
		y = m->v.xvec;
	    } else {
		ylist = m->v.ivec;
		p->err = aggregate_discrete_check(ylist, p->dset);
	    }

	    if (!p->err) {
		A = aggregate_by(x, y, xlist, ylist, fncall,
				 p->dset, &p->err);
	    }
	}
    } else if (f == F_SUBSTR) {
	if (l->t != STR) {
	    node_type_error(f, 1, STR, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (!scalar_node(r)) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    post_process = 0;
	    reset_p_aux(p, save_aux);
	    ret = aux_string_node(p);

	    if (ret != NULL) {
		int ini = node_get_int(m, p);
		int fin = node_get_int(r, p);

		if (!p->err) {
		    ret->v.str = gretl_substring(l->v.str, ini, fin, &p->err);
		}
	    }
	}
    } else if (f == F_MWEIGHTS) {
	if (!scalar_node(l)) {
	    node_type_error(f, 1, NUM, l, p);
	} else if (m->t != MAT) {
	    node_type_error(f, 2, MAT, m, p);
	} else if (!scalar_node(r) && r->t != STR) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    int length = node_get_int(l, p);
	    int method = node_get_midas_method(r, p);

	    if (!p->err) {
		A = midas_weights(length, m->v.m, method, &p->err);
	    }
	}
    } else if (f == F_MGRADIENT) {
	if (!scalar_node(l)) {
	    node_type_error(f, 1, NUM, l, p);
	} else if (m->t != MAT) {
	    node_type_error(f, 2, MAT, m, p);
	} else if (!scalar_node(r) && r->t != STR) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    int length = node_get_int(l, p);
	    int method = node_get_midas_method(r, p);

	    if (!p->err) {
		A = midas_gradient(length, m->v.m, method, &p->err);
	    }
	}
    }

    if (!p->err && post_process) {
	reset_p_aux(p, save_aux);
	ret = aux_matrix_node(p);
	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }
	    ret->v.m = A;
	}
    }

    return ret;
}

static NODE *eval_print_scan (NODE *l, NODE *m, NODE *r, int f, parser *p)
{
    NODE *ret;

    if (f == F_SPRINTF) {
	ret = aux_string_node(p);
    } else {
	ret = aux_scalar_node(p);
    }

    if (ret != NULL) {
	const char *fmt = m->v.str;
	const char *args = NULL;
	const char *lstr = NULL;
	int n = 0;

	if (l != NULL) {
	    /* sscanf only */
	    lstr = l->v.str;
	}

	if (!p->err) {
	    if (!null_or_empty(r)) {
		args = r->v.str;
	    }
	    if (f == F_SSCANF) {
		p->err = do_sscanf(lstr, fmt, args, p->dset, &n);
	    } else if (f == F_SPRINTF) {
		ret->v.str = do_sprintf_function(fmt, args, p->dset, &p->err);
	    } else {
		p->err = do_printf(fmt, args, p->dset, p->prn, &n);
	    }
	    if (f != F_SPRINTF && !p->err) {
		ret->v.xval = n;
	    }
	}
    }

    return ret;
}

static NODE *string_scrape_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL) {
	int comma = node_get_bool(r, p, 0);

	if (!p->err) {
	    ret->v.m = scrape_numerical_values(l->v.str, comma, &p->err);
	}
    }

    return ret;
}

static int x_to_period (double x, char c, int *julian, int *err)
{
    if (julian != NULL && c == 'y') {
	if (x < 0) {
	    *julian = 1;
	    x = -x;
	} else {
	    *julian = 0;
	}
    }
    
    if (na(x)) {
	return -1;
    } else if (x < 0 || fabs(x) > INT_MAX) {
	*err = E_INVARG;
	return -1;
    } else {
	int k = x;
	int ret = x;

	if (c == 'y' && k < 0) {
	    ret = -1;
	} else if (c == 'm' && (k < 1 || k > 12)) {
	    ret = -1;
	} else if (c == 'd' && (k < 1 || k > 31)) {
	    ret = -1;
	}

	if (ret <= 0) {
	    fprintf(stderr, "epochday: got %c = %d!\n", c, k);
	    *err = E_INVARG;
	}

	return ret;
    }
}

static NODE *eval_epochday (NODE *ny, NODE *nm, NODE *nd, parser *p)
{
    NODE *ret = NULL;
    NODE *nodes[3] = {ny, nm, nd};
    double *x[3] = {NULL, NULL, NULL};
    int ymd[3] = {-1, -1, -1};
    const char *code = "ymd";
    int julian = 0;
    double sval;
    int i;

    /* Policy: NA for y, m, or d will give an NA result;
       non-NA but out-of-bounds values for y, m or d produce
       an error.
    */

    for (i=0; i<3 && !p->err; i++) {
	if (scalar_node(nodes[i])) {
	    sval = node_get_scalar(nodes[i], p);
	    if (!p->err) {
		ymd[i] = x_to_period(sval, code[i], &julian, &p->err);
	    }
	} else if (nodes[i]->t == SERIES) {
	    x[i] = nodes[i]->v.xvec;
	} else {
	    node_type_error(F_EPOCHDAY, i+1, NUM, nodes[i], p);
	}
    }

    if (!p->err) {
	int y = ymd[0];
	int m = ymd[1];
	int d = ymd[2];

	/* From this point, -1 for y, m or d indicates either an NA
	   or the fact that a series and not a scalar was given for
	   that date-element.
	*/

	if (x[0] == NULL && x[1] == NULL && x[2] == NULL) {
	    /* y, m and d all given as scalars */
	    ret = aux_scalar_node(p);
	    if (!p->err) {
		if (y < 0 || m < 0 || d < 0) {
		    ret->v.xval = NADBL;
		} else {
		    if (julian) {
			ret->v.xval = epoch_day_from_julian_ymd(y, m, d);
		    } else {
			ret->v.xval = epoch_day_from_ymd(y, m, d);
		    }
		    if (ret->v.xval < 0) {
			p->err = E_INVARG;
		    }
		}
	    }
	} else {
	    ret = aux_series_node(p);
	    if (!p->err) {
		int t;

		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    if (x[0] != NULL) {
			y = x_to_period(x[0][t], 'y', &julian, &p->err);
		    }
		    m = (x[1] == NULL)? m : x_to_period(x[1][t], 'm', NULL, &p->err);
		    d = (x[2] == NULL)? d : x_to_period(x[2][t], 'd', NULL, &p->err);
		    if (p->err) {
			break;
		    } else if (y < 0 || m < 0 || d < 0) {
			/* got an NA somewhere */
			ret->v.xvec[t] = NADBL;
		    } else {
			if (julian) {
			    ret->v.xvec[t] = epoch_day_from_julian_ymd(y, m, d);
			} else {
			    ret->v.xvec[t] = epoch_day_from_ymd(y, m, d);
			}
			if (ret->v.xvec[t] < 0) {
			    p->err = E_INVARG;
			}
		    }
		}
	    }
	}
    }

    return ret;
}

/* Bessel function handler: the 'r' node can be of scalar, series or
   matrix type.  Right now, this only supports scalar order ('m'
   node).
*/

static NODE *eval_bessel_func (NODE *l, NODE *m, NODE *r, parser *p)
{
    char ftype;
    double v;
    NODE *ret = NULL;

    if (!starting(p) && r->t != SERIES) {
	return aux_any_node(p);
    }

    ftype = l->v.str[0];
    v = node_get_scalar(m, p);

    if (r->t == NUM) {
	ret = aux_scalar_node(p);
    	if (ret != NULL) {
	    double x = r->v.xval;

	    ret->v.xval = gretl_bessel(ftype, v, x, &p->err);
    	}
    } else if (r->t == MAT) {
	ret = aux_matrix_node(p);
	if (ret != NULL) {
	    const gretl_matrix *x = r->v.m;
	    int i, n = x->rows * x->cols;

	    if (node_allocate_matrix(ret, x->rows, x->cols, p)) {
		return NULL;
	    }

	    for (i=0; i<n && !p->err; i++) {
		ret->v.m->val[i] = gretl_bessel(ftype, v, x->val[i], &p->err);
	    }
	}
    } else if (r->t == SERIES) {
	ret = aux_series_node(p);
	if (ret != NULL) {
	    const double *x = r->v.xvec;
	    int t1 = (autoreg(p))? p->obs : p->dset->t1;
	    int t2 = (autoreg(p))? p->obs : p->dset->t2;
	    int t;

	    for (t=t1; t<=t2 && !p->err; t++) {
		ret->v.xvec[t] = gretl_bessel(ftype, v, x[t], &p->err);
	    }
	}
    }

    return ret;
}

/* Given an original value @x, see if it matches any of the @n0 values
   in @x0.  If so, return the substitute value from @x1, otherwise
   return the original.
*/

static double subst_val (double x, const double *x0, int n0,
			 const double *x1, int n1)
{
    int i;

    for (i=0; i<n0; i++) {
	if (x == x0[i]) {
	    return (n1 == 1)? *x1 : x1[i];
	} else if (isnan(x) && isnan(x0[i])) {
	    /* we'll count this as a match */
	    return (n1 == 1)? *x1 : x1[i];
	}
    }

    return x;
}

/* String search and replace: return a node containing a copy
   of the string on node @src in which all occurrences of
   the string on @n0 are replaced by the string on @n1.
   This is literal string replacement if @f is F_STRSUB,
   regular expression replacement if @f is F_REGSUB.
*/

static NODE *string_replace (NODE *src, NODE *n0, NODE *n1, int f,
			     parser *p)
{
    if (!starting(p)) {
	return aux_string_node(p);
    } else {
	NODE *ret = NULL;
	NODE *n[3] = {src, n0, n1};
	char const *S[3];
	int i;

	for (i=0; i<3; i++) {
	    /* all nodes must be of string type */
	    if (n[i]->t != STR) {
		node_type_error(f, i, STR, n[i], p);
		return NULL;
	    } else {
		S[i] = n[i]->v.str;
	    }
	}

	ret = aux_string_node(p);
	if (p->err) {
	    return NULL;
	}

	if (f == F_REGSUB) {
	    ret->v.str = gretl_regexp_replace(S[0], S[1], S[2], &p->err);
	} else {
	    ret->v.str = gretl_literal_replace(S[0], S[1], S[2], &p->err);
	}

	return ret;
    }
}

/* replace_value: non-interactive search-and-replace for series and
   matrices.  @src holds the series or matrix of which we want a
   modified copy; @n0 holds the value (or vector of values) to be
   replaced; and @n1 holds the replacement value(s). It would be nice
   to extend this to lists.
*/

static NODE *replace_value (NODE *src, NODE *n0, NODE *n1, parser *p)
{
    gretl_vector *vx0 = NULL;
    gretl_vector *vx1 = NULL;
    double x0 = 0, x1 = 0;
    int k0 = -1, k1 = -1;
    NODE *ret = NULL;

    if (!starting(p)) {
	return aux_any_node(p);
    }

    /* n0: the original value, to be replaced */
    if (n0->t == NUM) {
	x0 = n0->v.xval;
    } else if (n0->t == MAT) {
	vx0 = n0->v.m;
	if (gretl_is_null_matrix(vx0)) {
	    p->err = E_DATA;
	} else if ((k0 = gretl_vector_get_length(vx0)) == 0) {
	    p->err = E_NONCONF;
	}
    } else {
	node_type_error(F_REPLACE, 1, NUM, n0, p);
    }

    if (p->err) {
	return NULL;
    }

    /* n1: the replacement value */
    if (n1->t == NUM) {
	x1 = n1->v.xval;
    } else if (n1->t == MAT) {
	vx1 = n1->v.m;
	if (gretl_is_null_matrix(vx1)) {
	    p->err = E_DATA;
	} else if ((k1 = gretl_vector_get_length(vx1)) == 0) {
	    p->err = E_NONCONF;
	}
    } else {
	node_type_error(F_REPLACE, 2, NUM, n1, p);
    }

    if (!p->err) {
	if (n0->t == NUM && n1->t == MAT) {
	    /* can't replace scalar with vector */
	    p->err = E_TYPES;
	} else if (k0 > 0 && k1 > 0 && k0 != k1) {
	    /* if they're both vectors, they must be
	       the same length */
	    p->err = E_NONCONF;
	}
    }

    if (!p->err) {
	if (src->t == SERIES) {
	    ret = aux_series_node(p);
	} else if (src->t == MAT) {
	    ret = aux_matrix_node(p);
	} else {
	    node_type_error(F_REPLACE, 3, SERIES, src, p);
	}
    }

    if (!p->err) {
	double *px0 = (vx0 != NULL)? vx0->val : &x0;
	double *px1 = (vx1 != NULL)? vx1->val : &x1;
	double xt;
	int t;

	if (k0 < 0) k0 = 1;
	if (k1 < 0) k1 = 1;

	if (src->t == SERIES) {
	    for (t=p->dset->t1; t<=p->dset->t2; t++) {
		xt = src->v.xvec[t];
		ret->v.xvec[t] = subst_val(xt, px0, k0, px1, k1);
	    }
	} else if (src->t == MAT) {
	    gretl_matrix *m = src->v.m;
	    int n = gretl_matrix_rows(m) * gretl_matrix_cols(m);

	    ret->v.m = gretl_matrix_copy(m);
	    if (ret->v.m == NULL) {
		p->err = E_ALLOC;
	    } else {
		for (t=0; t<n; t++) {
		    xt = m->val[t];
		    ret->v.m->val[t] = subst_val(xt, px0, k0, px1, k1);
		}
	    }
	}
    }

    return ret;
}

static NODE *isoconv_node (NODE *t, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *e, *n = t->v.b1.b;
    NODE *ret = NULL;
    const double *x = NULL;
    double *ymd[3] = {NULL, NULL, NULL};
    int i, k = n->v.bn.n_nodes;

    if (p->dset == NULL) {
	p->err = E_NODATA;
	return NULL;
    }

    if (k < 3 || k > 4) {
	n_args_error(k, 4, t->t, p);
    } else {
	/* evaluate the first (series) argument */
	e = eval(n->v.bn.n[0], p);
	if (!p->err && e->t != SERIES) {
	    node_type_error(t->t, 1, SERIES, e, p);
	} else {
	    x = e->v.xvec + p->dset->t1;
	}
    }

    for (i=1; i<k && !p->err; i++) {
	/* the remaining args must be addresses of series */
	e = n->v.bn.n[i];
	if (i == 3 && null_or_empty(e)) {
	    ; /* OK for the last one to be omitted */
	} else if (e->t != U_ADDR) {
	    node_type_error(t->t, i+1, U_ADDR, e, p);
	} else {
	    e = e->v.b1.b;
	    if (e->t != SERIES) {
		node_type_error(t->t, i+1, SERIES, e, p);
	    } else {
		ymd[i-1] = p->dset->Z[e->vnum] + p->dset->t1;
	    }
	}
    }

    if (!p->err) {
	reset_p_aux(p, save_aux);
	ret = aux_scalar_node(p);
    }

    if (!p->err) {
	int n = sample_size(p->dset);

	ret->v.xval = iso_basic_to_extended(x, ymd[0], ymd[1], ymd[2], n);
    }

    return ret;
}

static GretlType lh_array_type (parser *p)
{
    GretlType gtype = p->lh.gtype;

    if (gtype == GRETL_TYPE_NONE) {
	/* shouldn't happen? */
	if (p->lh.t == ARRAY) {
	    gretl_array *a = gen_get_lhs_var(p, GRETL_TYPE_ARRAY);

	    gtype = gretl_array_get_type(a);
	}
    }

    if (gtype != GRETL_TYPE_NONE && !gretl_is_array_type(gtype)) {
	p->err = E_TYPES;
    }

    return gtype;
}

static int check_array_element_type (NODE *n, GretlType *pt, int i)
{
    GretlType t = *pt;
    int ok = 0;

    if (t == GRETL_TYPE_MATRICES) {
	ok = n->t == MAT;
    } else if (t == GRETL_TYPE_STRINGS) {
	ok = n->t == STR;
    } else if (t == GRETL_TYPE_BUNDLES) {
	ok = n->t == BUNDLE;
    } else if (t == GRETL_TYPE_LISTS) {
	ok = n->t == LIST;
    } else if (i == 0 && t == GRETL_TYPE_ANY) {
	/* The array type is not yet determinate: this is OK
	   only when we're looking at the first element: if
	   this element is of an acceptable type we use it
	   to determine the array type.
	*/
	*pt = 0;
	if (n->t == MAT) {
	    *pt = GRETL_TYPE_MATRICES;
	} else if (n->t == STR) {
	    *pt = GRETL_TYPE_STRINGS;
	} else if (n->t == BUNDLE) {
	    *pt = GRETL_TYPE_BUNDLES;
	} else if (n->t == LIST) {
	    *pt = GRETL_TYPE_LISTS;
	}
	if (*pt != 0) {
	    ok = 1;
	}
    }

    return ok ? 0 : E_TYPES;
}

static void node_nullify_ptr (NODE *n)
{
    if      (n->t == MAT)    n->v.m = NULL;
    else if (n->t == STR)    n->v.str = NULL;
    else if (n->t == BUNDLE) n->v.b = NULL;
    else if (n->t == LIST)   n->v.ivec = NULL;
    else if (n->t == ARRAY)  n->v.a = NULL;
    else if (n->t == SERIES) n->v.xvec = NULL;
}

/* serves retrieval of data for candidate array elements
   or bundle members
*/

static void *node_get_ptr (NODE *n, int f, parser *p, int *donate)
{
    void *ptr = NULL;
    int t = n->t;

    /* default to copying the node's data */
    *donate = 0;

    if (f == F_DEFBUNDLE) {
	/* specific to bundles */
	if (t == ARRAY) {
	    ptr = n->v.a;
	} else if (t == SERIES) {
	    ptr = n->v.xvec;
	} else if (t == NUM) {
	    ptr = &n->v.xval;
	} else if (scalar_matrix_node(n)) {
	    ptr = n->v.m->val;
	    t = NUM;
	}
    }

    if (ptr == NULL) {
	/* common to arrays and bundles */
	if (t == MAT) {
	    ptr = n->v.m;
	} else if (t == STR) {
	    ptr = n->v.str;
	} else if (t == BUNDLE) {
	    ptr = n->v.b;
	} else if (t == LIST) {
	    ptr = n->v.ivec;
	}
    }

    if (t == NUM) {
	*donate = 1;
    } else if (!reusable(p) && is_tmp_node(n)) {
	*donate = 1;
	node_nullify_ptr(n);
    }

    return ptr;
}

/* given an FARGS node, detect if the first argument
   is a pointer to bundle */

static int bundle_pointer_arg0 (NODE *t)
{
    NODE *n = t->v.b1.b;

    if (n->v.bn.n_nodes > 0) {
	NODE *n0 = n->v.bn.n[0];

	if (n0->t == U_ADDR && ubundle_node(n0->v.b1.b)) {
	    return 1;
	}
    }

    return 0;
}

/* evaluate a built-in function that has more than three arguments */

static NODE *eval_nargs_func (NODE *t, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *n = t->v.b1.b;
    NODE *ret = NULL;
    NODE *e = NULL;
    int i, k = n->v.bn.n_nodes;

    if (t->t == F_BKFILT) {
	const double *x = NULL;
	int bk[3] = {0};

	if (k < 1 || k > 4) {
	    n_args_error(k, 4, t->t, p);
	}

	/* evaluate the first (series) argument */
	e = eval(n->v.bn.n[0], p);
	if (!p->err && e->t != SERIES) {
	    node_type_error(t->t, 1, SERIES, e, p);
	}

	if (!p->err) {
	    x = e->v.xvec;
	}

	for (i=1; i<k && !p->err; i++) {
	    e = n->v.bn.n[i];
	    if (e->t == EMPTY) {
		; /* NULL arguments are OK */
	    } else {
		e = eval(n->v.bn.n[i], p);
		if (e == NULL) {
		    fprintf(stderr, "eval_nargs_func: failed "
			    "to evaluate arg %d\n", i);
		} else {
		    bk[i-1] = node_get_int(e, p);
		}
	    }
	}

	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_series_node(p);
	}

	if (!p->err) {
	    p->err = bkbp_filter(x, ret->v.xvec, p->dset, bk[0], bk[1], bk[2]);
	}
    } else if (t->t == F_FILTER) {
	const double *x = NULL;
	gretl_matrix *X = NULL;
	gretl_matrix *C = NULL;
	gretl_matrix *A = NULL;
	int freeA = 0, freeC = 0;
	double y0 = 0;

	if (k < 1 || k > 4) {
	    n_args_error(k, 4, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (i == 0) {
		/* the series or matrix to filter */
		if (e->t != SERIES && e->t != MAT) {
		   node_type_error(t->t, i+1, 0, e, p);
		} else if (e->t == SERIES) {
		   x = e->v.xvec;
		} else {
		   X = e->v.m;
		}
	    } else if (i == 1) {
		/* matrix for MA polynomial (but we'll take a scalar) */
		if (e->t != MAT && e->t != NUM && e->t != EMPTY) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else if (e->t == MAT) {
		    C = e->v.m;
		} else if (e->t == NUM) {
		    C = gretl_matrix_from_scalar(e->v.xval);
		    if (C == NULL) {
			p->err = E_ALLOC;
		    } else {
			freeC = 1;
		    }
		}
	    } else if (i == 2) {
		/* matrix for AR polynomial (but we'll take a scalar) */
		if (e->t != MAT && e->t != NUM && e->t != EMPTY) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else if (e->t == MAT) {
		    A = e->v.m;
		} else if (e->t == NUM) {
		    A = gretl_matrix_from_scalar(e->v.xval);
		    if (A == NULL) {
			p->err = E_ALLOC;
		    } else {
			freeA = 1;
		    }
		}
	    } else if (i == 3) {
		/* initial (scalar) value for output series */
		if (!scalar_node(e)) {
		    node_type_error(t->t, i+1, NUM, e, p);
		} else {
		    y0 = node_get_scalar(e, p);
		    if (!p->err && na(y0)) {
			p->err = E_MISSDATA;
		    }
		}
	    }
	}

	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    if (X != NULL) {
		/* matrix output wanted */
		ret = aux_matrix_node(p);
		if (!p->err) {
		    ret->v.m = filter_matrix(X, A, C, y0, &(p->err));
		}
	    } else if (x != NULL) {
		/* series output */
		ret = aux_series_node(p);
		if (!p->err) {
		    p->err = filter_series(x, ret->v.xvec, p->dset, A, C, y0);
		}
	    }
	}

	if (freeA) gretl_matrix_free(A);
	if (freeC) gretl_matrix_free(C);

    } else if (t->t == F_MCOVG) {
	gretl_matrix *X = NULL;
	gretl_vector *u = NULL;
	gretl_vector *w = NULL;
	int targ, maxlag = 0;

	if (k != 4) {
	    n_args_error(k, 4, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    targ = (i == 3)? NUM : MAT;
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if ((i == 1 || i == 2) && e->t == EMPTY) {
		; /* for u or w, NULL is acceptable */
	    } else if (e->t != targ) {
		node_type_error(t->t, i+1, targ, e, p);
	    } else if (i == 0) {
		X = e->v.m;
	    } else if (i == 1) {
		u = e->v.m;
	    } else if (i == 2) {
		w = e->v.m;
	    } else if (i == 3) {
		maxlag = e->v.xval;
	    }
	}

	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_matrix_node(p);
	}

	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }
	    ret->v.m = gretl_matrix_covariogram(X, u, w, maxlag, &p->err);
	}
    } else if (t->t == F_MOLS || t->t == F_MPOLS) {
	gretlopt opt = (t->t == F_MPOLS)? OPT_M : OPT_NONE;
	gretl_matrix *Y = NULL;
	gretl_matrix *X = NULL;
	const char *SU = NULL;
	const char *SV = NULL;
	int freeY = 0, freeX = 0;

	if (k < 2 || k > 4) {
	    n_args_error(k, 1, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (p->err) {
		break;
	    }
	    if (i == 0) {
		if (e->t == SERIES) {
		    Y = tmp_matrix_from_series(e, p);
		    freeY = 1;
		} else if (e->t != MAT) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else {
		    Y = e->v.m;
		}
	    } else if (i == 1) {
		if (e->t == SERIES) {
		    X = tmp_matrix_from_series(e, p);
		    freeX = 1;
		} else if (e->t != MAT) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else {
		    X = e->v.m;
		}
	    } else {
		if (e->t == EMPTY) {
		    ; /* OK */
		} else if (e->t != U_ADDR) {
		    node_type_error(t->t, i+1, U_ADDR, e, p);
		} else if (i == 2) {
		    SU = ptr_node_get_matrix_name(e, p);
		} else {
		    SV = ptr_node_get_matrix_name(e, p);
		}
	    }
	}

	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_matrix_node(p);
	}

	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }
	    ret->v.m = user_matrix_ols(Y, X, SU, SV, opt, &p->err);
	}

	if (freeY) gretl_matrix_free(Y);
	if (freeX) gretl_matrix_free(X);
    } else if (t->t == F_MRLS) {
	gretl_matrix *Y = NULL;
	gretl_matrix *X = NULL;
	gretl_matrix *R = NULL;
	gretl_matrix *Q = NULL;
	const char *SU = NULL;
	const char *SV = NULL;

	if (k < 4 || k > 6) {
	    n_args_error(k, 1, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (p->err) {
		break;
	    }
	    if (i < 4) {
		if (e->t != MAT) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else if (i == 0) {
		    Y = e->v.m;
		} else if (i == 1) {
		    X = e->v.m;
		} else if (i == 2) {
		    R = e->v.m;
		} else if (i == 3) {
		    Q = e->v.m;
		}
	    } else {
		if (e->t == EMPTY) {
		    ; /* OK */
		} else if (e->t != U_ADDR) {
		    node_type_error(t->t, i+1, U_ADDR, e, p);
		} else if (i == 4) {
		    SU = ptr_node_get_matrix_name(e, p);
		} else {
		    SV = ptr_node_get_matrix_name(e, p);
		}
	    }
	}

	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_matrix_node(p);
	}

	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }
	    ret->v.m = user_matrix_rls(Y, X, R, Q, SU, SV, &p->err);
	}
    } else if (t->t == F_NRMAX) {
	gretl_matrix *b = NULL;
	const char *sf = NULL;
	const char *sg = NULL;
	const char *sh = NULL;

	if (k < 2 || k > 4) {
	    n_args_error(k, 4, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (p->err) {
		break;
	    }
	    if (i == 0) {
		if (e->t == U_ADDR) {
		    e = e->v.b1.b;
		}
		if (e->t != MAT) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else {
		    b = e->v.m;
		}
	    } else if (i == 1) {
		if (e->t != STR) {
		    node_type_error(t->t, i+1, STR, e, p);
		} else {
		    sf = e->v.str;
		}
	    } else if (e->t == EMPTY) {
		; /* OK */
	    } else if (e->t != STR) {
		node_type_error(t->t, i+1, STR, e, p);
	    } else if (i == 2) {
		sg = e->v.str;
	    } else {
		sh = e->v.str;
	    }
	}

	if (!p->err) {
	    if (!gretl_vector_get_length(b)) {
		p->err = E_TYPES;
	    } else if (!is_function_call(sf) ||
		       (sg != NULL && !is_function_call(sg)) ||
		       (sh != NULL && !is_function_call(sh))) {
		p->err = E_TYPES;
	    }
	}

	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_scalar_node(p);
	}

	if (!p->err) {
	    int minimize = (t->flags & ALS_NODE)? 1 : 0;

	    ret->v.xval = user_NR(b, sf, sg, sh, p->dset,
				  minimize, p->prn, &p->err);
	}
    } else if (t->t == F_LOESS) {
	const double *y = NULL, *x = NULL;
	double bandwidth = 0.5;
	int poly_order = 1;
	gretlopt opt = OPT_NONE;

	if (k < 2 || k > 6) {
	    n_args_error(k, 5, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (p->err) {
		break;
	    }
	    if (i < 2) {
		if (e->t != SERIES) {
		    node_type_error(t->t, i+1, SERIES, e, p);
		} else if (i == 0) {
		    y = e->v.xvec;
		} else {
		    x = e->v.xvec;
		}
	    } else if (i == 2 || i == 3) {
		if (e->t != NUM && e->t != EMPTY) {
		    node_type_error(t->t, i+1, NUM, e, p);
		} else if (i == 2 && e->t == NUM) {
		    poly_order = node_get_int(e, p);
		} else if (e->t == NUM) {
		    bandwidth = e->v.xval;
		}
	    } else {
		if (e->t != EMPTY && e->t != NUM) {
		    node_type_error(t->t, i+1, NUM, e, p);
		} else {
		    int ival = node_get_int(e, p);

		    if (!p->err && ival != 0) {
			if (i == 4) {
			    opt |= OPT_R;
			} else {
			    opt |= OPT_O;
			}
		    }
		}
	    }
	}
	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_series_node(p);
	    if (ret != NULL) {
		p->err = gretl_loess(y, x, poly_order, bandwidth,
				     opt, p->dset, ret->v.xvec);
	    }
	}
    } else if (t->t == F_GHK) {
	gretl_matrix *C = NULL;
	gretl_matrix *A = NULL;
	gretl_matrix *B = NULL;
	gretl_matrix *U = NULL;
	const char *dP_name = NULL;

	if (k < 4 || k > 5) {
	    n_args_error(k, 5, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (i < 4 && e->t != MAT) {
		node_type_error(t->t, i+1, MAT, e, p);
	    } else if (i == 0) {
		C = e->v.m;
	    } else if (i == 1) {
		A = e->v.m;
	    } else if (i == 2) {
		B = e->v.m;
	    } else if (i == 3) {
		U = e->v.m;
	    } else {
		/* the optional last argument */
		if (e->t == EMPTY) {
		    ; /* OK */
		} else if (e->t != U_ADDR) {
		    node_type_error(t->t, i+1, U_ADDR, e, p);
		} else {
		    dP_name = ptr_node_get_matrix_name(e, p);
		}
	    }
	}
	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_matrix_node(p);
	}
	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }
	    if (dP_name == NULL) {
		ret->v.m = gretl_GHK(C, A, B, U, &p->err);
	    } else {
		ret->v.m = user_matrix_GHK(C, A, B, U, dP_name, &p->err);
	    }
	}
    } else if (t->t == F_QUADTAB) {
	int order = -1, method = 1;
	double a = NADBL;
	double b = NADBL;

	if (k < 1 || k > 4) {
	    n_args_error(k, 4, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (i == 0) {
		order = node_get_int(e, p);
	    } else if (!empty_or_num(e)) {
		node_type_error(t->t, i+1, NUM, e, p);
	    } else if (i == 1) {
		method = node_get_int(e, p);
	    } else if (i == 2) {
		a = e->v.xval;
	    } else {
		b = e->v.xval;
	    }
	}
	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_matrix_node(p);
	}
	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }
	    ret->v.m = gretl_quadrule_matrix_new(order, method,
						 a, b, &p->err);
	}
    } else if (t->t == F_QLRPVAL) {
	double X2 = NADBL;
	double p1 = 0, p2 = 0;
	int df = 0;

	if (k != 4) {
	    n_args_error(k, 4, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (i == 0) {
		X2 = node_get_scalar(e, p);
	    } else if (i == 1) {
		df = node_get_int(e, p);
	    } else if (i == 2) {
		p1 = node_get_scalar(e, p);
	    } else {
		p2 = node_get_scalar(e, p);
	    }
	}
	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_scalar_node(p);
	}
	if (!p->err) {
	    ret->v.xval = QLR_pval(X2, df, p1, p2);
	}
    } else if (t->t == F_BOOTCI) {
	int cnum = -1, method = 1, B = 0;
	int studentize = 1;
	double alpha = NADBL;

	if (k < 1 || k > 5) {
	    n_args_error(k, 5, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (i == 0) {
		cnum = node_get_int(e, p);
	    } else if (!empty_or_num(e)) {
		node_type_error(t->t, i+1, NUM, e, p);
	    } else if (i == 1) {
		B = node_get_int(e, p);
	    } else if (i == 2) {
		alpha = e->v.xval;
	    } else if (i == 3) {
		method = node_get_int(e, p);
	    } else {
		studentize = node_get_int(e, p);
	    }
	}
	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_matrix_node(p);
	}
	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }
	    ret->v.m = last_model_get_boot_ci(cnum, p->dset, B, alpha, method,
					      studentize, &p->err);
	}
    } else if (t->t == F_BOOTPVAL) {
	int cnum = -1, method = 1, B = 0;

	if (k < 1 || k > 3) {
	    n_args_error(k, 3, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (i == 0) {
		cnum = node_get_int(e, p);
	    } else if (!empty_or_num(e)) {
		node_type_error(t->t, i+1, NUM, e, p);
	    } else if (i == 1) {
		B = node_get_int(e, p);
	    } else {
		method = node_get_int(e, p);
	    }
	}
	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_scalar_node(p);
	}
	if (!p->err) {
	    ret->v.xval = last_model_get_boot_pval(cnum, p->dset, B,
						   method, &p->err);
	}
    } else if (t->t == F_MOVAVG) {
	const double *x = NULL;
	double d = 0, y0 = NADBL;
	int len = 0, ctrl = -9999;
	int EMA = 0;

	if (k < 2 || k > 4) {
	    n_args_error(k, 4, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    if (i > 1 && null_or_empty(n->v.bn.n[i])) {
		continue; /* OK */
	    }
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (i == 0) {
		if (e->t == SERIES) {
		    x = e->v.xvec;
		} else {
		    node_type_error(t->t, i+1, SERIES, e, p);
		}
	    } else if (i == 1) {
		d = node_get_scalar(e, p);
		if (d < 1.0 && d > 0.0) {
		    EMA = 1;
		} else if (d < 1.0) {
		    p->err = E_INVARG;
		} else {
		    len = node_get_int(e, p);
		}
	    } else if (i == 2) {
		ctrl = node_get_int(e, p);
	    } else {
		y0 = node_get_scalar(e, p);
	    }
	}
	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_series_node(p);
	}
	if (!p->err) {
	    if (ctrl == -9999) {
		/* set the respective defaults */
		ctrl = EMA ? 1 : 0;
	    }
	    if (EMA) {
		p->err = exponential_movavg_series(x, ret->v.xvec, p->dset,
						   d, ctrl, y0);
	    } else {
		p->err = movavg_series(x, ret->v.xvec, p->dset, len, ctrl);
	    }
	}
    } else if (t->t == HF_CLOGFI) {
	const char *dfname = NULL;
	gretl_matrix *z = NULL;
	int T = 0, K = 0;

	if (k < 3 || k > 4) {
	    n_args_error(k, 4, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (i == 0) {
		if (scalar_node(e)) {
		    T = node_get_int(e, p);
		} else {
		    node_type_error(t->t, 1, NUM, e, p);
		}
	    } else if (i == 1) {
		if (scalar_node(e)) {
		    K = node_get_int(e, p);
		} else {
		    node_type_error(t->t, 2, NUM, e, p);
		}
	    } else if (i == 2) {
		if (e->t == MAT) {
		    z = e->v.m;
		} else {
		    node_type_error(t->t, 3, MAT, e, p);
		}
	    } else if (i == 3) {
		/* optional matrix-pointer for storing the
		   derivative wrt z */
		if (e->t == EMPTY) {
		    ; /* OK */
		} else if (e->t != U_ADDR) {
		    node_type_error(t->t, 4, U_ADDR, e, p);
		} else {
		    dfname = ptr_node_get_matrix_name(e, p);
		}
	    }
	}

	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_scalar_node(p);
	}

	if (!p->err) {
	    ret->v.xval = clogit_fi(T, K, z, dfname, &p->err);
	}
    } else if (t->t == F_DEFARRAY) {
	GretlType gtype = lh_array_type(p);
	gretl_array *A = NULL;
	int anon = 0;
	void *ptr;

	if (gtype == GRETL_TYPE_NONE) {
	    /* type remains to be determined */
	    gtype = GRETL_TYPE_ANY;
	    anon = 1;
	}

	if (!p->err) {
	    A = gretl_array_new(gtype, 0, &p->err);
	}

	if (!p->err) {
	    for (i=0; i<k && !p->err; i++) {
		int donate = 0;

		e = eval(n->v.bn.n[i], p);
		if (!p->err) {
		    p->err = check_array_element_type(e, &gtype, i);
		}
		if (!p->err) {
		    if (anon) {
			gretl_array_set_type(A, gtype);
			anon = 0;
		    }
		    ptr = node_get_ptr(e, t->t, p, &donate);
		    if (donate) {
			/* copy not required */
			p->err = gretl_array_append_object(A, ptr, 0);
		    } else {
			p->err = gretl_array_append_object(A, ptr, 1);
		    }
		}
	    }
	}

	if (p->err) {
	    gretl_array_destroy(A);
	} else {
	    reset_p_aux(p, save_aux);
	    ret = aux_array_node(p);
	    if (ret != NULL) {
		ret->v.a = A;
	    }
	}
    } else if (t->t == F_DEFBUNDLE) {
	gretl_bundle *b = NULL;
	GretlType gtype;
	char *key = NULL;
	void *ptr;
	int donate;

	if (k % 2 != 0) {
	    p->err = E_PARSE;
	} else {
	    b = gretl_bundle_new();
	    if (b == NULL) {
		p->err = E_ALLOC;
	    }
	}

	if (!p->err) {
	    for (i=0; i<k && !p->err; i++) {
		ptr = NULL;
		e = eval(n->v.bn.n[i], p);
		if (p->err) {
		    break;
		}
		if (i == 0 || i % 2 == 0) {
		    /* we need a key string */
		    if (e->t == STR) {
			key = e->v.str;
		    } else {
			p->err = E_TYPES;
		    }
		} else {
		    /* we need some valid content */
		    gtype = gretl_type_from_gen_type(e->t);
		    if (type_can_be_bundled(gtype)) {
			int size = 0;

			if (e->t == SERIES) {
			    size = p->dset->n;
			} else if (scalar_matrix_node(e)) {
			    gtype = GRETL_TYPE_DOUBLE;
			}
			ptr = node_get_ptr(e, t->t, p, &donate);
			if (donate) {
			    gretl_bundle_donate_data(b, key, ptr, gtype, size);
			} else {
			    gretl_bundle_set_data(b, key, ptr, gtype, size);
			}
		    } else {
			p->err = E_TYPES;
		    }
		}
	    }
	}

	if (p->err) {
	    gretl_bundle_destroy(b);
	} else {
	    reset_p_aux(p, save_aux);
	    ret = aux_bundle_node(p);
	    if (ret != NULL) {
		ret->v.b = b;
	    }
	}
    } else if (t->t == F_DEFLIST) {
	int *li, *full_list = gretl_list_new(0);

	for (i=0; i<k && !p->err; i++) {
	    li = NULL;
	    e = eval(n->v.bn.n[i], p);
	    if (!p->err) {
		if (k == 1 && e->t == ARRAY) {
		    li = list_from_strings_array(e->v.a, p);
		} else if (ok_list_node(e)) {
		    li = node_get_list(e, p);
		} else {
		    p->err = E_TYPES;
		}
	    }
	    if (!p->err && li[0] > 0) {
		gretl_list_append_list(&full_list, li, &p->err);
	    }
	    free(li);
	}

	if (p->err) {
	    free(full_list);
	} else {
	    reset_p_aux(p, save_aux);
	    ret = aux_list_node(p);
	    if (ret != NULL) {
		ret->v.ivec = full_list;
	    }
	}
    }

    return ret;
}

static gretl_matrix *node_get_matrix_lenient (NODE *n,
					      int ok,
					      parser *p)
{
    gretl_matrix *m = NULL;

    if (n->t == NUM && (ok == NUM)) {
	m = gretl_matrix_from_scalar(n->v.xval);
    } else if (n->t == SERIES && (ok == SERIES)) {
	m = gretl_vector_from_series(n->v.xvec, p->dset->t1,
				     p->dset->t2);
    } else if (n->t == LIST && (ok == SERIES)) {
	m = gretl_matrix_data_subset(n->v.ivec, p->dset,
				     p->dset->t1, p->dset->t2,
				     M_MISSING_OK, &p->err);
    } else {
	p->err = E_TYPES;
    }

    if (!p->err && m == NULL) {
	p->err = E_ALLOC;
    }

    return m;
}

static gretl_bundle *get_kalman_bundle_arg (NODE *n, parser *p)
{
    gretl_bundle *b = NULL;
    NODE *e = NULL;

    if (n->t == FARGS) {
	/* multi-arguments node */
	e = n->v.bn.n[0];
	e = e->v.b1.b;
    } else if (n->t == U_ADDR) {
	e = n->v.b1.b;
    } else if (n->t == BUNDLE) {
	e = n;
    }

    if (e == NULL || e->t != BUNDLE) {
	p->err = E_TYPES;
    } else {
	b = e->v.b; /* get the actual bundle */
	if (gretl_bundle_get_type(b) != BUNDLE_KALMAN ||
	    gretl_bundle_get_private_data(b) == NULL) {
	    p->err = E_TYPES;
	    b = NULL;
	}
    }

    if (p->err) {
	gretl_errmsg_set("Argument 1 must point to a state-space bundle");
    }

    return b;
}

static NODE *eval_kalman_bundle_func (NODE *t, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *n = t->v.b1.b;
    NODE *ret = NULL;
    NODE *e = NULL;
    int i, k = n->v.bn.n_nodes;

    if (t->t == F_KSETUP) {
	gretl_matrix *M[5] = {NULL};
	int copy[5] = {0};

	if (k < 4) {
	    n_args_error(k, 4, t->t, p);
	}

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (!p->err) {
		if (e->t == MAT) {
		    M[i] = e->v.m;
		    if (is_tmp_node(e)) {
			e->v.m = NULL;
		    } else {
			copy[i] = 1;
		    }
		} else if (i == 0) {
		    /* obsy: accept series or list */
		    M[i] = node_get_matrix_lenient(e, SERIES, p);
		} else {
		    /* system matrices, state variance */
		    M[i] = node_get_matrix_lenient(e, NUM, p);
		}
	    }
	}

	if (!p->err) {
	    gretl_bundle *b = kalman_bundle_new(M, copy, k, &p->err);

	    if (!p->err) {
		ret = aux_bundle_node(p);
		if (ret != NULL) {
		    ret->v.b = b;
		}
	    }
	}
    } else if (t->t == F_KFILTER) {
	gretl_bundle *b = get_kalman_bundle_arg(n, p);

	if (!p->err && k != 1) {
	    n_args_error(k, 1, t->t, p);
	}
	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_scalar_node(p);
	}
	if (!p->err) {
	    ret->v.xval = kalman_bundle_run(b, p->prn, &p->err);
	}
    } else if (t->t == F_KDSMOOTH) {
	gretl_bundle *b = get_kalman_bundle_arg(n, p);
	int param = 1;
	int dkstyle = 0;

	if (!p->err) {
	    if (k == 2) {
		e = eval(n->v.bn.n[1], p);
		dkstyle = node_get_int(e, p);
	    } else if (k < 1 || k > 2) {
		n_args_error(k, 2, t->t, p);
	    }
	}
	if (!p->err) {
	    param += dkstyle != 0;
	    reset_p_aux(p, save_aux);
	    ret = aux_scalar_node(p);
	    if (!p->err) {
		ret->v.xval = kalman_bundle_smooth(b, param, p->prn);
	    }
	}
    } else if (t->t == F_KSMOOTH) {
	gretl_bundle *b = get_kalman_bundle_arg(n, p);

	if (!p->err && k != 1) {
	    n_args_error(k, 1, t->t, p);
	}
	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_scalar_node(p);
	    if (!p->err) {
		ret->v.xval = kalman_bundle_smooth(b, 0, p->prn);
	    }
	}
    } else if (t->t == F_KSIMUL) {
	/* we need a bundle pointer, a matrix,
	   and perhaps an optional boolean
	*/
	gretl_bundle *b = get_kalman_bundle_arg(n, p);
	gretl_matrix *U = NULL;
	int freeU = 0;
	int get_state = 0;

	if (!p->err && k != 2 && k != 3) {
	    n_args_error(k, 2, t->t, p);
	}

	for (i=1; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (!p->err && i == 1) {
		if (e->t == MAT) {
		    U = e->v.m;
		} else {
		    U = node_get_matrix_lenient(e, SERIES, p);
		    if (U != NULL) {
			freeU = 1;
		    }
		}
	    } else if (!p->err) {
		get_state = node_get_int(e, p);
	    }
	}

	if (!p->err) {
	    reset_p_aux(p, save_aux);
	    ret = aux_matrix_node(p);
	}

	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }
	    ret->v.m = kalman_bundle_simulate(b, U, get_state,
					      p->prn, &p->err);
	}

	if (freeU) gretl_matrix_free(U);
    }

    return ret;
}

static NODE *kalman_data_node (NODE *l, NODE *r, parser *p)
{
    NODE *save_aux = p->aux;
    gretl_bundle *b = get_kalman_bundle_arg(l, p);
    NODE *ret = NULL;

    if (!p->err) {
	reset_p_aux(p, save_aux);
	ret = aux_matrix_node(p);
    }

    if (!p->err) {
	if (ret->v.m != NULL) {
	    gretl_matrix_free(ret->v.m);
	}
	ret->v.m = kalman_bundle_simdata(b, r->v.m, p->prn, &p->err);
    }

    return ret;
}

/* Create a matrix using selected series, or a mixture of series and
   lists, or more than one list.  We proceed by setting up a "dummy"
   dataset and constructing a list that indexes into it.  (We can't
   use a regular list, in the general case, since some of the series
   may be temporary variables that are not part of the "real"
   dataset.)
*/

static gretl_matrix *assemble_matrix (GPtrArray *a, int nnodes, parser *p)
{
    NODE *n;
    gretl_matrix *m = NULL;
    const int *list;
    double **Z = NULL;
    int *dumlist;
    int i, j, k = 0;

#if EDEBUG
    fprintf(stderr, "assemble_matrix...\n");
#endif

    /* how many columns will we need? */
    for (i=0; i<nnodes; i++) {
	n = g_ptr_array_index(a, i);
	if (n->t == LIST) {
	    k += n->v.ivec[0];
	} else if (n->t == SERIES) {
	    k++;
	}
    }

    /* create dummy data array */
    Z = malloc(k * sizeof *Z);
    if (Z == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

#if EDEBUG
    fprintf(stderr, " got %d columns, Z at %p\n", k, (void *) Z);
#endif

    /* and a list associated with Z */
    dumlist = gretl_consecutive_list_new(0, k-1);
    if (dumlist == NULL) {
	p->err = E_ALLOC;
	free(Z);
	return NULL;
    }

    /* attach series pointers to Z */
    k = 0;
    for (i=0; i<nnodes; i++) {
	n = g_ptr_array_index(a, i);
	if (n->t == LIST) {
	    list = n->v.ivec;
	    for (j=1; j<=list[0]; j++) {
		Z[k++] = p->dset->Z[list[j]];
	    }
	} else if (n->t == SERIES) {
	    Z[k++] = n->v.xvec;
	}
    }

    if (!p->err) {
	DATASET dumset = {0};

	dumset.Z = Z;
	dumset.v = k;
	dumset.n = p->dset->n;
	dumset.t1 = p->dset->t1;
	dumset.t2 = p->dset->t2;

	m = real_matrix_from_list(dumlist, &dumset, p);
    }

    free(dumlist);
    free(Z);

    return m;
}

#define ok_matdef_sym(s) (s == NUM || s == SERIES || s == EMPTY || \
                          s == DUM || s == LIST)

/* composing a matrix from scalars, series or lists */

static NODE *matrix_def_node (NODE *nn, parser *p)
{
    GPtrArray *a;
    gretl_matrix *M = NULL;
    NODE *save_aux = p->aux;
    NODE *n, *ret = NULL;
    int k = nn->v.bn.n_nodes;
    int nnum = 0, nvec = 0;
    int dum = 0, nsep = 0;
    int nlist = 0;
    int seppos = -1;
    int i;

    if (autoreg(p)) {
	fprintf(stderr, "You can't define a matrix in this context\n");
	p->err = E_TYPES;
	return NULL;
    }

#if EDEBUG
    fprintf(stderr, "Processing MDEF...\n");
#endif

    a = g_ptr_array_sized_new(k);

    for (i=0; i<k && !p->err; i++) {
	n = eval(nn->v.bn.n[i], p);
	if (n == NULL && !p->err) {
	    p->err = E_UNSPEC; /* "can't happen" */
	}
	if (p->err) {
	    break;
	}
	if (!ok_matdef_sym(n->t) && !scalar_matrix_node(n)) {
	    fprintf(stderr, "matrix_def_node: node type %d: not OK\n", n->t);
	    p->err = E_TYPES;
	    break;
	}
	if (scalar_node(n)) {
	    nnum++;
	} else if (n->t == SERIES) {
	    nvec++;
	} else if (n->t == DUM) {
	    dum++;
	} else if (n->t == LIST) {
	    nlist++;
	} else if (n->t == EMPTY) {
	    if (nsep == 0) {
		seppos = i;
	    }
	    nsep++;
	}
	if (dum && k != 1) {
	    /* dummy must be singleton node */
	    p->err = E_TYPES;
	} else if ((nvec || nlist) && nnum) {
	    /* can't mix series/lists with scalars */
	    p->err = E_TYPES;
	} else if ((nvec || nlist) && nsep) {
	    /* can't have row separators in a matrix
	       composed of series or lists */
	    p->err = E_TYPES;
	}
	if (!p->err) {
	    g_ptr_array_add(a, n);
	}
    }

    if (!p->err) {
	if (nvec > 0 || nlist > 1) {
	    M = assemble_matrix(a, k, p);
	} else if (nnum > 0) {
	    M = matrix_from_scalars(a, k, nsep, seppos, p);
	} else if (nlist) {
	    M = matrix_from_list(g_ptr_array_index(a, 0), p);
	} else if (dum) {
	    n = g_ptr_array_index(a, 0);
	    if (n->v.idnum == DUM_DATASET) {
		M = matrix_from_list(NULL, p);
	    } else {
		pprintf(p->prn, "Wrong sort of dummy var\n");
		p->err = E_TYPES;
	    }
	} else {
	    /* empty matrix def */
	    M = gretl_null_matrix_new();
	}
    }

    if (a != NULL) {
	g_ptr_array_free(a, TRUE);
    }

    if (p->err) {
	if (M != NULL) {
	    gretl_matrix_free(M);
	}
    } else {
	reset_p_aux(p, save_aux);
	ret = aux_matrix_node(p);
	if (ret != NULL) {
	    ret->v.m = M;
	} else {
	    gretl_matrix_free(M);
	}
    }

    return ret;
}

static NODE *gen_series_from_string (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;
    gchar *line;
    int vnum = -1;
    int err = 0;

    line = g_strdup_printf("%s=%s", l->v.str, r->v.str);
    err = generate(line, p->dset, GRETL_TYPE_SERIES,
		   OPT_NONE, p->prn);

    if (!err) {
	vnum = current_series_index(p->dset, l->v.str);
    }

    ret = aux_scalar_node(p);
    if (ret != NULL) {
	ret->v.xval = vnum;
    }

    g_free(line);

    return ret;
}

static double *xvec_from_matrix (gretl_matrix *m,
				 parser *p,
				 int *subset,
				 int *err)
{
    double *ret = NULL;

    if (gretl_is_null_matrix(m) || m->cols != 1) {
	*err = E_TYPES;
    } else if (m->rows == p->dset->n) {
	ret = m->val;
    } else if (m->rows == sample_size(p->dset)) {
	*subset = 1;
	ret = m->val;
    } else {
	*err = E_TYPES;
    }

    return ret;
}

static NODE *gen_series_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (p->dset == NULL || p->dset->n == 0) {
	no_data_error(p);
    } else if (l->t == STR && r->t == STR) {
	return gen_series_from_string(l, r, p);
    } else if (l->t != STR || (r->t != SERIES && r->t != MAT)) {
	p->err = E_TYPES;
    } else {
	char *vname = l->v.str;
	int vnum = current_series_index(p->dset, vname);
	int subset = 0;
	double *xvec;
	int err = 0;

	if (r->t == SERIES) {
	    xvec = r->v.xvec;
	} else {
	    xvec = xvec_from_matrix(r->v.m, p, &subset, &err);
	}

	if (!err && vnum > 0) {
	    /* a series of this name already exists */
	    int t, i = 0;

	    for (t=p->dset->t1; t<=p->dset->t2; t++) {
		if (subset) {
		    p->dset->Z[vnum][t] = xvec[i++];
		} else {
		    p->dset->Z[vnum][t] = xvec[t];
		}
	    }
	} else if (!err) {
	    /* creating a new series */
	    GretlType ltype = user_var_get_type_by_name(vname);

	    if (ltype != GRETL_TYPE_NONE) {
		/* cannot overwrite a variable of another type */
		err = E_TYPES;
	    } else {
		err = check_varname(vname);
	    }
	    if (!err) {
		err = dataset_add_NA_series(p->dset, 1);
	    }
	    if (!err) {
		int t, i = 0, v = p->dset->v - 1;

		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    if (subset) {
			p->dset->Z[v][t] = xvec[i++];
		    } else {
			p->dset->Z[v][t] = xvec[t];
		    }
		}
	    }
	    if (!err) {
		vnum =  p->dset->v - 1;
		strcpy(p->dset->varname[vnum], vname);
	    }
	}

	ret = aux_scalar_node(p);
	if (ret != NULL) {
	    ret->v.xval = err ? -1 : vnum;
	}
    }

    return ret;
}

static NODE *gen_array_node (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (!empty_or_num(n)) {
	p->err = E_TYPES;
    } else if (p->lh.gtype == 0) {
	gretl_errmsg_set(_("array: no type was specified"));
	p->err = E_DATA;
    } else {
	int len = 0;

	if (n != NULL && n->t == NUM) {
	    len = node_get_int(n, p);
	}

	if (!p->err) {
	    ret = aux_array_node(p);
	    if (!p->err) {
		ret->v.a = gretl_array_new(p->lh.gtype, len, &p->err);
	    }
	}
    }

    return ret;
}

static NODE *get_series_stringvals (NODE *n, parser *p)
{
    NODE *ret = aux_array_node(p);

    if (!p->err) {
	int v = n->vnum;

	if (is_string_valued(p->dset, v)) {
	    int n_strs = 0;
	    char **S;

	    S = series_get_string_vals(p->dset, v, &n_strs);
	    ret->v.a = gretl_array_from_strings(S, n_strs, 1,
						&p->err);
	} else {
	    ret->v.a = gretl_array_new(GRETL_TYPE_STRINGS, 0, &p->err);
	}
    }

    return ret;
}

static NODE *stringify_series (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	gretl_array *A = r->v.a;
	char **S;
	int ns = 0;

	S = gretl_array_get_strings(A, &ns);

	if (S == NULL) {
	    p->err = E_DATA;
	} else {
	    /* flag an error on error? */
	    ret->v.xval = series_set_string_vals(p->dset, l->vnum, S, ns);
	    if (ret->v.xval != 0.0) {
		/* for now, yes */
		p->err = E_DATA;
	    }
	}
    }

    return ret;
}

enum {
    FORK_L,
    FORK_R,
    FORK_BOTH,
    FORK_NONE
};

/* Determine whether or not a series is constant in boolean terms,
   i.e. all elements zero, or all non-zero, over the relevant range.
   If so, return FORK_L (all 1) or FORK_R (all 0), othewise
   return FORK_UNK.
*/

static int vec_branch (const double *c, parser *p)
{
    int c1, t, t1, t2;
    int ret;

    t1 = (autoreg(p))? p->obs : p->dset->t1;
    t2 = (autoreg(p))? p->obs : p->dset->t2;

    c1 = (c[t1] != 0.0);
    ret = (c1)? FORK_L : FORK_R;

    for (t=t1; t<=t2; t++) {
	if (!xna(c[t])) {
	    if ((c1 && c[t] == 0) || (!c1 && c[t] != 0)) {
		ret = FORK_BOTH;
		break;
	    }
	}
    }

    return ret;
}

/* Given a series condition in a ternary "?" expression, return the
   evaluated counterpart.  We evaluate both forks and select based on
   the value of the condition at each observation.  We accept only
   scalar (NUM or 1x1 matrix) and series (SERIES) types on input, and
   always produce a SERIES type on output.
*/

static NODE *query_eval_series (const double *c, NODE *n, parser *p)
{
    NODE *l = NULL, *r = NULL, *ret = NULL;
    NODE *save_aux = p->aux;
    double *xvec = NULL, *yvec = NULL;
    double x = NADBL, y = NADBL;
    double xt, yt;
    int t, t1, t2;
    int branch;

    branch = vec_branch(c, p);

    if (autoreg(p) || branch != FORK_R) {
	l = eval(n->v.b3.m, p);
	if (p->err) {
	    return NULL;
	}
	if (l->t == SERIES) {
	    xvec = l->v.xvec;
	} else if (scalar_node(l)) {
	    x = node_get_scalar(l, p);
	} else {
	    p->err = E_TYPES;
	    return NULL;
	}
    }

    if (autoreg(p) || branch != FORK_L) {
	r = eval(n->v.b3.r, p);
	if (p->err) {
	    return NULL;
	}
	if (r->t == SERIES) {
	    yvec = r->v.xvec;
	} else if (scalar_node(r)) {
	    y = node_get_scalar(r, p);
	} else {
	    p->err = E_TYPES;
	    return NULL;
	}
    }

    reset_p_aux(p, save_aux);
    ret = aux_series_node(p);

    t1 = (autoreg(p))? p->obs : p->dset->t1;
    t2 = (autoreg(p))? p->obs : p->dset->t2;

    for (t=t1; t<=t2; t++) {
	if (xna(c[t])) {
	    ret->v.xvec[t] = NADBL;
	} else {
	    xt = (xvec != NULL)? xvec[t] : x;
	    yt = (yvec != NULL)? yvec[t] : y;
	    ret->v.xvec[t] = (c[t] != 0.0)? xt : yt;
	}
    }

    return ret;
}

/* The condition in the ternary query operator is a scalar,
   which has been evaluated to @x, and which must now be
   interpreted as a boolean.
*/

static NODE *query_eval_scalar (double x, NODE *n, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *l = NULL, *r = NULL, *ret = NULL;
    int indef = na(x) || isnan(x);
    int branch;

    branch = indef ? FORK_NONE : (x != 0 ? FORK_L : FORK_R);

    if (autoreg(p) || branch != FORK_R) {
	l = eval(n->v.b3.m, p);
	if (p->err) {
	    return NULL;
	}
    }

    if (autoreg(p) || branch != FORK_L) {
	r = eval(n->v.b3.r, p);
	if (p->err) {
	    return NULL;
	}
    }

    if (branch == FORK_NONE) {
	reset_p_aux(p, save_aux);
	ret = aux_scalar_node(p);
	if (ret != NULL) {
	    ret->v.xval = NADBL;
	}
    } else if (branch == FORK_L) {
	ret = l;
    } else if (branch == FORK_R) {
	ret = r;
    }

    return ret;
}

/* the condition in the ternary query operator is a matrix */

static NODE *query_eval_matrix (gretl_matrix *m, NODE *n, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *ret, *l, *r;
    gretl_matrix *mret;

    if (gretl_is_null_matrix(m)) {
	p->err = E_TYPES;
	return NULL;
    }

    l = eval(n->v.b3.m, p);

    if (!p->err) {
	r = eval(n->v.b3.r, p);
    }

    if (p->err) {
	return NULL;
    }

    if ((l->t != NUM && l->t != MAT) ||
	(r->t != NUM && r->t != MAT)) {
	p->err = E_TYPES;
	return NULL;
    }

    if (l->t == MAT && (l->v.m->cols != m->cols ||
			l->v.m->rows != m->rows)) {
	p->err = E_NONCONF;
	return NULL;
    } else if (r->t == MAT && (r->v.m->cols != m->cols ||
			       r->v.m->rows != m->rows)) {
	p->err = E_NONCONF;
	return NULL;
    }

    mret = gretl_matrix_copy(m);
    if (mret == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

    reset_p_aux(p, save_aux);
    ret = aux_matrix_node(p);

    if (!p->err) {
	double x, y;
	int j, i;

	for (j=0; j<m->cols; j++) {
	    for (i=0; i<m->rows; i++) {
		x = gretl_matrix_get(m, i, j);
		if (na(x)) x = M_NA;
		if (isnan(x)) {
		    gretl_matrix_set(mret, i, j, x);
		} else if (x != 0.0) {
		    y = l->t == NUM ? l->v.xval :
			gretl_matrix_get(l->v.m, i, j);
		    if (na(y)) y = M_NA;
		    gretl_matrix_set(mret, i, j, y);
		} else {
		    y = r->t == NUM ? r->v.xval :
			gretl_matrix_get(r->v.m, i, j);
		    if (na(y)) y = M_NA;
		    gretl_matrix_set(mret, i, j, y);
		}
	    }
	}
	ret->v.m = mret;
    }

    return ret;
}

/* Evaluate a ternary "query" expression: C ? X : Y.  The condition C
   may be a scalar, series or matrix.  The relevant sub-nodes of @t
   are named "l" (left, the condition), "m" and "r" (middle and right
   respectively, the two alternates).
*/

static NODE *eval_query (NODE *t, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *c, *ret = NULL;

    /* evaluate and check the condition */
    c = eval(t->v.b3.l, p);

#if EDEBUG
    fprintf(stderr, "eval_query: t=%p, l=%p, m=%p, r=%p\n",
	    (void *) t, (void *) t->v.b3.l, (void *) t->v.b3.m,
	    (void *) t->v.b3.r);
    if (c->t == NUM) {
	fprintf(stderr, " condition type=NUM, value=%g\n", c->v.xval);
    } else {
	fprintf(stderr, " condition type=%s\n", getsymb(c->t));
    }
#endif

    if (!p->err) {
	reset_p_aux(p, save_aux);
	if (c->t == NUM) {
	    ret = query_eval_scalar(c->v.xval, t, p);
	} else if (c->t == SERIES) {
	    ret = query_eval_series(c->v.xvec, t, p);
	} else if (c->t == MAT) {
	    if (gretl_matrix_is_scalar(c->v.m)) {
		ret = query_eval_scalar(c->v.m->val[0], t, p);
	    } else {
		ret = query_eval_matrix(c->v.m, t, p);
	    }
	} else {
	    /* invalid type for boolean condition */
	    p->err = E_TYPES;
	}
    }

#if EDEBUG
    fprintf(stderr, "eval_query return: ret = %p\n", (void *) ret);
#endif

    return ret;
}

#define dvar_scalar(i) (i > 0 && i < R_SCALAR_MAX)
#define dvar_series(i) (i > R_SCALAR_MAX && i < R_SERIES_MAX)
#define dvar_variant(i) (i > R_SERIES_MAX && i < R_MAX)

#define no_data(p) (p == NULL || p->n == 0)

double dvar_get_scalar (int i, const DATASET *dset)
{
    switch (i) {
    case R_NOBS:
	return (dset == NULL) ? NADBL :
	(dset->n == 0 ? 0 : sample_size(dset));
    case R_NVARS:
	return (dset == NULL)? NADBL : dset->v;
    case R_PD:
	return (no_data(dset))? NADBL : dset->pd;
    case R_T1:
	return (no_data(dset))? NADBL : dset->t1 + 1;
    case R_T2:
	return (no_data(dset))? NADBL : dset->t2 + 1;
    case R_TMAX:
	if (no_data(dset)) {
	    return NADBL;
	} else {
	    int tmax;

	    sample_range_get_extrema(dset, NULL, &tmax);
	    return tmax + 1;
	}
    case R_DATATYPE:
	return (no_data(dset))? NADBL : dataset_get_structure(dset);
    case R_TEST_PVAL:
	return get_last_pvalue();
    case R_TEST_STAT:
	return get_last_test_statistic();
    case R_TEST_LNL:
	return get_last_lnl();
    case R_TEST_BRK:
	return get_last_break();
    case R_STOPWATCH:
	return gretl_stopwatch();
    case R_WINDOWS:
#ifdef WIN32
	return 1;
#else
	return 0;
#endif
    case R_VERSION:
	return gretl_version_number(GRETL_VERSION);
    case R_ERRNO:
	return get_gretl_errno();
    case R_SEED:
	return gretl_rand_get_seed();
    case R_HUGE:
	return libset_get_double(CONV_HUGE);
    default:
	return NADBL;
    }
}

static int date_series_ok (const DATASET *dset)
{
    if (calendar_data(dset)) {
	return 1;
    } else if (quarterly_or_monthly(dset)) {
	return 1;
    } else if (annual_data(dset) || decennial_data(dset)) {
	return 1;
    } else if (dataset_has_panel_time(dset)) {
	return 1;
    } else {
	return 0;
    }
}

static int dvar_get_series (double *x, int i, const DATASET *dset)
{
    int t, YMD = calendar_data(dset);
    int err = 0;

    if (dset == NULL || dset->n == 0) {
	return E_NODATA;
    }

    if (i == R_OBSMIN && dset->pd < 2) {
	return E_PDWRONG;
    }

    if (i == R_OBSMIC && !YMD) {
	return E_PDWRONG;
    }

    if (i == R_PUNIT && !dataset_is_panel(dset)) {
	return E_PDWRONG;
    }

    if (i == R_DATES && !date_series_ok(dset)) {
	return E_PDWRONG;
    }

    if (i == R_OBSMAJ) {
	if (dset->pd == 1 && !dataset_is_time_series(dset)) {
	    i = R_INDEX;
	} else if (dataset_is_panel(dset)) {
	    i = R_PUNIT;
	}
    }

    if (i == R_INDEX) {
	for (t=0; t<dset->n; t++) {
	    x[t] = t + 1;
	}
    } else if (YMD && i != R_INDEX && i != R_DATES) {
	/* Watch out: we're handling most calendar-data cases
	   here, so we have to explicitly exclude cases that
	   require different treatment.
	*/
	char obs[12];
	int y, m, d;

	for (t=0; t<dset->n && !err; t++) {
	    ntodate(obs, t, dset);
	    if (sscanf(obs, YMD_READ_FMT, &y, &m, &d) != 3) {
		err = E_DATA;
	    } else if (i == R_OBSMAJ) {
		x[t] = y;
	    } else if (i == R_OBSMIN) {
		x[t] = m;
	    } else {
		x[t] = d;
	    }
	}
    } else if (i == R_PUNIT) {
	for (t=0; t<dset->n; t++) {
	    x[t] = t / dset->pd + 1;
	}
    } else if (i == R_OBSMAJ) {
	int maj;

	for (t=0; t<dset->n; t++) {
	    date_maj_min(t, dset, &maj, NULL);
	    x[t] = maj;
	}
    } else if (i == R_OBSMIN) {
	int min;

	for (t=0; t<dset->n; t++) {
	    date_maj_min(t, dset, NULL, &min);
	    x[t] = min;
	}
    } else if (i == R_DATES) {
	err = fill_dataset_dates_series(dset, x);
    } else {
	err = E_DATA;
    }

    return err;
}

static gretl_matrix *dvar_get_matrix (int i, int *err)
{
    gretl_matrix *m = NULL;

    switch (i) {
    case R_TEST_STAT:
	m = get_last_test_matrix(err);
	break;
    case R_TEST_PVAL:
	m = get_last_pvals_matrix(err);
	break;
    default:
	*err = E_DATA;
	break;
    }

    return m;
}

static gretl_matrix *dvar_get_submatrix (NODE *targ, int idx,
					 NODE *t, parser *p)
{
    NODE *r = eval(t->v.b2.r, p);
    gretl_matrix *M, *S = NULL;

    if (r == NULL || r->t != MSPEC) {
	if (!p->err) {
	    p->err = E_TYPES;
	}
	return NULL;
    }

    M = dvar_get_matrix(idx, &p->err);

    if (M != NULL) {
	S = matrix_get_submatrix(M, r->v.mspec, 0, &p->err);
	gretl_matrix_free(M);
    }

    return S;
}

static NODE *dollar_var_node (NODE *t, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	int idx, mslice = (t->t == DMSL);

	if (mslice) {
	    /* "slice" spec on right subnode, index on left */
	    idx = t->v.b2.l->v.idnum;
	} else {
	    idx = t->v.idnum;
	}

	if (mslice) {
	    ret = aux_matrix_node(p);
	    if (ret != NULL) {
		ret->v.m = dvar_get_submatrix(ret, idx, t, p);
	    }
	} else if (dvar_scalar(idx)) {
	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		ret->v.xval = dvar_get_scalar(idx, p->dset);
	    }
	} else if (dvar_series(idx)) {
	    ret = aux_series_node(p);
	    if (ret != NULL) {
		p->err = dvar_get_series(ret->v.xvec, idx, p->dset);
	    }
	} else if (dvar_variant(idx)) {
	    GretlType type = get_last_test_type();

	    if (type == GRETL_TYPE_MATRIX) {
		ret = aux_matrix_node(p);
		if (ret != NULL) {
		    ret->v.m = dvar_get_matrix(idx, &p->err);
		}
	    } else {
		/* scalar or none */
		ret = aux_scalar_node(p);
		if (ret != NULL) {
		    ret->v.xval = dvar_get_scalar(idx, p->dset);
		}
	    }
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

/* The incoming node here is binary: matrix ID on the left, and
   specification of a matrix sub-slice on the right.
*/

static gretl_matrix *
object_var_get_submatrix (NODE *targ, const char *oname, int idx,
			  NODE *t, parser *p, int needs_data)
{
    NODE *r = eval(t->v.b2.r, p);
    gretl_matrix *M, *S = NULL;

    if (r == NULL || r->t != MSPEC) {
	if (!p->err) {
	    p->err = E_TYPES;
	}
	return NULL;
    }

    if (needs_data) {
	M = saved_object_build_matrix(oname, idx, p->dset, &p->err);
    } else {
	M = saved_object_get_matrix(oname, idx, &p->err);
    }

    if (M != NULL) {
	S = matrix_get_submatrix(M, r->v.mspec, 0, &p->err);
	gretl_matrix_free(M);
    }

    return S;
}

static GretlType object_var_type (int idx, const char *oname,
				  int *needs_data)
{
    GretlType vtype = GRETL_TYPE_NONE;

    *needs_data = 0;

    if (model_data_scalar(idx)) {
	vtype = GRETL_TYPE_DOUBLE;
    } else if (model_data_series(idx)) {
	vtype = GRETL_TYPE_SERIES;
    } else if (model_data_matrix(idx)) {
	vtype = GRETL_TYPE_MATRIX;
    } else if (model_data_matrix_builder(idx)) {
	vtype = GRETL_TYPE_MATRIX;
	*needs_data = 1;
    } else if (model_data_list(idx)) {
	vtype = GRETL_TYPE_LIST;
    } else if (model_data_string(idx)) {
	vtype = GRETL_TYPE_STRING;
    } else if (idx == B_MODEL) {
	vtype = GRETL_TYPE_BUNDLE;
    }

    if (idx == M_UHAT || idx == M_YHAT || idx == M_SIGMA) {
	/* could be a matrix */
	vtype = saved_object_get_data_type(oname, idx);
    }

    return vtype;
}

/* For example, $coeff(sqft); or model1.$coeff(const) */

static NODE *dollar_str_node (NODE *t, MODEL *pmod, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	NODE *l = t->v.b2.l;
	NODE *r = t->v.b2.r;

	ret->v.xval = gretl_model_get_data_element(pmod, l->v.idnum, r->v.str,
						   p->dset, &p->err);

	if (na(ret->v.xval)) {
	    const char *s = get_string_by_name(r->v.str);

	    if (s != NULL) {
		p->err = 0;
		ret->v.xval = gretl_model_get_data_element(pmod, l->v.idnum, s,
							   p->dset, &p->err);
	    }
	}

	if (na(ret->v.xval)) {
	    p->err = E_INVARG;
	    pprintf(p->prn, _("'%s': invalid argument for %s()\n"),
		    r->v.str, mvarname(l->v.idnum));
	}
    }

    return ret;
}

/* Retrieve a data item from an object (typically, a model). We're not
   sure in advance here of the type of the data item (scalar, matrix,
   etc.): we look that up with object_var_type().

   This function handles two cases: the input @t may be a binary node
   with the name of an object on its left-hand subnode (as in
   "mymodel.$vcv"), or the object may be implicit (as when accessing
   data from the last model). In the latter case @t is itself the data
   item specification, while in the former the data item spec will be
   found on the right-hand subnode of @t.

   And there's another thing: in the case where the data item to be
   retrieved is a matrix, we handle the possibility that the user
   actually wants a sub-slice of that matrix. Handling that in
   this function may be the wrong thing to do -- perhaps it should be
   subsumed under the more general handling of <matrix>[<subspec>].
*/

static NODE *object_var_node (NODE *t, parser *p)
{
    NODE *ret = NULL;

#if EDEBUG
    fprintf(stderr, "object_var_node: t->t = %d\n", t->t);
#endif

    if (starting(p)) {
	NODE *r = NULL; /* the data spec node */
	const char *oname = NULL;
	GretlType vtype;
	int idx, mslice;
	int needs_data = 0;

	if (t->t == OVAR) {
	    /* objectname.<stuff> */
	    oname = t->v.b2.l->v.str;
	    r = t->v.b2.r;
	} else {
	    /* plain <stuff> */
	    r = t;
	}

	if (oname != NULL && gretl_get_object_by_name(oname) == NULL) {
	    gretl_errmsg_sprintf(_("%s: no such object\n"), oname);
	    p->err = E_UNKVAR;
	    return NULL;
	}

	if (oname != NULL && r->t == DMSTR) {
	    MODEL *pmod = get_model_by_name(oname);

	    if (pmod == NULL) {
		p->err = E_INVARG;
		return NULL;
	    } else {
		return dollar_str_node(r, pmod, p);
	    }
	}

	mslice = (r->t == DMSL);

	/* find the index which identifies the data item
	   the user wants */

	if (mslice) {
	    /* slice spec is on right subnode, index on left */
	    idx = r->v.b2.l->v.idnum;
	} else {
	    idx = r->v.idnum;
	}

	if (mslice && dvar_variant(idx)) {
	    /* we're in the wrong place */
	    return dollar_var_node(r, p);
	}

	vtype = object_var_type(idx, oname, &needs_data);

#if EDEBUG
	fprintf(stderr, "object_var_node: t->t = %d (%s), r->t = %d (%s)\n",
		t->t, getsymb(t->t), r->t, getsymb(r->t));
	fprintf(stderr, "idx = %d, vtype = %d, mslice = %d\n",
		idx, vtype, mslice);
#endif

	if (vtype == GRETL_TYPE_DOUBLE) {
	    ret = aux_scalar_node(p);
	} else if (vtype == GRETL_TYPE_SERIES) {
	    ret = aux_series_node(p);
	} else if (vtype == GRETL_TYPE_LIST) {
	    ret = aux_list_node(p);
	} else if (vtype == GRETL_TYPE_STRING) {
	    ret = aux_string_node(p);
	} else if (vtype == GRETL_TYPE_BUNDLE) {
	    ret = aux_bundle_node(p);
	} else {
	    ret = aux_matrix_node(p);
	}

	if (ret == NULL) {
	    return ret;
	} else if (vtype == GRETL_TYPE_DOUBLE) {
	    ret->v.xval = saved_object_get_scalar(oname, idx, p->dset,
						  &p->err);
	} else if (vtype == GRETL_TYPE_SERIES) {
	    p->err = saved_object_get_series(ret->v.xvec, oname, idx, p->dset);
	} else if (vtype == GRETL_TYPE_LIST) {
	    ret->v.ivec = saved_object_get_list(oname, idx, &p->err);
	} else if (vtype == GRETL_TYPE_STRING) {
	    ret->v.str = saved_object_get_string(oname, idx, p->dset,
						 &p->err);
	} else if (vtype == GRETL_TYPE_BUNDLE) {
	    ret->v.b = bundle_from_model(NULL, p->dset, &p->err);
	} else if (mslice) {
	    /* the right-hand subnode needs more work */
	    ret->v.m = object_var_get_submatrix(ret, oname, idx, r, p, needs_data);
	} else if (vtype == GRETL_TYPE_MATRIX) {
	    if (needs_data) {
		ret->v.m = saved_object_build_matrix(oname, idx,
						     p->dset, &p->err);
	    } else {
		ret->v.m = saved_object_get_matrix(oname, idx, &p->err);
	    }
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *wildlist_node (NODE *n, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	int *list = varname_match_list(p->dset, n->v.str,
				       &p->err);

	ret->v.ivec = list;
    }

    return ret;
}

static NODE *ellipsis_list_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	int *list = ellipsis_list(p->dset, l->vnum, r->vnum,
				  &p->err);

	ret->v.ivec = list;
    }

    return ret;
}

static NODE *list_join_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	int *L1 = node_get_list(l, p);
	int *L2 = node_get_list(r, p);

	if (!p->err) {
	    ret->v.ivec = gretl_lists_join_with_separator(L1, L2);
	    if (ret->v.ivec == NULL) {
		p->err = E_ALLOC;
	    }
	}

	free(L1);
	free(L2);
    }

    return ret;
}

static NODE *two_scalars_func (NODE *l, NODE *r, int t, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	double xl = node_get_scalar(l, p);
	double xr = node_get_scalar(r, p);

	if (!na(xl) && !na(xr)) {
	    if (t == F_XMIN) {
		ret->v.xval = (xl < xr)? xl : xr;
	    } else if (t == F_XMAX) {
		ret->v.xval = (xl > xr)? xl : xr;
	    } else if (t == F_RANDINT) {
		int k;

		p->err = gretl_rand_int_minmax(&k, 1, xl, xr);
		if (!p->err) {
		    ret->v.xval = k;
		}
	    }
	}
    }

    return ret;
}

static NODE *kpss_crit_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	int T = node_get_int(l, p);
	int trend = node_get_int(r, p);

	if (!p->err) {
	    ret->v.m = kpss_critvals(T, trend, &p->err);
	}
    }

    return ret;
}

static NODE *scalar_postfix_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	double x = n->v.xval;

	ret->v.xval = x;
	if (n->t == NUM_P) {
	    p->err = node_replace_scalar(n, x + 1.0);
	} else {
	    p->err = node_replace_scalar(n, x - 1.0);
	}
    }

    return ret;
}

static int series_calc_nodes (NODE *l, NODE *r)
{
    int ret = 0;

    if (l->t == SERIES) {
	ret = (r->t == SERIES || r->t == NUM || scalar_matrix_node(r));
    } else if (r->t == SERIES) {
	ret = scalar_node(l);
    }

    return ret;
}

static int cast_series_to_list (parser *p, NODE *n, short f)
{
    if (p->targ == LIST && useries_node(n)) {
	return (f == F_LOG || f == F_DIFF ||
		f == F_LDIFF || f == F_SDIFF ||
		f == F_ODEV);
    } else {
	return 0;
    }
}

static NODE *lhs_terminal_node (NODE *t, NODE *l, NODE *r,
				parser *p)
{
    /* Pass through eval'd @l and @r subnodes, but don't eval
       the parent @t itself */
    NODE *ret = aux_b2_node(p);

    ret->t = t->t;   /* transcribe type */
    ret->v.b2.l = l; /* evaluated left-hand */
    ret->v.b2.r = r; /* evaluated right-hand */

    /* prevent double-freeing of children @l and @r */
    ret->flags |= LHT_NODE;

    return ret;
}

/* reattach_series: on successive executions of a given
   compiled "genr", the "xvec" pointer recorded on a
   SERIES node will have become invalid if:

   (1) the series has been renamed,
   (2) the dataset has been differently sub-sampled, or
   (3) the series has been deleted (should be impossible).

   In case (2) the ID number of the series should still
   be valid, so it ought to be sufficient to reconnect the
   xvec pointer as in the second branch below. In case (1),
   however, it's necessary to look up the ID number by name
   again (as things stand).

   Note that n->v.xvec will be NULL when this function is
   reached only in case a genr is attached to a loop that
   is saved across function calls -- then the pointer is
   reset to NULL on each call to the function (but not
   on each iteration of the loop itself).

   The "1" below, in the first condition, is a safety-first
   measure that should be removable if the renaming case
   is handled effectively by calling get_loop_renaming():
   this returns non-zero if there's any renaming of series
   going on within a loop on the current execution stack.
*/

static void reattach_series (NODE *n, parser *p)
{
    if (/* 1 || */ n->v.xvec == NULL || get_loop_renaming()) {
	/* do a full reset */
	n->vnum = current_series_index(p->dset, n->vname);
	if (n->vnum < 0) {
	    gretl_errmsg_sprintf("'%s': not a series", n->vname);
	    p->err = E_DATA;
	} else {
	    n->v.xvec = p->dset->Z[n->vnum];
	}
    } else if (n->vnum >= 0 && n->vnum < p->dset->v) {
	/* handles the case of the dataset moving due to
	   sub-sampling, provided that no series are deleted
	   (the name -> ID mapping remains unchanged)
	*/
 	n->v.xvec = p->dset->Z[n->vnum];
    }
}

static void reattach_data_error (NODE *n, parser *p)
{
    char msg[256];

    sprintf(msg, "'%s': expected %s", n->vname, getsymb(n->t));

    if (n->uv == NULL) {
	strcat(msg, " but name look-up failed");
	p->err = E_DATA;
    } else {
	const char *s = gretl_type_get_name(n->uv->type);
	gchar *add;

	if (n->uv->ptr == NULL) {
	    add = g_strdup_printf(", found %s with NULL data", s);
	    p->err = E_DATA;
	} else {
	    add = g_strdup_printf(" but found %s", s);
	    p->err = E_TYPES;
	}
	strcat(msg, add);
	gretl_errmsg_set(msg);
	g_free(add);
    }
}

static void node_reattach_data (NODE *n, parser *p)
{
    if (n->t == SERIES) {
	reattach_series(n, p);
    } else {
	GretlType type = 0;
	void *data = NULL;

	if (n->uv == NULL) {
	    n->uv = get_user_var_by_name(n->vname);
	}

	if (n->uv != NULL) {
	    data = n->uv->ptr;
	    type = n->uv->type;
	}

	if (data == NULL) {
	    reattach_data_error(n, p);
	} else if (uscalar_node(n)) {
	    if (type == GRETL_TYPE_DOUBLE) {
		n->v.xval = *(double *) data;
#if ONE_BY_ONE_CAST
	    } else if (type == GRETL_TYPE_MATRIX) {
		/* allow type-mutation */
		n->t = MAT;
		n->v.m = (gretl_matrix *) data;
#endif
	    } else {
		reattach_data_error(n, p);
	    }
	} else if (n->t == MAT && type == GRETL_TYPE_MATRIX) {
	    n->v.m = data;
	} else if (n->t == LIST && type == GRETL_TYPE_LIST) {
	    n->v.ivec = data;
	} else if (n->t == BUNDLE && type == GRETL_TYPE_BUNDLE) {
	    n->v.b = data;
	} else if (n->t == STR && type == GRETL_TYPE_STRING) {
	    n->v.str = data;
	} else if (n->t == ARRAY && type == GRETL_TYPE_ARRAY) {
	    n->v.a = data;
	} else {
	    reattach_data_error(n, p);
	}
    }
}

static void node_type_error (int ntype, int argnum, int goodt,
			     NODE *bad, parser *p)
{
    const char *nstr;

    if (ntype == 0) {
	p->err = E_TYPES;
	return;
    }

    if (ntype == LAG) {
	nstr = (goodt == NUM)? "lag order" : "lag variable";
    } else {
	nstr = getsymb(ntype);
    }

    if (goodt == 0 && bad != NULL) {
	if (bad->t == EMPTY) {
	    pprintf(p->prn, _("%s: insufficient arguments"), nstr);
	} else if (argnum <= 0) {
	    pprintf(p->prn, _("%s: invalid argument type %s"),
		    nstr, typestr(bad->t));
	} else {
	    pprintf(p->prn, _("%s: argument %d: invalid type %s"),
		    nstr, argnum, typestr(bad->t));
	}
	pputc(p->prn, '\n');
	p->err = E_TYPES;
	return;
    }

    if (ntype < OP_MAX) {
	if (argnum <= 0) {
	    pprintf(p->prn, _("%s: operand should be %s"),
		    nstr, typestr(goodt));
	} else {
	    pprintf(p->prn, _("%s: operand %d should be %s"),
		    nstr, argnum, typestr(goodt));
	}
    } else {
	if (argnum <= 0) {
	    pprintf(p->prn, _("%s: argument should be %s"),
		    nstr, typestr(goodt));
	} else {
	    pprintf(p->prn, _("%s: argument %d should be %s"),
		    nstr, argnum, typestr(goodt));
	}
    }

    if (bad != NULL) {
	pprintf(p->prn, _(", is %s"), typestr(bad->t));
    }
    pputc(p->prn, '\n');

    if (!strcmp(nstr, "&")) {
	pputs(p->prn, "(for logical AND, please use \"&&\")\n");
    } else if (!strcmp(nstr, "|")) {
	pputs(p->prn, "(for logical OR, please use \"||\")\n");
    }

    p->err = E_TYPES;
}

static int node_is_true (NODE *n, parser *p)
{
    double x = node_get_scalar(n, p);

    return !na(x) && x != 0.0;
}

static int node_is_false (NODE *n, parser *p)
{
    return (node_get_scalar(n, p) == 0.0);
}

#define eval_left(t)   (evalb1(t) || evalb2(t) || evalb3(t))
#define eval_middle(t) (evalb3(t))
#define eval_right(t)  (evalb2(t) || evalb3(t))

/* Fetch the appropriate child node for evaluation based
   on the structure of the parent, @t, and the index @i,
   which indicates position 0, 1 or 2.
*/

static NODE *input_node (NODE *t, int i)
{
    if (evalb1(t->t)) {
	return t->v.b1.b;
    } else if (evalb2(t->t)) {
	return (i==0)? t->v.b2.l : t->v.b2.r;
    } else {
	/* 3-place node */
	return (i==0)? t->v.b3.l : (i==1)? t->v.b3.m : t->v.b3.r;
    }
}

/* core function: evaluate the parsed syntax tree */

static NODE *eval (NODE *t, parser *p)
{
    NODE *l = NULL, *m = NULL, *r = NULL;
    NODE *ret = NULL;

    if (t == NULL) {
	return NULL;
    }

#if EDEBUG
    if (t->vname != NULL) {
	fprintf(stderr, "eval: incoming node %p ('%s', vname=%s)\n",
		(void *) t, getsymb(t->t), t->vname);
    } else {
	fprintf(stderr, "eval: incoming node %p ('%s')\n",
		(void *) t, getsymb(t->t));
    }
#endif

    if (!p->err && eval_left(t->t)) {
	l = eval(input_node(t, 0), p);
	if (l == NULL && !p->err) {
	    p->err = 1;
	}
    }

    if (!p->err && eval_middle(t->t)) {
	if (m_return(t->t)) {
	    m = t->v.b3.m;
	} else {
	    m = eval(t->v.b3.m, p);
	    if (m == NULL && !p->err) {
		p->err = 1;
	    }
	}
    }

    if (!p->err && eval_right(t->t)) {
	if (r_return(t->t)) {
	    r = input_node(t, 2);
	} else {
	    if (t->t == B_AND || t->t == B_OR) {
		/* logical operators: avoid redundant evaluation */
		if (l != NULL && l->t == NUM) {
		    if ((t->t == B_AND && node_is_false(l, p)) ||
			(t->t == B_OR && node_is_true(l, p))) {
			/* no need to evaluate @r */
			ret = l;
			if (is_aux_node(ret)) {
			    /* mark return node as proxy */
			    ret->flags |= PRX_NODE;
			}
			goto finish;
		    }
		}
	    }
	    if (r == NULL && !p->err) {
		r = eval(input_node(t, 2), p);
	    }
	}
    }

    if (p->err) {
	goto bailout;
    }

    /* establish convenience pointer */
    p->aux = t->aux;

    switch (t->t) {
    case MSPEC:
    case EMPTY:
	ret = t;
	break;
    case DBUNDLE:
	ret = dollar_bundle_node(t, p);
	break;
    case NUM:
    case SERIES:
    case MAT:
    case STR:
    case LIST:
    case BUNDLE:
    case ARRAY:
	if (compiled(p) && starting(p) && uvar_node(t)) {
	    node_reattach_data(t, p);
	}
	ret = t;
	break;
    case UNDEF:
	ret = maybe_rescue_undef_node(t, p);
	break;
    case U_ADDR:
	if (!uvar_node(t->v.b1.b)) {
	    p->err = E_DATA;
	} else if (compiled(p) && starting(p)) {
	    node_reattach_data(t->v.b1.b, p);
	}
	ret = t;
	break;
    case WLIST:
	ret = wildlist_node(t, p);
	break;
    case DUM:
	if (t->v.idnum == DUM_DATASET) {
	    ret = dataset_list_node(p);
	} else if (t->v.idnum == DUM_TREND) {
	    ret = trend_node(p);
	} else {
	    /* otherwise treat as terminal */
	    ret = t;
	}
	break;
    case NUM_P:
    case NUM_M:
	if (compiled(p) && starting(p)) {
	    node_reattach_data(t, p);
	}
	ret = scalar_postfix_node(t, p);
	break;
    case FARGS:
	/* will be evaluated in context */
	ret = t;
	break;
    case B_ADD:
    case B_SUB:
    case B_MUL:
    case B_DIV:
    case B_MOD:
    case B_POW:
    case B_AND:
    case B_OR:
    case B_EQ:
    case B_NEQ:
    case B_GT:
    case B_LT:
    case B_GTE:
    case B_LTE:
	/* arithmetic and logical binary operators: be as
	   flexible as possible with regard to argument types
	*/
	if (t->t == B_ADD && l->t == STR && r->t == NUM) {
	    ret = string_offset(l, r, p);
	} else if ((t->t == B_EQ || t->t == B_NEQ) && l->t == STR && r->t == STR) {
	    ret = compare_strings(l, r, t->t, p);
	} else if (l->t == NUM && r->t == NUM) {
	    ret = scalar_calc(l, r, t->t, p);
	} else if (l->t == BUNDLE && r->t == BUNDLE) {
	    ret = bundle_op(l, r, t->t, p);
	} else if (l->t == ARRAY && r->t == ARRAY) {
	    ret = array_op(l, r, t->t, p);
	} else if (stringvec_node(l) && stringvec_node(r)) {
	    ret = stringvec_calc(l, r, t->t, p);
	} else if (series_calc_nodes(l, r)) {
	    ret = series_calc(l, r, t->t, p);
	} else if (l->t == MAT && r->t == MAT) {
	    if (bool_comp(t->t)) {
		ret = matrix_bool(l, r, t->t, p);
	    } else {
		ret = matrix_matrix_calc(l, r, t->t, p);
	    }
	} else if ((l->t == MAT && r->t == NUM) ||
		   (l->t == NUM && r->t == MAT)) {
	    ret = matrix_scalar_calc(l, r, t->t, p);
	} else if ((l->t == MAT && r->t == SERIES) ||
		   (l->t == SERIES && r->t == MAT)) {
	    ret = matrix_series_calc(l, r, t->t, p);
	} else if (t->t >= B_EQ && t->t <= B_NEQ &&
		   ((l->t == SERIES && r->t == STR) ||
		    (l->t == STR && r->t == SERIES))) {
	    ret = series_string_calc(l, r, t->t, p);
	} else if ((t->t == B_AND || t->t == B_OR || t->t == B_SUB) &&
		   ok_list_node(l) && ok_list_node(r)) {
	    ret = list_list_op(l, r, t->t, p);
	} else if (t->t == B_POW && ok_list_node(l) && ok_list_node(r)) {
	    ret = list_list_op(l, r, t->t, p);
	} else if ((t->t == B_EQ || t->t == B_NEQ) &&
		   ((l->t == NUM && r->t == STR) ||
		    (l->t == STR && r->t == NUM))) {
	    ret = num_string_comp(l, r, t->t, p);
	} else if (bool_comp(t->t)) {
	    if (ok_list_node(l) && (r->t == NUM || r->t == SERIES)) {
		ret = list_bool_comp(l, r, t->t, 0, p);
	    } else if (ok_list_node(r) && (l->t == NUM || l->t == SERIES)) {
		ret = list_bool_comp(r, l, t->t, 1, p);
	    } else if (ok_list_node(l) && ok_list_node(r)) {
		ret = list_list_comp(r, l, t->t, p);
	    } else {
		p->err = E_TYPES;
	    }
	} else if ((t->t == B_ADD || t->t == B_SUB) &&
		   l->t == SERIES && ok_list_node(r)) {
	    ret = series_list_calc(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case B_TRMUL:
	/* matrix on left, otherwise be flexible */
	if (l->t == MAT && r->t == MAT) {
	    ret = matrix_matrix_calc(l, r, t->t, p);
	} else if (l->t == MAT && r->t == SERIES) {
	    ret = matrix_series_calc(l, r, t->t, p);
	} else if ((l->t == MAT && r->t == NUM) ||
		   (l->t == NUM && r->t == MAT)) {
	    ret = matrix_scalar_calc(l, r, t->t, p);
	} else if (l->t == MAT && r->t == EMPTY) {
	    ret = matrix_transpose_node(l, p);
	} else if (l->t == NUM && r->t == EMPTY) {
	    ret = l;
	} else {
	    p->err = E_TYPES;
	}
	break;
    case B_DOTMULT:
    case B_DOTDIV:
    case B_DOTPOW:
    case B_DOTADD:
    case B_DOTSUB:
    case B_DOTEQ:
    case B_DOTGT:
    case B_DOTLT:
    case B_DOTGTE:
    case B_DOTLTE:
    case B_DOTNEQ:
	/* matrix-matrix or matrix-scalar binary operators:
	   in addition we permit scalar-scalar to allow for
	   the possibility that results that could be taken
	   to be 1 x 1 matrix results have been registered
	   internally as scalars.
	*/
	if ((l->t == MAT && r->t == MAT) ||
	    (l->t == MAT && r->t == NUM) ||
	    (l->t == NUM && r->t == MAT) ||
	    (l->t == NUM && r->t == NUM)) {
	    ret = matrix_matrix_calc(l, r, t->t, p);
	} else if ((l->t == MAT && r->t == SERIES) ||
		   (l->t == SERIES && r->t == MAT)) {
	    ret = matrix_series_calc(l, r, t->t, p);
	} else {
	    node_type_error(t->t, (l->t == MAT)? 2 : 1,
			    MAT, (l->t == MAT)? r : l, p);
	}
	break;
    case B_HCAT:
    case B_VCAT:
    case F_QFORM:
    case F_HDPROD:
    case F_CMULT:
    case F_CDIV:
    case F_MRSEL:
    case F_MCSEL:
    case F_DSUM:
    case B_LDIV:
    case B_KRON:
	/* matrix-only binary operators (but promote scalars) */
	if ((l->t == MAT || l->t == NUM) &&
	    (r->t == MAT || r->t == NUM)) {
	    ret = matrix_matrix_calc(l, r, t->t, p);
	} else if (t->t == B_HCAT && l->t == STR && r->t == STR) {
	    /* exception: string concatenation */
	    ret = two_string_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, (l->t == MAT)? 2 : 1,
			    MAT, (l->t == MAT)? r : l, p);
	}
	break;
    case B_ELLIP:
	/* list-making ellipsis */
	if (useries_node(l) && useries_node(r)) {
	    ret = ellipsis_list_node(l, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case B_JOIN:
	/* list join with separator */
	if (ok_list_node(l) && ok_list_node(r)) {
	    ret = list_join_node(l, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_MSORTBY:
    case HF_CXTRACT:
	/* matrix on left, scalar on right */
	if (l->t == MAT && scalar_node(r)) {
	    ret = matrix_scalar_func(l, r, t->t, p);
	} else if (l->t == MAT) {
	    node_type_error(t->t, 2, NUM, r, p);
	} else {
	    node_type_error(t->t, 1, MAT, l, p);
	}
	break;
    case F_LLAG:
	if ((scalar_node(l) || l->t == MAT) && m->t != MAT && ok_list_node(m)) {
	    ret = list_make_lags(l, m, r, p);
	} else if ((scalar_node(l) || l->t == MAT) && m->t == MAT) {
	    ret = matrix_make_lags(l, m, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_HFLAG:
	if (scalar_node(l) && scalar_node(m) && ok_list_node(r)) {
	    ret = hf_list_make_lags(l, m, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_HFLIST:
	if (l->t == MAT && scalar_node(m) && r->t == STR) {
	    ret = hf_list_node(l, m, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_HFDIFF:
    case F_HFLDIFF:
	if (ok_list_node(l) && empty_or_num(r)) {
	    ret = apply_list_func(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case U_NEG:
    case U_POS:
    case U_NOT:
    case F_ABS:
    case F_TOINT:
    case F_CEIL:
    case F_FLOOR:
    case F_ROUND:
    case F_SIN:
    case F_COS:
    case F_TAN:
    case F_ASIN:
    case F_ACOS:
    case F_ATAN:
    case F_SINH:
    case F_COSH:
    case F_TANH:
    case F_ASINH:
    case F_ACOSH:
    case F_ATANH:
    case F_LOG:
    case F_LOG10:
    case F_LOG2:
    case F_EXP:
    case F_SQRT:
    case F_CNORM:
    case F_DNORM:
    case F_QNORM:
    case F_LOGISTIC:
    case F_GAMMA:
    case F_LNGAMMA:
    case F_DIGAMMA:
    case F_INVMILLS:
    case F_EASTER:
	/* functions taking one argument, any type */
	if (l->t == NUM) {
	    ret = apply_scalar_func(l, t->t, p);
	} else if (l->t == SERIES) {
	    if (cast_series_to_list(p, l, t->t)) {
		ret = apply_list_func(l, NULL, t->t, p);
	    } else {
		ret = apply_series_func(l, t->t, p);
	    }
	} else if (l->t == MAT) {
	    ret = apply_matrix_func(l, t->t, p);
	} else if (ok_list_node(l) && t->t == F_LOG) {
	    ret = apply_list_func(l, NULL, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_DUMIFY:
    case F_CDUMIFY:
	/* series or list argument wanted */
	if (ok_list_node(l)) {
	    if (t->t == F_CDUMIFY) {
		ret = cdummify_func(l, p);
	    } else {
		ret = dummify_func(l, r, p);
	    }
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_GETINFO:
	/* named series (or ID) argument wanted */
	if (useries_node(l) || l->t == NUM) {
	    ret = get_series_info(l, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_SEASONALS:
	/* two optional args: int, bool */
	if (!null_or_empty(l) && l->t != NUM) {
	    node_type_error(t->t, 1, NUM, l, p);
	} else if (!null_or_empty(r) && r->t != NUM) {
	    node_type_error(t->t, 2, NUM, r, p);
	} else {
	    ret = seasonals_node(l, r, p);
	}
	break;
    case F_MISSZERO:
    case F_ZEROMISS:
	/* one series or scalar argument needed */
	if (l->t == SERIES || l->t == MAT) {
	    ret = apply_series_func(l, t->t, p);
	} else if (l->t == NUM) {
	    ret = apply_scalar_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, SERIES, l, p);
	}
	break;
    case F_MISSING:
    case F_DATAOK:
	/* series, scalar or list argument needed */
	if (l->t == MAT) {
	    if (t->t == F_DATAOK) {
		ret = matrix_to_matrix_func(l, NULL, t->t, p);
	    } else {
		node_type_error(t->t, 0, SERIES, l, p);
	    }
	} else if (l->t == SERIES) {
	    ret = apply_series_func(l, t->t, p);
	} else if (l->t == NUM) {
	    ret = apply_scalar_func(l, t->t, p);
	} else if (ok_list_node(l)) {
	    ret = list_ok_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, SERIES, l, p);
	}
	break;
    case F_ISNAN:
	/* scalar or matrix */
	if (scalar_node(l)) {
	    ret = scalar_isnan_node(l, p);
	} else if (l->t == MAT) {
	    ret = matrix_isnan_node(l, p);
	} else {
	    node_type_error(t->t, 0, NUM, l, p);
	}
	break;
    case F_SLEEP:
	if (scalar_node(l)) {
	    ret = sleep_node(l, p);
	} else {
	    node_type_error(t->t, 0, NUM, l, p);
	}
	break;
    case LAG:
	if (p->targ == LIST) {
	    ret = get_lag_list(l, r, p);
	    break;
	}
	/* Falls through. */
    case F_LJUNGBOX:
    case F_POLYFIT:
	/* series on left, scalar on right */
	if (l->t != SERIES) {
	    node_type_error(t->t, 1, SERIES, l, p);
	} else if (!scalar_node(r)) {
	    node_type_error(t->t, 2, NUM, r, p);
	} else if (t->t == LAG) {
	    ret = series_lag(l, r, p);
	} else if (t->t == F_LJUNGBOX) {
	    ret = series_ljung_box(l, r, p);
	} else if (t->t == F_POLYFIT) {
	    ret = series_polyfit(l, r, p);
	}
	break;
    case OBS:
	if (l->t != SERIES) {
	    node_type_error(t->t, 1, SERIES, l, p);
	} else if (!scalar_node(r) && r->t != STR) {
	    node_type_error(t->t, 2, NUM, r, p);
	} else if (t->flags & LHT_NODE) {
	    ret = lhs_terminal_node(t, l, r, p);
	} else {
	    ret = series_obs(l, r, p);
	}
	break;
    case MSL:
	/* user matrix plus subspec */
	if (t->flags & LHT_NODE) {
	    ret = lhs_terminal_node(t, l, r, p);
	} else {
	    ret = submatrix_node(l, r, p);
	}
	break;
    case OSL:
	/* object plus subspec */
	if (t->flags & LHT_NODE) {
	    ret = lhs_terminal_node(t, l, r, p);
	} else {
	    ret = subobject_node(l, r, p);
	}
	break;
    case MSLRAW:
	/* unevaluated matrix subspec */
	ret = mspec_node(l, r, p);
	break;
    case SUBSL:
    case B_RANGE:
	/* matrix sub-slice, x:y, or lag range, 'p to q' */
	ret = process_subslice(l, r, p);
	break;
    case BMEMB:
    case F_INBUNDLE:
	/* name of bundle plus string */
	if (l->t == BUNDLE && r->t == STR) {
	    if (t->t == BMEMB) {
		if (t->flags & LHT_NODE) {
		    ret = lhs_terminal_node(t, l, r, p);
		} else {
		    ret = get_bundle_member(l, r, p);
		}
	    } else {
		ret = test_bundle_key(l, r, p);
	    }
	} else if (l->t == BUNDLE) {
	    node_type_error(t->t, 1, STR, r, p);
	} else {
	    node_type_error(t->t, 0, BUNDLE, l, p);
	}
	break;
    case DBMEMB:
	/* name of $ bundle plus string */
	if ((l->t == DBUNDLE || l->t == BUNDLE) && r->t == STR) {
	    ret = get_bundle_member(l, r, p);
	} else if (l->t == DBUNDLE || l->t == BUNDLE) {
	    node_type_error(t->t, 1, STR, r, p);
	} else {
	    node_type_error(t->t, 0, DBUNDLE, l, p);
	}
	break;
    case F_CURL:
	ret = curl_bundle_node(l, p);
	break;
    case F_SVM:
	if (l->t == LIST && m->t == BUNDLE) {
	    ret = svm_predict_node(l, m, r, p);
	} else if (l->t == LIST) {
	    node_type_error(t->t, 1, BUNDLE, m, p);
	} else {
	    node_type_error(t->t, 0, LIST, l, p);
	}
	break;
    case F_TYPESTR:
	/* numerical type code to string */
	if (scalar_node(l)) {
	    ret = type_string_node(l, p);
	} else {
	    node_type_error(t->t, 0, NUM, l, p);
	}
	break;
    case F_LDIFF:
    case F_SDIFF:
    case F_ODEV:
	if (l->t == SERIES && cast_series_to_list(p, l, t->t)) {
	    ret = apply_list_func(l, NULL, t->t, p);
	} else if (l->t == SERIES || (t->t != F_ODEV && l->t == MAT)) {
	    ret = series_series_func(l, r, t->t, p);
	} else if (ok_list_node(l)) {
	    ret = apply_list_func(l, NULL, t->t, p);
	} else {
	    node_type_error(t->t, 0, SERIES, l, p);
	}
	break;
    case F_DROPCOLL:
	/* list argument is required on left, optional scalar
	   on the right */
	if (l->t == LIST) {
	    ret = apply_list_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, 0, LIST, l, p);
	}
	break;
    case F_HPFILT:
    case F_FRACDIFF:
    case F_FRACLAG:
    case F_BOXCOX:
    case F_PNOBS:
    case F_PMIN:
    case F_PMAX:
    case F_PSUM:
    case F_PMEAN:
    case F_PXSUM:
    case F_PXNOBS:
    case F_PSD:
    case F_DESEAS:
    case F_TRAMOLIN:
	/* series argument needed */
	if (l->t == SERIES || l->t == MAT) {
	    ret = series_series_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, 0, SERIES, l, p);
	}
	break;
    case F_FREQ:
	/* series -> matrix */
	if (l->t == SERIES || l->t == MAT) {
	    ret = series_matrix_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, SERIES, l, p);
	}
	break;
    case F_PSHRINK:
	if (l->t == SERIES) {
	    ret = do_panel_shrink(l, p);
	} else {
	    node_type_error(t->t, 0, SERIES, l, p);
	}
	break;
    case F_PEXPAND:
	if (l->t == MAT) {
	    ret = do_panel_expand(l, p);
	} else {
	    node_type_error(t->t, 0, MAT, l, p);
	}
	break;
    case F_CUM:
    case F_DIFF:
    case F_RESAMPLE:
    case F_RANKING:
	/* series or matrix argument */
	if (l->t == SERIES) {
	    if (cast_series_to_list(p, l, t->t)) {
		ret = apply_list_func(l, NULL, t->t, p);
	    } else {
		ret = series_series_func(l, r, t->t, p);
	    }
	} else if (l->t == MAT) {
	    ret = matrix_to_matrix_func(l, r, t->t, p);
	} else if (t->t == F_DIFF && ok_list_node(l)) {
	    ret = apply_list_func(l, NULL, t->t, p);
	} else if (t->t == F_RESAMPLE && ok_list_node(l) &&
		   null_or_empty(r)) {
	    ret = apply_list_func(l, NULL, t->t, p);
	} else {
	    node_type_error(t->t, 0, SERIES, l, p);
	}
	break;
    case F_SORT:
    case F_DSORT:
    case F_VALUES:
    case F_UNIQ:
    case F_PERGM:
    case F_IRR:
	/* series or vector argument needed */
	if (l->t == SERIES || l->t == MAT) {
	    if (t->t == F_PERGM) {
		ret = pergm_node(l, r, p);
	    } else if (t->t == F_VALUES || t->t == F_UNIQ) {
		ret = vector_values(l, t->t, p);
	    } else if (t->t == F_IRR) {
		ret = do_irr(l, p);
	    } else {
		ret = vector_sort(l, t->t, p);
	    }
	} else {
	    node_type_error(t->t, 0, SERIES, l, p);
	}
	break;
    case F_SUM:
    case F_SUMALL:
    case F_MEAN:
    case F_SD:
    case F_VCE:
    case F_SST:
    case F_SKEWNESS:
    case F_KURTOSIS:
    case F_MIN:
    case F_MAX:
    case F_MEDIAN:
    case F_GINI:
    case F_NOBS:
    case F_T1:
    case F_T2:
	/* functions taking series arg, returning scalar */
	if (l->t == SERIES || l->t == MAT) {
	    ret = series_scalar_func(l, t->t, p);
	} else if ((t->t == F_MEAN || t->t == F_SD ||
		    t->t == F_VCE || t->t == F_MIN ||
		    t->t == F_MAX || t->t == F_SUM ||
		    t->t == F_MEDIAN)
		   && ok_list_node(l)) {
	    /* list -> series also acceptable for these cases */
	    ret = list_to_series_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, SERIES, l, p);
	}
	break;
    case F_ECDF:
    case F_NORMTEST:
	/* series or vector (plus optional string arg for,
	   normtest); returns matrix */
	if (l->t != SERIES && l->t != MAT) {
	    p->err = E_TYPES;
	} else if (null_or_empty(r) || r->t == STR) {
	    ret = series_matrix_node(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_LRVAR:
    case F_NPV:
    case F_ISCONST:
    case F_ISDUMMY:
	/* takes series and scalar arg, returns scalar */
	if (l->t == SERIES || l->t == MAT) {
	    if (t->t == F_ISCONST || t->t == F_ISDUMMY ) {
		ret = isconst_or_dum_node(l, r, p, t->t);
	    } else if (scalar_node(r)) {
		if (t->t == F_QUANTILE && l->t == MAT) {
		    ret = matrix_quantiles_node(l, r, p);
		} else {
		    ret = series_scalar_scalar_func(l, r, t->t, p);
		}
	    } else {
		node_type_error(t->t, 2, NUM, r, p);
	    }
	} else {
	    node_type_error(t->t, 1, SERIES, l, p);
	}
	break;
    case F_QUANTILE:
	if (l->t == SERIES) {
	    if (scalar_node(r)) {
		ret = series_scalar_scalar_func(l, r, t->t, p);
	    } else {
		node_type_error(t->t, 2, NUM, r, p);
	    }
	} else if (l->t == MAT) {
	    if (r->t == MAT || scalar_node(r)) {
		ret = matrix_quantiles_node(l, r, p);
	    } else {
		node_type_error(t->t, 2, MAT, r, p);
	    }
	} else {
	    node_type_error(t->t, 1, (r->t == MAT)? MAT : SERIES,
			    l, p);
	}
	break;
    case F_RUNIFORM:
    case F_RNORMAL:
	/* functions taking zero or two scalars as args */
	if (scalar_node(l) && scalar_node(r)) {
	    ret = series_fill_func(l, r, t->t, p);
	} else if (l->t == EMPTY && r->t == EMPTY) {
	    ret = series_fill_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, (l->t == NUM)? 2 : 1,
			    NUM, (l->t == NUM)? r : l, p);
	}
	break;
    case F_COV:
    case F_COR:
    case F_FCSTATS:
    case F_NAALEN:
    case F_KMEIER:
	/* functions taking two series/vectors as args */
	if ((l->t == SERIES || l->t == MAT) &&
	    (r->t == SERIES || r->t == MAT)) {
	    ret = series_2_func(l, r, t->t, p);
	} else if ((l->t == SERIES || l->t == MAT) &&
		   null_or_empty(r) &&
		   (t->t == F_NAALEN || t->t == F_KMEIER)) {
	    ret = series_2_func(l, NULL, t->t, p);
	} else if ((l->t == SERIES || l->t == MAT) &&
		   r->t == LIST && t->t == F_FCSTATS) {
	    ret = series_2_func(l, r, t->t, p);
	} else {
	    p->err = E_INVARG;
	}
	break;
    case F_NPCORR:
	/* two series or vectors, plus optional control string */
	if ((l->t == SERIES || l->t == MAT) &&
	    (m->t == SERIES || m->t == MAT) &&
	    (null_or_empty(r) || r->t == STR)) {
	    ret = npcorr_node(l, m, r, p);
	} else {
	    p->err = E_INVARG;
	}
	break;
    case F_MXTAB:
	/* functions taking two series or matrices as args and returning
	   a matrix */
	if ((l->t == SERIES && r->t == SERIES) || (l->t == MAT && r->t == MAT)) {
	    ret = mxtab_func(l, r, p);
	} else {
	    node_type_error(t->t, (l->t == SERIES)? 2 : 1,
			    SERIES, (l->t == SERIES)? r : l, p);
	}
	break;
    case F_SORTBY:
	/* takes two series as args, returns series */
	if (l->t == SERIES && r->t == SERIES) {
	    ret = series_sort_by(l, r, p);
	} else {
	    node_type_error(t->t, (l->t == SERIES)? 2 : 1,
			    SERIES, (l->t == SERIES)? r : l, p);
	}
	break;
    case F_IMAT:
    case F_ZEROS:
    case F_ONES:
    case F_MUNIF:
    case F_MNORM:
	/* matrix-creation functions */
	if (scalar_node(l) && (r == NULL || scalar_node(r))) {
	    ret = matrix_fill_func(l, r, t->t, p);
	} else if (!scalar_node(l)) {
	    node_type_error(t->t, 1, NUM, l, p);
	} else {
	    node_type_error(t->t, 2, NUM, r, p);
	}
	break;
    case F_SUMC:
    case F_SUMR:
    case F_PRODC:
    case F_PRODR:
    case F_MEANC:
    case F_MEANR:
    case F_SDC:
    case F_MCOV:
    case F_MCORR:
    case F_CDEMEAN:
    case F_CHOL:
    case F_PSDROOT:
    case F_INV:
    case F_INVPD:
    case F_GINV:
    case F_DIAG:
    case F_TRANSP:
    case F_MREVERSE:
    case F_VEC:
    case F_VECH:
    case F_UNVECH:
    case F_UPPER:
    case F_LOWER:
    case F_NULLSPC:
    case F_MEXP:
    case F_MINC:
    case F_MAXC:
    case F_MINR:
    case F_MAXR:
    case F_IMINC:
    case F_IMAXC:
    case F_IMINR:
    case F_IMAXR:
    case F_FFT:
    case F_FFTI:
    case F_POLROOTS:
    case HF_CINV:
    case HF_CFFT:
    case HF_CTRAN:
	/* matrix -> matrix functions */
	if (l->t == MAT || l->t == NUM) {
	    ret = matrix_to_matrix_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, 1, MAT, l, p);
	}
	break;
    case F_ROWS:
    case F_COLS:
    case F_DET:
    case F_LDET:
    case F_TRACE:
    case F_NORM1:
    case F_INFNORM:
    case F_RCOND:
    case F_CNUMBER:
    case F_RANK:
	/* matrix -> scalar functions */
	if (l->t == MAT || l->t == NUM) {
	    ret = matrix_to_scalar_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, MAT, l, p);
	}
	break;
    case F_MREAD:
    case F_BREAD:
	if (l->t != STR) {
	    node_type_error(t->t, 1, STR, l, p);
	} else if (!empty_or_num(r)) {
	    node_type_error(t->t, 2, NUM, r, p);
	} else {
	    ret = read_object_func(l, r, t->t, p);
	}
	break;
    case F_QR:
    case F_EIGSYM:
    case F_EIGGEN:
    case HF_CEIGH:
	/* matrix -> matrix functions, with indirect return */
	if (l->t != MAT) {
	    node_type_error(t->t, 1, MAT, l, p);
	} else if (r->t != U_ADDR && r->t != EMPTY) {
	    node_type_error(t->t, 2, U_ADDR, r, p);
	} else {
	    ret = matrix_to_matrix2_func(l, r, t->t, p);
	}
	break;
    case HF_CMATRIX:
    case HF_CMMULT:
    case HF_CHPROD:
	if (l->t == MAT && r->t == MAT) {
	    ret = matrix_to_matrix_func(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_FDJAC:
    case F_NUMHESS:
	/* matrix, fncall, optional scalar */
	if (l->t == MAT && m->t == STR) {
	    ret = numeric_jacobian_or_hessian(l, m, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_MWRITE:
    case HF_CPRINTF:
	/* matrix, with string as second arg */
	if (t->t == F_MWRITE && l->t == MAT &&
	    m->t == STR && empty_or_num(r)) {
	    ret = matrix_file_write(l, m, r, p);
	} else if (t->t == HF_CPRINTF && empty_or_string(r)) {
	    ret = cmatrix_printf(l, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_BWRITE:
	/* bundle, with string as second arg */
	if (l->t == BUNDLE && m->t == STR && empty_or_num(r)) {
	    ret = bundle_file_write(l, m, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_BFGSMAX:
	/* matrix-pointer, plus one or two string args */
	if ((l->t == U_ADDR || l->t == MAT) && m->t == STR) {
	    ret = BFGS_maximize(l, m, r, p, t);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_BFGSCMAX:
	ret = BFGS_constrained_max(t, p);
	break;
    case F_SIMANN:
    case F_NMMAX:
    case F_GSSMAX:
	/* matrix-pointer, plus string and int args */
	if ((l->t == U_ADDR || l->t == MAT) && m->t == STR) {
	    ret = deriv_free_node(l, m, r, p, t);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_IMHOF:
	/* matrix, scalar as second arg */
	if (l->t == MAT && scalar_node(r)) {
	    ret = matrix_imhof(l, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_COLNAMES:
    case F_ROWNAMES:
	/* matrix, list, string or strings array as second arg */
	if (l->t == MAT && (ok_list_node(r) || r->t == STR || r->t == ARRAY)) {
	    ret = matrix_add_names(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_COLNAME:
	/* matrix, scalar as second arg */
	if (l->t == MAT && scalar_node(r)) {
	    ret = matrix_get_colname(l, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_XMIN:
    case F_XMAX:
    case F_RANDINT:
    case F_KPSSCRIT:
	/* two scalars */
	if (scalar_node(l) && scalar_node(r)) {
	    if (t->t == F_KPSSCRIT) {
		ret = kpss_crit_node(l, r, p);
	    } else {
		ret = two_scalars_func(l, r, t->t, p);
	    }
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_MSHAPE:
    case F_SVD:
    case F_TRIMR:
    case F_TOEPSOLV:
    case F_CORRGM:
    case F_SEQ:
    case F_REPLACE:
    case F_STRNCMP:
    case F_WEEKDAY:
    case F_DAYSPAN:
    case F_SMPLSPAN:
    case F_MONTHLEN:
    case F_EPOCHDAY:
    case F_KDENSITY:
    case F_SETNOTE:
    case F_BWFILT:
    case F_CHOWLIN:
    case F_VARSIMUL:
    case F_IRF:
    case F_STRSUB:
    case F_REGSUB:
    case F_MLAG:
    case F_EIGSOLVE:
    case F_NADARWAT:
    case F_PRINCOMP:
    case F_HALTON:
    case F_AGGRBY:
    case F_IWISHART:
    case F_SUBSTR:
    case F_MWEIGHTS:
    case F_MGRADIENT:
	/* built-in functions taking three args */
	if (t->t == F_REPLACE) {
	    ret = replace_value(l, m, r, p);
	} else if (t->t == F_STRSUB || t->t == F_REGSUB) {
	    ret = string_replace(l, m, r, t->t, p);
	} else if (t->t == F_EPOCHDAY) {
	    ret = eval_epochday(l, m, r, p);
	} else {
	    ret = eval_3args_func(l, m, r, t->t, p);
	}
	break;
    case F_PRINTF:
    case F_SPRINTF:
	if (l->t == STR && empty_or_string(r)) {
	    ret = eval_print_scan(NULL, l, r, t->t, p);
	} else {
	    node_type_error(t->t, 0, STR, NULL, p);
	}
	break;
    case F_SSCANF:
	if (l->t == STR && m->t == STR && r->t == STR) {
	    ret = eval_print_scan(l, m, r, t->t, p);
	} else {
	    node_type_error(t->t, 0, STR, NULL, p);
	}
	break;
    case F_STRSCRAPE:
	if (l->t == STR && empty_or_num(r)) {
	    ret = string_scrape_node(l, r, p);
	} else if (l->t != STR) {
	    node_type_error(t->t, 0, STR, NULL, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_BESSEL:
	/* functions taking one char, one scalar/series and one
	   matrix/series/scalar as args */
	if (l->t != STR) {
	    node_type_error(t->t, 1, STR, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(t->t, 2, NUM, m, p);
	} else if (r->t != NUM && r->t != SERIES && r->t != MAT) {
	    node_type_error(t->t, 3, NUM, r, p);
	} else {
	    ret = eval_bessel_func(l, m, r, p);
	}
	break;
    case F_BKFILT:
    case F_MOLS:
    case F_MPOLS:
    case F_MRLS:
    case F_FILTER:
    case F_MCOVG:
    case F_NRMAX:
    case F_LOESS:
    case F_GHK:
    case F_QUADTAB:
    case F_QLRPVAL:
    case F_BOOTCI:
    case F_BOOTPVAL:
    case F_MOVAVG:
    case F_DEFARRAY:
    case F_DEFBUNDLE:
    case F_DEFLIST:
    case HF_CLOGFI:
	/* built-in functions taking more than three args */
	ret = eval_nargs_func(t, p);
	break;
    case F_KSETUP:
    case F_KFILTER:
    case F_KSMOOTH:
    case F_KSIMUL:
    case F_KDSMOOTH:
	if (t->t == F_KSETUP || bundle_pointer_arg0(t)) {
	    ret = eval_kalman_bundle_func(t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_KSIMDATA:
	if (r->t == MAT) {
	    ret = kalman_data_node(l, r, p);
	} else {
	    node_type_error(t->t, 2, MAT, r, p);
	}
	break;
    case F_ISOCONV:
	ret = isoconv_node(t, p);
	break;
    case OVAR:
    case MVAR:
    case DMSL:
	/* variable "under" user-defined object, or last object */
	ret = object_var_node(t, p);
	break;
    case DMSTR:
	ret = dollar_str_node(t, NULL, p);
	break;
    case DVAR:
	/* dataset "dollar" variable */
	ret = dollar_var_node(t, p);
	break;
    case MDEF:
	/* matrix definition */
	ret = matrix_def_node(t, p);
	break;
    case F_OBSNUM:
    case F_ISDISCR:
    case F_ISNULL:
    case F_ISSTRING:
    case F_STRLEN:
    case F_NLINES:
    case F_REMOVE:
	if (l->t == STR) {
	    ret = object_status(l, t->t, p);
	} else {
	    node_type_error(t->t, 1, STR, l, p);
	}
	break;
    case F_TYPEOF:
    case F_EXISTS:
	if (l->t == STR) {
	    ret = object_status(l, t->t, p);
	} else {
	    ret = generic_typeof_node(l, p);
	}
	break;
    case F_NELEM:
	ret = n_elements_node(l, p);
	break;
    case F_INLIST:
	if (ok_list_node(l)) {
	    ret = in_list_node(l, r, p);
	} else {
	    node_type_error(t->t, 1, LIST, l, p);
	}
	break;
    case F_PDF:
    case F_CDF:
    case F_INVCDF:
    case F_CRIT:
    case F_PVAL:
    case F_RANDGEN:
    case F_MRANDGEN:
    case F_RANDGEN1:
    case F_URCPVAL:
	if (t->v.b1.b->t == FARGS) {
	    if (t->t == F_URCPVAL) {
		ret = eval_urcpval(t, p);
	    } else {
		ret = eval_pdist(t, p);
	    }
	} else {
	    node_type_error(t->t, 0, FARGS, t->v.b1.b, p);
	}
	break;
    case CON:
	/* built-in constant */
	ret = retrieve_const(t, p);
	break;
    case UFUN:
	ret = eval_ufunc(t, p);
	break;
#ifdef USE_RLIB
    case RFUN:
	ret = eval_Rfunc(t, p);
	break;
#endif
    case QUERY:
	ret = eval_query(t, p);
	break;
    case B_LCAT:
	/* list concatenation */
	if (ok_list_node(l) && ok_list_node(r)) {
	    ret = eval_lcat(l, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_SQUARE:
	if (ok_list_node(l) && empty_or_num(r)) {
	    ret = apply_list_func(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_WMEAN:
    case F_WVAR:
    case F_WSD:
	/* two lists -> series */
	if (ok_list_node(l) && ok_list_node(r)) {
	    ret = list_list_series_func(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_LINCOMB:
	/* list + matrix -> series */
	if (ok_list_node(l) && r->t == MAT) {
	    ret = lincomb_func(l, r, NULL, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_MLINCOMB:
	/* list + matrix + int -> series */
	if (ok_list_node(l) && m->t == MAT &&
	    (scalar_node(r) || r->t == STR)) {
	    ret = lincomb_func(l, m, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_ARGNAME:
    case F_BACKTICK:
    case F_STRSTRIP:
    case F_FIXNAME:
	if (l->t == STR) {
	    ret = single_string_func(l, t->t, p);
	} else if (t->t == F_ARGNAME && (uscalar_node(l) || useries_node(l))) {
	    ret = argname_from_uvar(l, p);
	} else {
	    node_type_error(t->t, 0, STR, l, p);
	}
	break;
    case F_READFILE:
	if (l->t == STR) {
	    ret = readfile_node(l, r, p);
	} else {
	    node_type_error(t->t, 1, STR, l, p);
	}
	break;
    case F_GETENV:
    case F_NGETENV:
	if (l->t == STR) {
	    ret = do_getenv(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, STR, l, p);
	}
	break;
    case F_OBSLABEL:
	if (l->t == NUM || l->t == MAT) {
	    ret = int_to_string_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, NUM, l, p);
	}
	break;
    case F_VARNAME:
	if (l->t == NUM || l->t == MAT || l->t == SERIES) {
	    ret = int_to_string_func(l, t->t, p);
	} else if (l->t == LIST) {
	    ret = list_to_string_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, NUM, l, p);
	}
	break;
    case F_VARNAMES:
	if (l->t == LIST) {
	    ret = list_to_string_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, LIST, l, p);
	}
	break;
    case F_VARNUM:
    case F_TOLOWER:
    case F_TOUPPER:
	if (l->t == STR) {
	    if (t->t == F_TOLOWER || t->t == F_TOUPPER) {
		ret = one_string_func(l, t->t, p);
	    } else {
		ret = varnum_node(l, p);
	    }
	} else {
	    node_type_error(t->t, 0, STR, l, p);
	}
	break;
    case F_STRSTR:
    case F_JSONGET:
    case F_XMLGET:
	if (l->t == STR && r->t == STR) {
	    ret = two_string_func(l, r, t->t, p);
	} else if (l->t == STR && r->t == ARRAY && t->t == F_XMLGET) {
	    ret = two_string_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, (l->t == STR)? 2 : 1,
			    STR, (l->t == STR)? r : l, p);
	}
	break;
    case F_STRSPLIT:
	if (l->t == STR && empty_or_num(r)) {
	    ret = strsplit_node(l, r, p);
	} else {
	    node_type_error(t->t, (l->t == STR)? 2 : 1,
			    (l->t == STR)? NUM : STR,
			    (l->t == STR)? r : l, p);
	}
	break;
    case F_GETLINE:
	if (l->t == STR && (r->t == STR || null_or_empty(r))) {
	    ret = getline_node(l, r, p);
	} else {
	    node_type_error(t->t, (l->t == STR)? 2 : 1,
			    STR, (l->t == STR)? r : l, p);
	}
	break;
    case F_ERRMSG:
	if (l->t == NUM || l->t == EMPTY) {
	    ret = errmsg_node(l, p);
	} else {
	    node_type_error(t->t, 0, NUM, l, p);
	}
	break;
    case F_ISODATE:
    case F_JULDATE:
	ret = isodate_node(l, r, t->t, p);
	break;
    case F_ATOF:
	if (l->t == STR) {
	    ret = atof_node(l, p);
	} else {
	    node_type_error(t->t, 0, STR, l, p);
	}
	break;
    case F_MPI_RECV:
	ret = mpi_transfer_node(l, NULL, NULL, t->t, p);
	break;
    case F_MPI_SEND:
    case F_BCAST:
    case F_ALLREDUCE:
	if (t->t == F_ALLREDUCE && r->t != STR) {
	    node_type_error(t->t, 2, STR, r, p);
	} else {
	    ret = mpi_transfer_node(l, r, NULL, t->t, p);
	}
	break;
    case F_REDUCE:
    case F_SCATTER:
	if (m->t != STR) {
	    node_type_error(t->t, 2, STR, m, p);
	} else if (!empty_or_num(r)) {
	    node_type_error(t->t, 3, NUM, r, p);
	} else {
	    ret = mpi_transfer_node(l, m, r, t->t, p);
	}
	break;
    case F_GENSERIES:
	ret = gen_series_node(l, r, p);
	break;
    case F_ARRAY:
	ret = gen_array_node(l, p);
	break;
    case F_STRVALS:
	if (!useries_node(l)) {
	    node_type_error(t->t, 0, USERIES, l, p);
	} else {
	    ret = get_series_stringvals(l, p);
	}
	break;
    case F_STRINGIFY:
	if (!useries_node(l) || r->t != ARRAY) {
	    int l_ok = useries_node(l);

	    node_type_error(t->t, l_ok ? 2 : 1,
			    l_ok ? ARRAY : USERIES,
			    l_ok ? r : l, p);
	} else {
	    ret = stringify_series(l, r, p);
	}
	break;
    default:
	printf("eval: weird node %s (t->t = %d)\n", getsymb(t->t),
	       t->t);
	p->err = E_PARSE;
	break;
    }

 finish:

    if (!p->err && ret != NULL && ret != t && is_aux_node(ret)) {
	p->err = attach_aux_node(t, ret, p);
    }

 bailout:

#if EDEBUG
    fprintf(stderr, "eval (t->t = %03d, %s): returning NODE at %p, err %d\n",
	    t->t, getsymb(t->t), (void *) ret, p->err);
    if (t->t == SERIES)
	fprintf(stderr, " (SERIES node, xvec at %p, vnum = %d)\n",
		(void *) t->v.xvec, t->vnum);
#endif

    return ret;
}

#if !AUX_NODES_DEBUG

/* non-debugging variant: easier to see what's going on */

static inline int attach_aux_node (NODE *t, NODE *ret, parser *p)
{
    if (t->aux == NULL) {
	t->aux = ret;
	ret->refcount += 1;
    } else if (t->aux != ret) {
	if (t->t == QUERY || t->t == B_AND || t->t == B_OR) {
	    /* OK to switch aux node in these cases */
	    free_node(t->aux, p);
	    t->aux = ret;
	    ret->refcount += 1;
	} else if (is_proxy_node(t->aux)) {
	    /* an extension to the above */
	    free_node(t->aux, p);
	    t->aux = ret;
	    ret->refcount += 1;
	} else {
	    /* otherwise if we're trying to switch aux node,
	       something must have gone wrong
	    */
	    fprintf(stderr, "! node %s already has aux node %s attached\n",
		    getsymb(t->t), getsymb(t->aux->t));
	    return E_DATA;
	}
    }

    return 0;
}

#else

/* variant with lots of debugging spew */

static inline int attach_aux_node (NODE *t, NODE *ret, parser *p)
{
    if (t->aux == NULL) {
	fprintf(stderr, "++ attach aux node %p (%s) to node %p (%s)\n",
		ret, getsymb(ret->t), t, getsymb(t->t));
	if (ret->refcount > 0) {
	    fprintf(stderr, "   note: refcount on %p = %d\n",
		    (void *) ret, ret->refcount);
	}
	t->aux = ret;
	ret->refcount += 1;
    } else if (t->aux != ret) {
	if (t->t == QUERY || t->t == B_AND || t->t == B_OR) {
	    /* the result node may switch in these cases (only?) */
	    fprintf(stderr, "boolean: freeing %p\n", (void *) t->aux);
	    free_node(t->aux, p);
	    fprintf(stderr, "boolean: attaching %p\n", (void *) ret);
	    t->aux = ret;
	    ret->refcount += 1;
	} else {
	    fprintf(stderr, "!! node %p (%s) has aux node %p (%s) attached\n"
		    "   not attaching ret at %p (%s)\n", t, getsymb(t->t),
		    t->aux, getsymb(t->aux->t), ret, getsymb(ret->t));
	    return E_DATA;
	}
    }

    return 0;
}

#endif

static int more_input (const char *s)
{
    while (*s) {
	if (!isspace((unsigned char) *s)) {
	    return 1;
	}
	s++;
    }

    return 0;
}

/* get the next input character for the lexer */

int parser_getc (parser *p)
{
#if EDEBUG > 1
    fprintf(stderr, "parser_getc: src='%s'\n", p->point);
#endif

    p->ch = 0;

    if (more_input(p->point)) {
	p->ch = *p->point;
	p->point += 1;
    }

#if EDEBUG > 1
    if (p->ch) {
	fprintf(stderr, "parser_getc: returning '%c'\n", p->ch);
    }
#endif

    return p->ch;
}

/* advance the read position by n characters */

void parser_advance (parser *p, int n)
{
    p->point += n;
    p->ch = *p->point;
    p->point += 1;
}

/* throw back the last-read character */

void parser_ungetc (parser *p)
{
    p->point -= 1;
    p->ch = *(p->point - 1);
}

/* Look ahead for the first occurrence of a given character in
   the remaining input stream; return its 0-based index or
   -1 if not found.
*/

int parser_char_index (parser *p, int c)
{
    int i;

    for (i=0; p->point[i] != '\0'; i++) {
	if (p->point[i] == c) {
	    return i;
	}
    }

    return -1;
}

/* Look ahead to the next non-space character in the
   parser stream and return it; if @skip then start at
   offset 1 beyond the current p->point.
*/

int parser_next_nonspace_char (parser *p, int skip)
{
    int i, offset = skip ? 1 : 0;

    if (*p->point == '\0') {
	return 0;
    }

    for (i=offset; p->point[i] != '\0'; i++) {
	if (!isspace(p->point[i])) {
	    return p->point[i];
	}
    }

    return 0;
}

/* gets the next character in the input stream */

int parser_next_char (parser *p)
{
    return p->point[0];
}

/* For error reporting: print the input up to the current
   parse point, unless it's not valid UTF-8. Return 0
   if the input is printed OK, otherwise non-zero.
*/

int parser_print_input (parser *p)
{
    int len = p->point - p->input;
    char *s = gretl_strndup(p->input, len);
    int err = 0;

    if (s != NULL) {
	if (g_utf8_validate(s, -1, NULL)) {
	    pprintf(p->prn, "> %s\n", s);
	} else {
	    err = 1;
	}
	free(s);
    } else {
	err = 1;
    }

    return err;
}

/* "pretty print" syntactic nodes and symbols */

static void printsymb (int symb, const parser *p)
{
    pputs(p->prn, getsymb(symb));
}

static void printnode (NODE *t, parser *p, int value)
{
    if (t == NULL) {
	pputs(p->prn, "NULL");
    } else if (!value && useries_node(t)) {
	pprintf(p->prn, "%s", p->dset->varname[t->vnum]);
    } else if (!value && uscalar_node(t)) {
	pprintf(p->prn, "%s", t->vname);
    } else if (t->t == NUM) {
	if (na(t->v.xval)) {
	    pputs(p->prn, "NA");
	} else {
	    pprintf(p->prn, "%.8g", t->v.xval);
	}
    } else if (t->t == SERIES) {
	const double *x = t->v.xvec;
	int i, j = 1;

	if (p->lh.vnum > 0 && p->lh.vnum < p->dset->v) {
	    pprintf(p->prn, "%s\n", p->dset->varname[p->lh.vnum]);
	}

	for (i=p->dset->t1; i<=p->dset->t2; i++, j++) {
	    if (na(x[i])) {
		pputs(p->prn, "NA");
	    } else {
		pprintf(p->prn, "%g", x[i]);
	    }
	    if (j % 8 == 0) {
		pputc(p->prn, '\n');
	    } else if (i < p->dset->t2) {
		pputc(p->prn, ' ');
	    }
	}
    } else if (t->t == MAT) {
	if (t->vname != NULL) {
	    pputs(p->prn, t->vname);
	} else {
	    gretl_matrix_print_to_prn(t->v.m, NULL, p->prn);
	}
    } else if (t->t == BUNDLE) {
	gretl_bundle_print(t->v.b, p->prn);
    } else if (t->t == DBUNDLE) {
	pputs(p->prn, bvarname(t->v.idnum));
    } else if (t->t == ARRAY) {
	gretl_array_print(t->v.a, p->prn);
    } else if (t->t == UOBJ) {
	pprintf(p->prn, "%s", t->v.str);
    } else if (t->t == DVAR) {
	pputs(p->prn, dvarname(t->v.idnum));
    } else if (t->t == MVAR) {
	pputs(p->prn, mvarname(t->v.idnum));
    } else if (t->t == CON) {
	pputs(p->prn, constname(t->v.idnum));
    } else if (t->t == DUM) {
	pputs(p->prn, dumname(t->v.idnum));
    } else if (binary_op(t->t)) {
	pputc(p->prn, '(');
	printnode(t->v.b2.l, p, 0);
	printsymb(t->t, p);
	printnode(t->v.b2.r, p, 0);
	pputc(p->prn, ')');
    } else if (t->t == MSL || t->t == DMSL) {
	printnode(t->v.b2.l, p, 0);
	pputc(p->prn, '[');
	printnode(t->v.b2.r, p, 0);
	pputc(p->prn, ']');
    } else if (t->t == MSLRAW) {
	pputs(p->prn, "MSLRAW");
    } else if (t->t == SUBSL) {
	pputs(p->prn, "SUBSL");
    } else if (t->t == OVAR) {
	printnode(t->v.b2.l, p, 0);
	pputc(p->prn, '.');
	printnode(t->v.b2.r, p, 0);
    } else if (func1_symb(t->t)) {
	printsymb(t->t, p);
	pputc(p->prn, '(');
	printnode(t->v.b1.b, p, 0);
	pputc(p->prn, ')');
    } else if (unary_op(t->t)) {
	printsymb(t->t, p);
	printnode(t->v.b1.b, p, 0);
    } else if (func2_symb(t->t)) {
	printsymb(t->t, p);
	pputc(p->prn, '(');
	printnode(t->v.b2.l, p, 0);
	if (t->v.b2.r->t != EMPTY) {
	    pputc(p->prn, ',');
	}
	printnode(t->v.b2.r, p, 0);
	pputc(p->prn, ')');
    } else if (t->t == STR) {
	pprintf(p->prn, "%s", t->v.str);
    } else if (t->t == PTR) {
	pprintf(p->prn, "%s", t->vname);
    } else if (t->t == MDEF) {
	pprintf(p->prn, "{ MDEF }");
    } else if (t->t == DMSTR || t->t == UFUN) {
	printnode(t->v.b2.l, p, 0);
	pputc(p->prn, '(');
	printnode(t->v.b2.r, p, 0);
	pputc(p->prn, ')');
    } else if (t->t == LISTVAR) {
	pprintf(p->prn, "%s.%s", t->v.b2.l->v.str, t->v.b2.r->v.str);
    } else if (t->t == LIST) {
	pputs(p->prn, "LIST");
    } else if (t->t == LAG) {
	pputs(p->prn, "LAG");
    } else if (t->t != EMPTY) {
	pputs(p->prn, "weird tree - ");
	printsymb(t->t, p);
    }
}

/* which modified assignment operators of the type '+='
   will we accept, when generating various types of
   result? */
#define ok_matrix_op(o) (o == B_ASN  || o == B_DOTASN || \
			 o == B_ADD  || o == B_SUB || \
			 o == B_MUL  || o == B_DIV || \
			 o == B_HCAT || o == B_VCAT)
#define ok_list_op(o) (o == B_ASN || o == B_ADD || o == B_SUB)
#define ok_string_op(o) (o == B_ASN || o == B_ADD || \
			 o == B_HCAT || o == INC)
#define ok_array_op(o) (o == B_ASN || o == B_ADD)

struct mod_assign {
    int c;
    int op;
};

/* supported "inflections" of assignment */

struct mod_assign m_assign[] = {
    { '+', B_ADD },
    { '-', B_SUB },
    { '*', B_MUL },
    { '/', B_DIV },
    { '%', B_MOD},
    { '^', B_POW },
    { '~', B_HCAT },
    { '|', B_VCAT },
    { '.', B_DOTASN },
    { 0, 0}
};

/* read operator from formula: this is either
   simple assignment or something like '+=' */

static int get_op (char *s)
{
    if (s[0] == '=') {
	s[1] = '\0';
	return B_ASN;
    }

    if (!strcmp(s, "++")) {
	return INC;
    }

    if (!strcmp(s, "--")) {
	return DEC;
    }

    if (s[1] == '=') {
	int i;

	for (i=0; m_assign[i].c; i++) {
	    if (s[0] == m_assign[i].c) {
		return m_assign[i].op;
	    }
	}
    }

    return 0;
}

static char *get_opstr (int op)
{
    static char opstr[3] = {0};

    if (op == B_ASN) {
	return "=";
    } else if (op == INC) {
	return "++";
    } else if (op == DEC) {
	return "--";
    } else {
	int i;

	for (i=0; m_assign[i].c; i++) {
	    if (op == m_assign[i].op) {
		opstr[0] = m_assign[i].c;
		opstr[1] = '=';
		return opstr;
	    }
	}
	return "??";
    }
}

/* implement the declaration of new variables */

static void do_declaration (parser *p)
{
    char **S = NULL;
    int i, v, n;

    n = check_declarations(&S, p);

    if (n == 0) {
	return;
    }

    for (i=0; i<n && !p->err; i++) {
	if (S[i] != NULL) {
	    if (p->targ == SERIES) {
		p->err = dataset_add_NA_series(p->dset, 1);
		if (!p->err) {
		    v = p->dset->v - 1;
		    strcpy(p->dset->varname[v], S[i]);
		}
	    } else {
		GretlType type = 0;

		if (p->targ == MAT) {
		    type = GRETL_TYPE_MATRIX;
		} else if (p->targ == NUM) {
		    type = GRETL_TYPE_DOUBLE;
		} else if (p->targ == STR) {
		    type = GRETL_TYPE_STRING;
		} else if (p->targ == BUNDLE) {
		    type = GRETL_TYPE_BUNDLE;
		} else if (p->targ == LIST) {
		    type = GRETL_TYPE_LIST;
		} else if (p->targ == ARRAY) {
		    type = p->lh.gtype;
		} else {
		    p->err = E_DATA;
		}
		if (!p->err) {
		    p->err = create_user_var(S[i], type);
		}
	    }
	}
    }

    strings_array_free(S, n);
}

/* The expression supplied for evaluation does not contain an '=':
   can we interpret it as an implicit request to print the value
   of an existing variable?
*/

static void parser_try_print (parser *p, const char *s)
{
    if (p->lh.t != 0 && p->lh.expr == NULL) {
	p->flags |= P_DISCARD;
	p->point = s;
    } else {
	p->err = E_EQN;
    }
}

/* Here we try to parse out the LHS of the statement
   and also the operator. If we find a unitary LHS
   (simply an indentifier) we write it into p->lh.name,
   but if we find a compound LHS (such as a sub-matrix
   specification) we save it as p->lh.expr. The
   content of @ps is advanced to the first position
   beyond the operator.
*/

static int extract_lhs_and_op (const char **ps, parser *p,
			       char *opstr)
{
    const char *s = *ps;
    int quoted = 0;
    int i, n = 0;
    int err = 0;

#if LHDEBUG
    fprintf(stderr, "extract: input='%s'\n", s);
#endif

    if (p->targ != UNK && strchr(s, '=') == NULL) {
	/* we got a type specification but no assignment,
	   so should be variable declaration(s) ?
	*/
	p->flags |= P_DECL;
	p->lh.expr = gretl_strdup(s);
	goto done;
    }

    /* Count bytes preceding first unquoted '='. Note that
       the "unquoted" condition is required only because
       a string-literal bundle key might contain an equals
       sign, as in b["foo=bar"] = ...
    */
    for (i=0; s[i] != '\0'; i++) {
	if (s[i] == '"') {
	    quoted = !quoted;
	} else if (!quoted && s[i] == '=') {
	    break;
	}
	n++;
    }

    if (n > 0) {
	char *lhs = NULL;
	int lhlen = n;

	if (s[n] == '=') {
	    /* we actually reached an '=' */
	    if (strspn(s + n - 1, "+-*/%^~|.") == 1) {
		/* preceded by a modifier: inflected assignment */
		lhlen--;
		opstr[0] = s[n-1];
		opstr[1] = '=';
	    } else {
		/* no: straight assignment */
		opstr[0] = '=';
	    }
	    n++; /* add one for '=' */
	}

	if (lhlen > 0) {
	    lhs = gretl_strndup(s, lhlen);
	    tailstrip(lhs);
	    lhlen = strlen(lhs);
	}

	if (opstr[0] == '\0' && lhlen > 2) {
	    /* check for postfix operator */
	    char *test = lhs + lhlen - 2;

	    if (!strcmp(test, "++") || !strcmp(test, "--")) {
		strcpy(opstr, test);
		*test = '\0';
		lhlen -= 2;
	    }
	}
	    
	if (lhlen > 0) {
	    if (lhlen == gretl_namechar_spn(lhs)) {
		/* a straight identifier? */
		if (lhlen >= VNAMELEN) {
		    pprintf(p->prn, _("'%s': name is too long (max %d characters)\n"),
			    lhs, VNAMELEN - 1);
		    err = E_PARSE;
		} else {
		    strcpy(p->lh.name, lhs);
		}
	    } else if ((p->flags & P_PRIV) && (*lhs == '$' || *lhs == '_') &&
		       gretl_namechar_spn(lhs + 1) == lhlen - 1) {
		/* "private" genr of the form $foo=expr or _foo=expr */
		strcpy(p->lh.name, lhs);
	    } else {
		/* treat as an expression to be evaluated */
		p->lh.expr = lhs;
		lhs = NULL; /* protect against freeing */
	    }
	} else {
	    /* nothing relevant found */
	    err = E_PARSE;
	}

	if (!err && opstr[0] != '\0') {
	    p->op = get_op(opstr);
	}

	free(lhs);
	*ps = s + n;
    }

 done:

#if LHDEBUG
    fprintf(stderr, "extract: name='%s', expr='%s', op='%s', err=%d,\n s='%s'\n",
	    p->lh.name, p->lh.expr ? p->lh.expr : "nil", opstr, err, *ps);
#endif    

    return err;
}

static void maybe_do_type_errmsg (const char *name, int t)
{
    const char *tstr = NULL;

    if (t == NUM) {
	tstr = "scalar";
    } else if (t == SERIES) {
	tstr = "series";
    } else if (t == MAT) {
	tstr = "matrix";
    } else if (t == STR) {
	tstr = "string";
    } else if (t == BUNDLE) {
	tstr = "bundle";
    } else if (t == LIST) {
	tstr = "list";
    } else if (t == ARRAY) {
	tstr = "array";
    }

    if (tstr != NULL) {
	if (name != NULL && *name != '\0') {
	    gretl_errmsg_sprintf(_("The variable %s is of type %s, "
				   "not acceptable in context"),
				 name, tstr);
	} else {
	    gretl_errmsg_sprintf(_("A variable of type %s is not "
				   "acceptable in context"), tstr);
	}
    }
}

static int overwrite_type_check (parser *p)
{
    int err = 0;

    if (p->targ != p->lh.t) {
	/* don't overwrite one type with another */
	maybe_do_type_errmsg(p->lh.name, p->lh.t);
	err = E_TYPES;
    }

    return err;
}

static int overwrite_const_check (const char *s)
{
    return object_is_const(s) ? overwrite_err(s) : 0;
}

/* Check that we're not trying to modify a const object
   via a compound LHS expression; and while we're at
   it, check whether we should be generating a list
   (an element of an array of lists).
*/

static int compound_const_check (NODE *lhs, parser *p)
{
    NODE *n = lhs;
    int i = 0, err = 0;

    while (n->t == MSL || n->t == OBS || n->t == BMEMB || n->t == OSL) {
	n = n->v.b2.l;
	if (i == 0 && lhs->t == OSL && n->t == ARRAY) {
	    if (gretl_array_get_type(n->v.a) == GRETL_TYPE_LISTS) {
		p->flags |= P_LISTDEF;
	    }
	}
	i++;
    }

    /* do we have a const object at the tip of the tree? */
    if (n->vname != NULL) {
	err = overwrite_const_check(n->vname);
    }

    return err;
}

static int ok_array_decl (parser *p, const char *s)
{
    p->lh.gtype = 0;

    if (!strncmp(s, "strings ", 8)) {
	p->lh.gtype = GRETL_TYPE_STRINGS;
    } else if (!strncmp(s, "matrices ", 9)) {
	p->lh.gtype = GRETL_TYPE_MATRICES;
    } else if (!strncmp(s, "bundles ", 8)) {
	p->lh.gtype = GRETL_TYPE_BUNDLES;
    } else if (!strncmp(s, "lists ", 6)) {
	p->lh.gtype = GRETL_TYPE_LISTS;
    }

    return p->lh.gtype != 0;
}

/* Given an existing LHS variable, whose type is recorded in
   p->lh.t, check that the specified operator is supported
   for the type. Return error code if not.
*/

static int check_operator_validity (parser *p, const char *opstr)
{
    if (p->lh.t == MAT && !ok_matrix_op(p->op)) {
	/* matrices: we accept only a limited range of
	   modified assignment operators */
	gretl_errmsg_sprintf(_("'%s' : not implemented for matrices"), opstr);
	return E_PARSE;
    } else if (p->lh.t == LIST && !ok_list_op(p->op)) {
	/* lists: same story as matrices */
	gretl_errmsg_sprintf(_("'%s' : not implemented for lists"), opstr);
	return E_PARSE;
    } else if (p->lh.t == STR && !ok_string_op(p->op)) {
	/* strings: ditto */
	gretl_errmsg_sprintf(_("'%s' : not implemented for strings"), opstr);
	return E_PARSE;
    } else if (p->lh.t == ARRAY && !ok_array_op(p->op)) {
	/* arrays: ditto */
	gretl_errmsg_sprintf(_("'%s' : not implemented for arrays"), opstr);
	return E_PARSE;
    } else if (p->lh.t == BUNDLE && p->op != B_ASN) {
	/* bundles: no modified assignment (yet) */
	gretl_errmsg_sprintf(_("'%s' : not implemented for this type"), opstr);
	return E_PARSE;
    } else if (p->lh.t != MAT && (p->op == B_VCAT || p->op == B_DOTASN)) {
	/* vertical concat: only OK for matrices */
	gretl_errmsg_sprintf(_("'%s' : only defined for matrices"), opstr);
	return E_PARSE;
    } else if (p->lh.t != MAT && p->lh.t != STR && p->op == B_HCAT) {
	/* horizontal concat: only OK for matrices, strings */
	gretl_errmsg_sprintf(_("'%s' : not implemented for this type"), opstr);
	return E_PARSE;
    } else if (p->lh.t == SERIES && is_string_valued(p->dset, p->lh.vnum)) {
	/* string-valued series: do not overwrite wholesale */
	gretl_errmsg_set("Cannot overwrite entire string-valued series");
	return E_TYPES;
    }

    /* otherwise OK? */
    return 0;
}

/* Do we have an inline type specification preceding the
   statement proper? In most cases we shouldn't, since it
   will already have been handled by the tokenizer (and
   the type will now be recorded in p->targ). But we allow
   for finding a typespec here in case of "genrs" within
   nls/mle/gmm blocks, where the statement bypasses the
   regular tokenizer. (FIXME?)
*/

static void check_for_inline_typespec (const char **ps, parser *p)
{
    const char *s = *ps;
    
    if (!strncmp(s, "scalar ", 7)) {
	p->targ = NUM;
	s += 7;
    } else if (!strncmp(s, "series ", 7)) {
	p->targ = SERIES;
	s += 7;
    } else if (!strncmp(s, "matrix ", 7)) {
	p->targ = MAT;
	s += 7;
    } else if (!strncmp(s, "list ", 5)) {
	p->targ = LIST;
	s += 5;
    } else if (!strncmp(s, "string ", 7)) {
	p->targ = STR;
	s += 7;
    } else if (!strncmp(s, "bundle ", 7)) {
	p->targ = BUNDLE;
	s += 7;
    } else if (ok_array_decl(p, s)) {
	p->targ = ARRAY;
	s += strcspn(s, " ") + 1;
    }

    /* advance pointer */
    *ps = s;
}

/* Check @p->lh.name for the name of an existing series or 
   user_var of some kind. If found, record the relevant
   info in p->lh.t, and also either p->lh.vnum (series) or
   p->lh.uv (other types).
*/

static int check_existing_lhs_type (parser *p, int *newvar)
{
    user_var *uvar;
    int v, err = 0;

    v = current_series_index(p->dset, p->lh.name);
    if (v >= 0) {
	p->lh.vnum = v;
	p->lh.t = SERIES;
	*newvar = 0;
	return 0;
    }
    
    uvar = get_user_var_by_name(p->lh.name);

    if (uvar != NULL) {
	GretlType vtype = uvar->type;

	p->lh.uv = uvar;
	*newvar = 0;

	if (vtype == GRETL_TYPE_MATRIX) {
	    p->lh.t = MAT;
	} else if (vtype == GRETL_TYPE_DOUBLE) {
	    if (uvar->flags & UV_NODECL) {
		if (p->targ == UNK) {
		    p->flags |= P_NODECL;
		} else {
		    uvar->flags &= ~UV_NODECL;
		}
	    }
	    p->lh.t = NUM;
	} else if (vtype == GRETL_TYPE_LIST) {
	    p->lh.t = LIST;
	} else if (vtype == GRETL_TYPE_STRING) {
	    p->lh.t = STR;
	} else if (vtype == GRETL_TYPE_BUNDLE) {
	    p->lh.t = BUNDLE;
	} else if (vtype == GRETL_TYPE_ARRAY) {
	    p->lh.gtype = gretl_array_get_type(uvar->ptr);
	    p->lh.t = ARRAY;
	}
    }

    return err;
}

/* pre-process a "genr" statement */

static void gen_preprocess (parser *p, int flags)
{
    const char *s = p->input;
    char opstr[3] = {0};
    int newvar = 1;

    while (isspace(*s)) s++;

    /* skip leading command word, if any */
    if (!strncmp(s, "genr ", 5)) {
	s += 5;
    } else if (!strncmp(s, "print ", 6)) {
	/* allow this within (e.g.) mle block */
	p->flags |= P_DISCARD;
	s += 6;
    }

    while (isspace(*s)) s++;

    if (p->targ == UNK) {
	check_for_inline_typespec(&s, p);
    } else if (gretl_array_type(p->targ)) {
	/* record a plural type spec such as "matrices" under
	   the "lh" member of @p.
	*/
	p->lh.gtype = p->targ;
	p->targ = ARRAY;
    }

    /* check for types that cannot be generated in the
       absence of a dataset */
    if ((p->targ == SERIES || p->targ == LIST) &&
	(p->dset == NULL || p->dset_n == 0)) {
	no_data_error(p);
	return;
    }

    if (p->flags & P_DISCARD) {
	/* doing a simple "eval" */
	p->point = s;
	return;
    }

    /* extract LHS expression and operator, and test for a declaration */
    p->err = extract_lhs_and_op(&s, p, opstr);
    if (p->err || (p->flags & P_DECL)) {
	return;
    }

    /* record next read position */
    p->point = s;

    if (p->lh.expr != NULL) {
	/* create syntax tree for the LHS expression */
	const char *savepoint = p->point;

	p->point = p->lh.expr;
	p->ch = parser_getc(p);
	lex(p);
	p->lhtree = expr(p);
#if LHDEBUG
	fprintf(stderr, "parsed lhtree, err=%d\n", p->err);
	print_tree(p->lhtree, p, 0);
#endif
	p->point = savepoint;
	p->ch = 0;
	if (!p->err) {
	    p->err = compound_const_check(p->lhtree, p);
	}
	if (p->err) {
	    return;
	} else {
	    goto get_rhs;
	}
    }

    /* find out if the LHS var already exists, and if
       so, what type it is */
    if (!p->err) {
	p->err = check_existing_lhs_type(p, &newvar);
    }

#if LHDEBUG
    fprintf(stderr, "newvar=%d, err=%d\n", newvar, p->err);
#endif

    if (p->err) {
	return;
    }    

    if (newvar) {
	/* new variable: check name for legality */
	if (!(flags & P_PRIV)) {
	    p->err = check_varname(p->lh.name);
	}
    } else {
	/* pre-existing var: check for const-ness */
	p->err = overwrite_const_check(p->lh.name);
    }

    if (p->err) {
	return;
    }

    if (p->lh.t != 0) {
	if (p->targ == UNK) {
	    /* when a result type is not specified, set this
	       from existing LHS variable, if present
	    */
	    p->targ = p->lh.t;
	} else if (overwrite_type_check(p)) {
	    /* don't overwrite one type with another */
	    p->err = E_TYPES;
	    return;
	}
    }

 get_rhs:

    /* advance past white space */
    while (isspace(*s)) s++;
    p->point = p->rhs = s;

    if (p->lh.expr != NULL) {
	goto alt_set_targ;
    }

    /* expression ends here with no operator: a call to print? */
    if (*s == '\0' && p->op == 0) {
	parser_try_print(p, p->lh.name);
	return;
    }

    /* if the LHS variable does not already exist, then
       we can't do '+=' or anything of that sort, only
       simple assignment, B_ASN
    */
    if (newvar && p->op != B_ASN) {
	undefined_symbol_error(p->lh.name, p);
	return;
    }

    if (p->op) {
	p->err = check_operator_validity(p, opstr);
	if (p->err) {
	    return;
	}
    }

 alt_set_targ:

    if (p->targ == UNK && *p->rhs == '{') {
	/* if the target type is still unknown and the RHS
	   expression is wrapped in '{' and '}', make the target
	   a matrix
	*/
	p->targ = MAT;
    } else if (p->targ == LIST) {
	/* flag list target to parses */
	p->flags |= P_LISTDEF;
    }
}

/* tests for saving variable */

static int matrix_may_be_masked (const gretl_matrix *m, int n,
				 parser *p)
{
    int mt1 = gretl_matrix_get_t1(m);
    int mt2 = gretl_matrix_get_t2(m);
    int fullrows = mt2 - mt1 + 1;
    int nobs = get_matrix_mask_nobs();

    if (n == nobs && fullrows > n) {
	p->flags |= P_MMASK;
	return 1;
    } else {
	return 0;
    }
}

/* check whether a matrix result can be assigned to a series
   on return */

static int series_compatible (const gretl_matrix *m, parser *p)
{
    int n = gretl_vector_get_length(m);
    int mt2 = gretl_matrix_get_t2(m);
    int T = sample_size(p->dset);
    int ok = 0;

    if (mt2 > 0) {
	int mt1 = gretl_matrix_get_t1(m);

	if (n == mt2 - mt1 + 1) {
	    /* sample is recorded on matrix */
	    ok = 1;
	} else if (matrix_may_be_masked(m, n, p)) {
	    ok = 1;
	}
    } else if (n == T) {
	/* length matches current sample */
	ok = 1;
    } else if (n == p->dset->n) {
	/* length matches full series length */
	ok = 1;
    } else if (n == 1) {
	/* scalar: can be expanded */
	ok = 1;
    }

    return ok;
}

#define GENMISS_VERBOSE 0

/* below: note that the node we're checking is p->ret, which may
   differ in type from the target to which it will be
   assigned/converted
*/

static void gen_check_errvals (parser *p)
{
    NODE *r = p->ret;

    if (r == NULL || (r->t == SERIES && r->v.xvec == NULL)) {
	return;
    }

    if (p->targ == MAT) {
	/* the matrix target is handled separately */
	return;
    }

    if (r->t == NUM) {
	if (!isfinite(r->v.xval)) {
#if SCALARS_ENSURE_FINITE
# if GENMISS_VERBOSE
	    fprintf(stderr, "scalar return = %g\n", r->v.xval);
# endif
	    r->v.xval = NADBL;
	    set_gretl_warning(W_GENMISS);
#else
	    set_gretl_warning(W_GENNAN);
#endif
	}
    } else if (r->t == SERIES) {
	int t;

	for (t=p->dset->t1; t<=p->dset->t2; t++) {
	    if (!isfinite(r->v.xvec[t])) {
#if SERIES_ENSURE_FINITE
# if GENMISS_VERBOSE
		fprintf(stderr, "series return: xvec[%d] = %g\n", t, r->v.xvec[t]);
# endif
		r->v.xvec[t] = NADBL;
		set_gretl_warning(W_GENMISS);
#else
		set_gretl_warning(W_GENNAN);
#endif
		break;
	    }
	}
    } else if (r->t == MAT) {
	/* note: but p->targ != MAT */
	const gretl_matrix *m = r->v.m;
	int i, k = gretl_matrix_rows(m) * gretl_matrix_cols(m);

	if (p->targ == NUM && k == 1) {
	    if (!isfinite(m->val[0])) {
#if SCALARS_ENSURE_FINITE
		m->val[0] = NADBL;
		set_gretl_warning(W_GENMISS);
#else
		set_gretl_warning(W_GENNAN);
#endif
	    }
	} else {
#if MATRIX_NA_CHECK
	    /* convert any NAs to NaNs */
	    for (i=0; i<k; i++) {
		if (na(m->val[i])) {
		    m->val[i] = M_NA;
		    set_gretl_warning(W_GENNAN);
		} else if (!isfinite(m->val[i])) {
		    set_gretl_warning(W_GENNAN);
		}
	    }
#endif
	}
    }
}

static gretl_matrix *series_to_matrix (const double *x,
				       parser *p,
				       int *prechecked)
{
    int i, n = sample_size(p->dset);
    gretl_matrix *v;
    double xi;

    v = gretl_column_vector_alloc(n);

    if (v == NULL) {
	p->err = E_ALLOC;
    } else {
	for (i=0; i<n; i++) {
	    xi = x[i + p->dset->t1];
	    if (na(xi)) {
		set_gretl_warning(W_GENNAN);
		xi = M_NA;
	    }
	    v->val[i] = xi;
	}
	gretl_matrix_set_t1(v, p->dset->t1);
	gretl_matrix_set_t2(v, p->dset->t2);
    }

    return v;
}

static gretl_matrix *retrieve_matrix_result (parser *p,
					     int *prechecked)
{
    NODE *r = p->ret;
    gretl_matrix *m = NULL;

#if EDEBUG
    fprintf(stderr, "retrieve_matrix_result: r->t = %d\n", r->t);
#endif

    if (r->t == NUM) {
	m = gretl_matrix_from_scalar(r->v.xval);
	if (m == NULL) {
	    p->err = E_ALLOC;
	} else if (na(r->v.xval)) {
	    set_gretl_warning(W_GENNAN);
	    if (prechecked != NULL) {
		*prechecked = 1;
	    }
	}
    } else if (r->t == SERIES) {
	m = series_to_matrix(r->v.xvec, p, prechecked);
    } else if (r->t == LIST) {
	m = gretl_list_to_vector(r->v.ivec, &p->err);
    } else if (r->t == MAT && is_tmp_node(r)) {
	/* result matrix is newly allocated, steal it */
#if EDEBUG
	fprintf(stderr, "matrix result (%p) is tmp, stealing it\n",
		(void *) r->v.m);
#endif
	m = r->v.m;
	r->v.m = NULL; /* avoid double-freeing */
    } else if (r->t == MAT) {
	/* r->v.m is an existing user matrix (or bundled matrix):
	   must make a copy to keep pointers distinct
	*/
	m = gretl_matrix_copy(r->v.m);
#if EDEBUG
	fprintf(stderr, "matrix result (%p) is pre-existing, copied to %p\n",
		(void *) r->v.m, (void *) m);
#endif
	if (m == NULL) {
	    p->err = E_ALLOC;
	}
	if (!p->err && prechecked != NULL) {
	    *prechecked = 1;
	}
    } else {
	fprintf(stderr, "Looking for matrix, but r->t = %d\n", r->t);
	p->err = E_TYPES;
    }

    return m;
}

/* Check to see if the existing LHS matrix is of the
   same dimensions as the RHS result */

static int LHS_matrix_reusable (parser *p, gretl_matrix **pm)
{
    gretl_matrix *m = gen_get_lhs_var(p, GRETL_TYPE_MATRIX);
    int ok = 0;

    if (m == NULL) {
	return 0;
    } else if (p->ret->t == NUM) {
	ok = (m->rows == 1 && m->cols == 1);
    } else if (p->ret->t == SERIES) {
	int T = sample_size(p->dset);

	ok = (m->rows == T && m->cols == 1);
    } else if (p->ret->t == MAT) {
	gretl_matrix *rm = p->ret->v.m;

	ok = (rm != NULL &&
	      m->rows == rm->rows &&
	      m->cols == rm->cols);
    }

    *pm = m;

    return ok;
}

/* Generating a matrix, and there's a pre-existing LHS matrix:
   we re-use the left-hand side matrix if possible.
*/

static gretl_matrix *assign_to_matrix (parser *p, int *prechecked)
{
    gretl_matrix *m = NULL;
    double x;

    if (LHS_matrix_reusable(p, &m)) {
	/* The result is of the same dimensions as the LHS matrix:
	   this means that we don't need to construct an RHS
	   matrix if it doesn't already exist as such, nor do we
	   need to copy it if it does already exist.
	*/
#if EDEBUG
	fprintf(stderr, "assign_to_matrix: reusing LHS\n");
#endif
	if (p->ret->t == NUM) {
	    /* using RHS scalar */
	    x = p->ret->v.xval;
	    m->val[0] = na(x)? M_NA : x;
	    *prechecked = 1;
	} else if (p->ret->t == SERIES) {
	    /* using RHS series */
	    int i, s = p->dset->t1;

	    for (i=0; i<m->rows; i++) {
		x = p->ret->v.xvec[s++];
		m->val[i] = na(x)? M_NA : x;
	    }
	    *prechecked = 1;
	} else {
	    /* using RHS matrix: just copy data across */
	    p->err = gretl_matrix_copy_data(m, p->ret->v.m);
	}
    } else {
	/* Dimensions diff: replace the LHS matrix */
#if EDEBUG
	fprintf(stderr, "assign_to_matrix: replacing\n");
#endif
	m = retrieve_matrix_result(p, prechecked);
	if (!p->err) {
	    p->err = gen_replace_lhs(p, GRETL_TYPE_MATRIX, m);
	}
    }

    return m;
}

/* Assigning to an existing (whole) LHS matrix, but using '+='
   or some such modified/inflected assignment. Note that
   save_generated_var() is the only caller.
*/

static gretl_matrix *assign_to_matrix_mod (gretl_matrix *m1,
					   parser *p,
					   int *prechecked)
{
    gretl_matrix *m2 = NULL;

    if (m1 == NULL) {
	p->err = E_DATA;
    }

    if (!p->err) {
	if (p->op == B_DOTASN) {
	    p->err = dot_assign_to_matrix(m1, p, prechecked);
	    m2 = m1; /* no change in matrix pointer */
	} else {
	    gretl_matrix *tmp = retrieve_matrix_result(p, prechecked);

	    if (tmp != NULL) {
		p->err = real_matrix_calc(m1, tmp, p->op, &m2);
		gretl_matrix_free(tmp);
	    }
	}
    }

    return m2;
}

static void do_array_append (parser *p)
{
    gretl_array *A = NULL;
    GretlType atype;
    NODE *rhs = p->ret;
    void *ptr = NULL;

    A = gen_get_lhs_var(p, GRETL_TYPE_ARRAY);

    if (A == NULL) {
	p->err = E_DATA;
	return;
    }

    atype = gretl_array_get_content_type(A);

    if (atype == GRETL_TYPE_STRING && rhs->t == STR) {
	ptr = rhs->v.str;
    } else if (atype == GRETL_TYPE_MATRIX && rhs->t == MAT) {
	ptr = rhs->v.m;
    } else if (atype == GRETL_TYPE_BUNDLE && rhs->t == BUNDLE) {
	ptr = rhs->v.b;
    } else if (atype == GRETL_TYPE_LIST && rhs->t == LIST) {
	ptr = rhs->v.ivec;
    } else if (rhs->t == ARRAY) {
	/* special: not actually appending an _element_;
	   stick rhs array onto end of lhs array
	*/
	p->err = gretl_array_append_array(A, rhs->v.a);
    } else {
	p->err = E_TYPES;
    }

    if (!p->err && ptr != NULL) {
	int copy = !is_tmp_node(rhs);

	p->err = gretl_array_append_object(A, ptr, copy);
	if (!copy && !p->err) {
	    rhs->v.ptr = NULL;
	}
    }
}

static int create_or_edit_string (parser *p)
{
    const char *src = NULL;
    const char *orig = NULL;
    char *newstr = NULL;
    user_var *uvar;

    if (p->ret->t == NUM) {
	/* OK only in case of "+=" */
	if (p->op != B_ADD) {
	    p->err = E_TYPES;
	    return p->err;
	}
    } else if (p->ret->t == EMPTY) {
	src = "";
    } else {
	src = p->ret->v.str;
    }

#if EDEBUG
    fprintf(stderr, "edit_string: src='%s'\n", src);
#endif

    uvar = p->lh.uv;

    if (uvar != NULL) {
	orig = uvar->ptr;
    } else if (p->op != B_ASN) {
	/* without an existing LHS string we can only assign */
	p->err = E_DATA;
	return p->err;
    }

    if (p->ret->t == NUM) {
	/* taking an offset into an existing string */
	int len = strlen(orig);
	int adj = p->ret->v.xval;

	if (adj < 0) {
	    p->err = E_DATA;
	} else if (adj == 0) {
	    ; /* no-op */
	} else {
	    src = (adj < len)? (orig + adj) : "";
	    newstr = gretl_strdup(src);
	    if (newstr == NULL) {
		p->err = E_ALLOC;
	    } else {
		gen_replace_lhs(p, GRETL_TYPE_STRING, newstr);
	    }
	}
    } else if (src == NULL) {
	; /* no-op -- e.g. argname() didn't get anything */
    } else if (p->op == B_ASN) {
	/* simple assignment */
	newstr = gretl_strdup(src);
	if (newstr == NULL) {
	    p->err = E_ALLOC;
	} else if (uvar == NULL) {
	    gen_add_uvar(p, GRETL_TYPE_STRING, newstr);
	} else {
	    gen_replace_lhs(p, GRETL_TYPE_STRING, newstr);
	}
    } else if (p->op == B_HCAT || p->op == B_ADD) {
	/* string concatenation */
	if (*src == '\0') {
	    ; /* no-op */
	} else {
	    newstr = malloc(strlen(orig) + strlen(src) + 1);
	    if (newstr == NULL) {
		p->err = E_ALLOC;
	    } else {
		strcpy(newstr, orig);
		strcat(newstr, src);
		gen_replace_lhs(p, GRETL_TYPE_STRING, newstr);
	    }
	}
    }

    return p->err;
}

static int create_or_edit_list (parser *p)
{
    int *list = node_get_list(p->ret, p); /* note: copied */

#if EDEBUG
    printlist(list, "incoming list in edit_list()");
#endif

    if (!p->err) {
	if (p->lh.t != LIST) {
	    /* no pre-existing LHS list: must be simple assignment */
	    p->err = remember_list(list, p->lh.name, NULL);
	} else if (p->op == B_ASN || p->op == B_ADD || p->op == B_SUB) {
	    /* replace, append or subtract list members */
	    p->err = gen_edit_list(p, list, p->op);
	} else {
	    p->err = E_TYPES;
	}
    }

    if (!p->err) {
	set_dataset_is_changed();
    }

    free(list);

    return p->err;
}

#define arrayable_type(t) (t == MAT || t == STR || t == BUNDLE || \
			   t == LIST)

#define ok_return_type(t) (t == NUM || t == SERIES || t == MAT ||	\
			   t == LIST || t == DUM || t == EMPTY ||	\
                           t == STR || t == BUNDLE || t == ARRAY ||	\
			   t == U_ADDR || t == DBUNDLE)

/* Note: we're doing this only in relation to "primary" types
   (excluding bundle members, array elements, matrix sub-
   specs).
*/

static int gen_check_return_type (parser *p)
{
    NODE *r = p->ret;
    int err = 0;

    if (r == NULL) {
	fprintf(stderr, "gen_check_return_type: p->ret = NULL!\n");
	return E_DATA;
    }

#if EDEBUG
    fprintf(stderr, "gen_check_return_type: targ=%s; ret at %p, type %s\n",
	    getsymb(p->targ), (void *) r, getsymb(r->t));
#endif

    if (!ok_return_type(r->t)) {
	return E_TYPES;
    }

    if (r->t == SERIES && r->v.xvec == NULL) {
	fprintf(stderr, "got SERIES return with xvec = NULL!\n");
	return E_DATA;
    }

    if (p->targ == NUM) {
	if (r->t == NUM || scalar_matrix_node(r)) {
	    ; /* scalar or 1 x 1 matrix: OK */
	} else if (r->t == MAT && (p->flags & P_NODECL)) {
	    ; /* morphing to matrix may be OK */
	} else {
	    err = E_TYPES;
	}
    } else if (p->targ == SERIES) {
	/* result must be scalar, series, or conformable matrix */
	if (r->t == NUM || r->t == SERIES) {
	    ; /* OK */
	} else if (r->t == MAT) {
	    if (!series_compatible(r->v.m, p)) {
		err = E_TYPES;
	    }
	} else {
	    err = E_TYPES;
	}
    } else if (p->targ == MAT) {
	; /* no-op: handled later */
    } else if (p->targ == LIST) {
	if (r->t != EMPTY && !ok_list_node(r)) {
	    err = E_TYPES;
	}
    } else if (p->targ == STR) {
	if (r->t != EMPTY && r->t != STR && r->t != NUM) {
	    err = E_TYPES;
	}
    } else if (p->targ == BUNDLE) {
	if (p->op == B_ASN) {
	    /* plain assignment: bundle or null */
	    if (r->t != BUNDLE && r->t != DBUNDLE && r->t != EMPTY) {
		err = E_TYPES;
	    }
	} else {
	    /* the only other assignment possibility is "+=",
	       in which case we'll only accept a bundle
	    */
	    if (r->t != BUNDLE) {
		err = E_TYPES;
	    }
	}
    } else if (p->targ == ARRAY) {
	if (p->op == B_ASN) {
	    /* plain assignment: array or null */
	    if (r->t != ARRAY && r->t != EMPTY) {
		err = E_TYPES;
	    }
	} else {
	    /* arrays: the only other assignment possibility is "+=",
	       in which case we'll only accept an array or an
	       object which matches the content type on the left
	       (but the matching check is deferred)
	    */
	    if (r->t != ARRAY && !arrayable_type(r->t)) {
		err = E_TYPES;
	    }
	}
    }

    if (err == E_TYPES) {
	maybe_do_type_errmsg(p->lh.name, p->lh.t);
    }

#if EDEBUG
    fprintf(stderr, "gen_check_return_type: returning with p->err = %d\n",
	    err);
#endif

    return err;
}

/* allocate storage if saving a series to the dataset:
   lh.vnum <= 0 means that the LHS series does not already
   exist
*/

static int gen_allocate_storage (parser *p)
{
    if (p->lh.vnum <= 0) {
	if (p->dset == NULL || p->dset->Z == NULL) {
	    p->err = E_DATA;
	} else {
	    p->err = dataset_add_NA_series(p->dset, 1);
	    if (!p->err) {
		p->lh.vnum = p->dset->v - 1;
	    }
	}
    }

    return p->err;
}

#if SERIES_ENSURE_FINITE

static void series_make_finite (double *x, int n)
{
    int i;

    /* This may be questionable (converting NaNs and infinities to
       NA), but is necessary for backward compatibility (e.g. when
       taking logs of a series that contains zeros).
    */

    for (i=0; i<n; i++) {
	if (!isfinite(x[i])) {
	    x[i] = NADBL;
	    set_gretl_warning(W_GENMISS);
	}
    }
}

#endif

static void align_matrix_to_series (double *y, const gretl_matrix *m,
				    parser *p)
{
    const gretl_matrix *mask = get_matrix_mask();
    int t, s = 0;

    if (mask == NULL || mask->rows != p->dset->n) {
	p->err = E_DATA;
	return;
    }

    for (t=0; t<p->dset->n; t++) {
	if (mask->val[t] != 0.0) {
	    if (t >= p->dset->t1 && t <= p->dset->t2) {
		y[t] = xy_calc(y[t], m->val[s], p->op, SERIES, p);
	    }
	    s++;
	}
    }
}

static int assign_null_to_bundle (parser *p)
{
    gretl_bundle *b;
    int err = 0;

    if (p->lh.t == BUNDLE) {
	b = gen_get_lhs_var(p, GRETL_TYPE_BUNDLE);
	gretl_bundle_void_content(b);
    } else {
	b = gretl_bundle_new();
	if (b == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gen_add_uvar(p, GRETL_TYPE_BUNDLE, b);
	}
    }

    return err;
}

static int assign_null_to_array (parser *p)
{
    gretl_array *a;
    int err = 0;

    if (p->lh.t == ARRAY) {
	a = gen_get_lhs_var(p, GRETL_TYPE_ARRAY);
	gretl_array_void_content(a);
    } else {
	a = gretl_array_new(p->lh.gtype, 0, &err);
	if (!err) {
	    err = gen_add_uvar(p, p->lh.gtype, a);
	}
    }

    return err;
}

/* apply postfix '++' or '--' to LHS scalar, or '++' to
   LHS string (only) */

static int do_incr_decr (parser *p)
{
    if (p->lh.uv != NULL && p->lh.uv->type == GRETL_TYPE_DOUBLE) {
	double x = *(double *) p->lh.uv->ptr;

	if (!na(x)) {
	    x += (p->op == INC)? 1 : -1;
	    *(double *) p->lh.uv->ptr = x;
	}
    } else if (p->lh.uv != NULL && p->lh.uv->type == GRETL_TYPE_STRING) {
	if (p->op == DEC) {
	    p->err = E_TYPES;
	} else {
	    char *s = p->lh.uv->ptr;

	    if (*s != '\0') {
		char *smod = gretl_strdup(s + 1);

		gen_replace_lhs(p, GRETL_TYPE_STRING, smod);
	    }
	}
    } else {
	p->err = E_TYPES;
    }

    return p->err;
}

static int save_generated_var (parser *p, PRN *prn)
{
    NODE *r = p->ret;
    double **Z = NULL;
    double x;
    int no_decl = 0;
    int t, v = 0;

#if EDEBUG
    fprintf(stderr, "save (%s): '%s'\n  callcount=%d\n"
	    "lh.t=%s, targ=%s, no_decl=%d, r->t=%s\n",
	    p->lhtree != NULL ? "compound" : "unitary",
	    p->lh.name, p->callcount, getsymb(p->lh.t),
	    getsymb(p->targ), (p->flags & P_NODECL)? 1 : 0,
	    (r == NULL)? "none" : getsymb(r->t));
#endif    

    if (p->lhtree != NULL) {
	/* handle compound target first */
	int compound_t;
	
	p->lhtree->flags |= LHT_NODE;
	p->flags |= P_START;
#if LHDEBUG	
	fprintf(stderr, "*** eval lhtree -> lhres ***\n");
#endif
	p->lhres = eval(p->lhtree, p);
#if LHDEBUG	
	if (p->lhres != NULL) {
	    print_tree(p->lhres, p, 0);
	    fprintf(stderr, "*** lhtree post-eval ***\n");
	    print_tree(p->lhtree, p, 0);
	}
#endif
	if (p->err) {
	    return p->err;
	}
	compound_t = p->lhres->t;
#if LHDEBUG
	fprintf(stderr, "save_generated_var: type = %s\n",
		getsymb(compound_t));
#endif
	if (compound_t == BMEMB) {
	    p->err = set_bundle_value(p->lhres, r, p);
	} else if (compound_t == MSL) {
	    p->err = set_matrix_value(p->lhres, r, p);
	} else if (compound_t == OBS) {
	    p->err = set_series_obs_value(p->lhres, r, p);
	} else if (compound_t == OSL) {
	    NODE *lh1 = p->lhres->v.b2.l;

	    if (lh1->t == ARRAY) {
		p->err = set_array_value(p->lhres, r, p);
	    } else if (lh1->t == LIST) {
		p->err = set_list_value(p->lhres, r, p);
	    } else if (lh1->t == STR) {
		p->err = set_string_value(p->lhres, r, p);
	    } else if (lh1->t == MAT) {
		/* is this ever reached?? */
		p->err = set_matrix_value(p->lhres, r, p);
	    }
	} else {
	    gretl_errmsg_set(_("Invalid left-hand side expression"));
	    p->err = E_TYPES;
	}
	return p->err; /* done */
    }

    if (p->op == INC || p->op == DEC) {
	return do_incr_decr(p);
    }    

    if (p->callcount < 2) {
	/* first exec: test for type mismatch errors */
	p->err = gen_check_return_type(p);
	if (p->err) {
	    return p->err;
	}
    }

#if ONE_BY_ONE_CAST
    if (p->targ == UNK) {
	if (scalar_matrix_node(r)) {
	    /* "cast" 1 x 1 matrix to scalar */
	    no_decl = 1;
	    p->targ = NUM;
	    p->flags |= P_NODECL;
	} else {
	    p->targ = r->t;
	}
    } else if (p->targ == NUM && r->t == MAT && (p->flags & P_NODECL)) {
	/* We're looking at a @targ that was previously
	   set to NUM by the "auto-cast" mechanism: allow
	   it to morph to matrix if need be.
	*/
	if (scalar_matrix_node(r)) {
	    ; /* not a problem */
	} else if (p->lh.t == 0) {
	    /* no pre-existing scalar var */
	    p->targ = MAT;
	} else if (p->lh.t == NUM) {
	    /* type-convert existing scalar */
	    p->err = gretl_scalar_convert_to_matrix(p->lh.uv);
	    if (!p->err) {
		p->targ = MAT;
	    }
	}
    }
#else
    if (p->targ == UNK) {
	p->targ = r->t;
    }
#endif

#if EDEBUG
    fprintf(stderr, "after preliminaries: targ=%s, op='%s'\n",
	    getsymb(p->targ), getsymb(p->op));
#endif

    /* allocate dataset storage, if needed */
    if (p->targ == SERIES) {
	gen_allocate_storage(p);
	if (p->err) {
	    return p->err;
	}
    }

    if (p->dset != NULL && p->dset->Z != NULL) {
	Z = p->dset->Z;
	v = p->lh.vnum;
    }

    /* put the generated data into place */

    if (p->targ == NUM) {
	if (p->lh.t == NUM) {
	    /* modifying existing scalar */
	    x = gretl_scalar_get_value(p->lh.name, NULL);
	    if (r->t == NUM) {
		x = xy_calc(x, r->v.xval, p->op, NUM, p);
	    } else if (scalar_matrix_node(r)) {
		x = xy_calc(x, r->v.m->val[0], p->op, NUM, p);
	    } else {
		p->err = E_TYPES;
	    }
	    if (!p->err) {
		p->err = gen_replace_scalar(p, x);
	    }
	} else {
	    /* a new scalar */
	    if (r->t == NUM) {
		x = r->v.xval;
	    } else if (scalar_matrix_node(r)) {
		x = r->v.m->val[0];
	    } else {
		p->err = E_TYPES;
	    }
#if SCALARS_ENSURE_FINITE
	    if (!p->err && !isfinite(x)) {
		x = NADBL;
	    }
#endif
	    if (!p->err) {
		if (no_decl) {
		    p->err = gretl_scalar_add_mutable(p->lh.name, x);
		} else {
		    p->err = gretl_scalar_add(p->lh.name, x);
		}
	    }
	}
    } else if (p->targ == SERIES) {
	/* writing a series */
	if (r->t == NUM) {
	    for (t=p->dset->t1; t<=p->dset->t2; t++) {
		Z[v][t] = xy_calc(Z[v][t], r->v.xval, p->op, SERIES, p);
	    }
	} else if (r->t == SERIES) {
	    const double *x = r->v.xvec;
	    int t1 = p->dset->t1;

	    if (autoreg(p) && p->op == B_ASN) {
		while (xna(x[t1]) && t1 <= p->dset->t2) {
		    t1++;
		}
	    }
	    if (p->op == B_ASN) {
		/* avoid multiple calls to xy_calc */
		if (Z[v] != x) {
		    size_t sz = (p->dset->t2 - t1 + 1) * sizeof *x;

		    memcpy(Z[v] + t1, x + t1, sz);
		}
	    } else {
		for (t=t1; t<=p->dset->t2; t++) {
		    Z[v][t] = xy_calc(Z[v][t], x[t], p->op, SERIES, p);
		}
	    }
	} else if (r->t == MAT) {
	    const gretl_matrix *m = r->v.m;
	    int k = gretl_vector_get_length(m);
	    int mt1 = gretl_matrix_get_t1(m);
	    int s;

	    if (p->flags & P_MMASK) {
		/* result needs special alignment */
		align_matrix_to_series(Z[v], m, p);
	    } else if (k == 1) {
		/* result is effectively a scalar */
		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    Z[v][t] = xy_calc(Z[v][t], m->val[0], p->op, SERIES, p);
		}
	    } else if (k == p->dset->n) {
		/* treat result as full-length series */
		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    Z[v][t] = xy_calc(Z[v][t], m->val[t], p->op, SERIES, p);
		}
	    } else if (k == sample_size(p->dset)) {
		/* treat as series of current sample length
		   2016-03-03: removed clause "&& mt1 == 0"
		*/
		for (t=p->dset->t1, s=0; t<=p->dset->t2; t++, s++) {
		    Z[v][t] = xy_calc(Z[v][t], m->val[s], p->op, SERIES, p);
		}
	    } else if (mt1 > 0) {
		/* align using matrix "t1" value */
		for (t=mt1; t<mt1 + k && t<=p->dset->t2; t++) {
		    if (t >= p->dset->t1) {
			Z[v][t] = xy_calc(Z[v][t], m->val[t - mt1], p->op,
					  SERIES, p);
		    }
		}
	    }
	}
#if SERIES_ENSURE_FINITE
	if (!p->err) {
	    series_make_finite(Z[v], p->dset->n);
	}
#endif
	strcpy(p->dset->varname[v], p->lh.name);
#if EDEBUG
	fprintf(stderr, "var %d: gave generated series the name '%s'\n",
		v, p->lh.name);
	fprintf(stderr, " value[1] = %g\n", p->dset->Z[v][1]);
#endif
	if (!p->err) {
	    /* (probably) changed or added a series */
	    set_dataset_is_changed();
	}
    } else if (p->targ == MAT) {
	/* we're writing a matrix */
	gretl_matrix *m = NULL;
	int prechecked = 0;

	if (p->lh.uv == NULL) {
	    /* there's no pre-existing left-hand side matrix */
	    m = retrieve_matrix_result(p, &prechecked);
	    if (!p->err) {
		p->err = gen_add_uvar(p, GRETL_TYPE_MATRIX, m);
	    }	    
	} else if (p->op == B_ASN) {
	    /* uninflected assignment to an existing matrix */
	    m = assign_to_matrix(p, &prechecked);
	} else {
	    /* inflected assignment to entire existing matrix */
	    gretl_matrix *m1 = gen_get_lhs_var(p, GRETL_TYPE_MATRIX);
	    
	    m = assign_to_matrix_mod(m1, p, &prechecked);
	    if (!p->err) {
		p->err = gen_replace_lhs(p, GRETL_TYPE_MATRIX, m);
	    }
	}
	/* note: for use by genr_get_output_matrix() */
	p->lh.mret = m;
#if MATRIX_NA_CHECK
	if (!p->err && !prechecked && m != NULL &&
	    gretl_matrix_xna_check(m)) {
	    set_gretl_warning(W_GENNAN);
	    prechecked = 1;
	}
#endif
    } else if (p->targ == LIST) {
	create_or_edit_list(p);
    } else if (p->targ == STR) {
	create_or_edit_string(p);
    } else if (p->targ == BUNDLE) {
	if (r->t == EMPTY) {
	    /* as in "bundle b = null" */
	    p->err = assign_null_to_bundle(p);
	} else {
	    /* full assignment of RHS bundle */
	    gretl_bundle *b;

	    if (r->t == DBUNDLE) {
		b = bvar_get_bundle(r, p);
	    } else if (is_tmp_node(r) || (p->flags & P_UFRET)) {
		/* grabbing r->v.b is OK */
		b = r->v.b;
	    } else {
		/* we need to make a copy */
		b = gretl_bundle_copy(r->v.b, &p->err);
	    }

	    if (!p->err) {
		p->err = gen_add_or_replace(p, GRETL_TYPE_BUNDLE, b);
		if (!p->err && r->t != DBUNDLE && b == r->v.b) {
		    /* avoid destroying the assigned bundle */
		    r->v.b = NULL;
		}
	    }
	}
    } else if (p->targ == ARRAY) {
	if (p->op != B_ASN) {
	    do_array_append(p);
	} else if (r->t == EMPTY) {
	    /* as in, e.g., "strings A = null" */
	    p->err = assign_null_to_array(p);
	} else {
	    /* full assignment of RHS array */
	    GretlType atype = p->lh.gtype > 0 ? p->lh.gtype :
		gretl_array_get_type(r->v.a);
	    gretl_array *a;

	    if (is_tmp_node(r) || (p->flags & P_UFRET)) {
		/* grabbing r->v.a is OK */
		a = r->v.a;
	    } else {
		/* we need to make a copy */
		a = gretl_array_copy(r->v.a, &p->err);
	    }

	    if (!p->err) {
		p->err = gen_add_or_replace(p, atype, a);
		if (!p->err && a == r->v.a) {
		    /* avoid destroying the assigned array */
		    r->v.a = NULL;
		}
	    }
	}
    }

#if EDEBUG
    fprintf(stderr, "save_generated_var: returning p->err = %d\n",
	    p->err);
#endif

    return p->err;
}

static void maybe_update_lhs_uvar (parser *p, GretlType *type)
{
    if (p->targ == SERIES) {
	/* targetting a series */
	int v = p->lh.vnum;

	if (v <= 0 || v >= p->dset->v) {
	    p->lh.vnum = current_series_index(p->dset, p->lh.name);
	}
	if (p->lh.vnum < 0) {
	    p->lh.vnum = 0;
	}
	return;
    }

    if (p->lh.uv == NULL) {
	p->lh.uv = get_user_var_by_name(p->lh.name);
    }

    if (p->lh.uv != NULL) {
	*type = p->lh.uv->type;
    }

    switch (*type) {
    case GRETL_TYPE_DOUBLE:
	p->lh.t = NUM;
	break;
    case GRETL_TYPE_MATRIX:
	p->lh.t = MAT;
	if (p->targ == NUM) {
	    p->targ = MAT;
	}
	break;
    case GRETL_TYPE_LIST:
	p->lh.t = LIST;
	break;
    case GRETL_TYPE_STRING:
	p->lh.t = STR;
	break;
    case GRETL_TYPE_BUNDLE:
	p->lh.t = BUNDLE;
	break;
    case GRETL_TYPE_ARRAY:
	p->lh.t = ARRAY;
	break;
    default:
	p->lh.t = 0;
	break;
    }
}

static void parser_reinit (parser *p, DATASET *dset, PRN *prn)
{
    /* flags that should be reinstated if they were
       set at compile time, or in previous execution
    */
    int saveflags[] = {
	P_NATEST, P_AUTOREG, P_SLAVE,
	P_DISCARD, P_NODECL, P_LISTDEF,
	0
    };
    int i, prevflags = p->flags;
    GretlType lhtype = 0;
    int dset_n;

    if (p->callcount > 1) {
	/* the flags should basically have stabilized by now */
	p->flags |= P_START;
	p->flags &= ~P_DELTAN;
    } else {
	p->flags = (P_START | P_PRIV | P_EXEC);

	for (i=0; saveflags[i] > 0; i++) {
	    if (prevflags & saveflags[i]) {
		p->flags |= saveflags[i];
	    }
	}
    }

    p->dset = dset;
    p->prn = prn;

    p->obs = 0;
    p->sym = 0;
    p->ch = 0;
    p->xval = 0.0;
    p->idnum = 0;
    p->idstr = NULL;
    p->data = NULL;

    p->ret = NULL;
    p->lhres = NULL;
    p->err = 0;

#if EDEBUG
    fprintf(stderr, "parser_reinit: targ=%s, lhname='%s', op='%s', "
	    "callcount=%d, compiled=%d\n",
	    getsymb(p->targ), p->lh.name, getsymb(p->op),
	    p->callcount, compiled(p));
#endif

    if (*p->lh.name != '\0') {
	maybe_update_lhs_uvar(p, &lhtype);
    }

    /* allow for change in length of dataset */
    dset_n = dset != NULL ? dset->n : 0;
    if (dset_n != p->dset_n) {
	p->dset_n = dset_n;
	p->flags |= P_DELTAN;
    }
}

static void parser_init (parser *p, const char *str,
			 DATASET *dset, PRN *prn,
			 int flags, int targtype)
{
    p->point = p->rhs = p->input = str;
    p->dset = dset;
    p->dset_n = dset != NULL ? dset->n : 0;
    p->prn = prn;
    p->flags = flags | P_START;
    p->targ = targtype;
    p->op = 0;

    p->lhtree = NULL;
    p->lhres = NULL;
    p->tree = NULL;
    p->ret = NULL;

    /* left-hand side info */
    p->lh.t = 0;
    p->lh.name[0] = '\0';
    p->lh.label = NULL;
    p->lh.vnum = 0;
    p->lh.uv = NULL;
    p->lh.expr = NULL;
    p->lh.gtype = 0;
    p->lh.mret = NULL;

    /* auxiliary apparatus */
    p->aux = NULL;

    p->callcount = 0;
    p->obs = 0;
    p->sym = 0;
    p->upsym = 0;
    p->ch = 0;
    p->xval = 0.0;
    p->idnum = 0;
    p->idstr = NULL;
    p->err = 0;

    if (p->input == NULL) {
	p->err = E_DATA;
	return;
    }

    if (p->flags & P_VOID) {
        p->flags |= P_DISCARD;
    } else if (p->targ == UNK || !(p->flags & P_ANON)) {
	gen_preprocess(p, flags);
    } else if (p->targ == LIST) {
	p->flags |= P_LISTDEF;
    }

    if (!p->err) {
	p->ch = parser_getc(p);
    }
}

/* called from genmain.c (only!) */

void gen_save_or_print (parser *p, PRN *prn)
{
    if (p->flags & P_DISCARD) {
	if (p->ret == NULL) {
	    return;
	} else if (p->ret->t == MAT) {
	    gretl_matrix_print_to_prn(p->ret->v.m, p->lh.name, p->prn);
	} else if (p->ret->t == LIST) {
	    if (p->lh.name[0] != '\0') {
		gretl_list_print(get_list_by_name(p->lh.name),
				 p->dset, p->prn);
	    } else {
		gretl_list_print(p->ret->v.ivec, p->dset, p->prn);
	    }
	} else if (p->ret->t == STR) {
	    if (p->lh.name[0] != '\0') {
		pprintf(p->prn, "%s\n", gen_get_lhs_var(p, GRETL_TYPE_STRING));
	    } else {
		pprintf(p->prn, "%s\n", p->ret->v.str);
	    }
	} else if (p->ret->t == BUNDLE) {
	    gretl_bundle_print(p->ret->v.b, prn);
	} else if (p->ret->t == ARRAY) {
	    gretl_array_print(p->ret->v.a, prn);
	} else {
	    /* scalar, series */
	    printnode(p->ret, p, 1);
	    pputc(p->prn, '\n');
	}
    } else if (p->flags & P_DECL) {
	do_declaration(p);
    } else {
	save_generated_var(p, prn);
    }
}

void gen_cleanup (parser *p)
{
    int save = (p->flags & (P_COMPILE | P_EXEC));

#if EDEBUG
    fprintf(stderr, "gen cleanup on %p: save = %d\n",
	    p, save ? 1 : 0);
#endif

    if (p->lh.label != NULL) {
	free(p->lh.label);
	p->lh.label = NULL;
    }

    if (p->err && (p->flags & P_COMPILE)) {
	save = 0;
    }

    if (!save) {
	if (p->lhtree != NULL) {
	    if (p->lhtree != p->lhres) {
		/* we have to scrub the LHT_NODE flag on p->lhtree,
		   or else its children will not get freed and we'll
		   leak memory
		*/
		p->lhtree->flags &= ~LHT_NODE;
		rndebug(("freeing p->lhtree %p\n", (void *) p->lhtree));
		free_tree(p->lhtree, p, FR_LHTREE);
	    }
	    if (p->lhres != NULL) {
		rndebug(("freeing p->lhres %p\n", (void *) p->lhres));
		free_tree(p->lhres, p, FR_LHRES);
	    }
	}

	if (p->tree != p->ret) {
	    rndebug(("freeing p->tree %p\n", (void *) p->tree));
	    free_tree(p->tree, p, FR_TREE);
	}

	rndebug(("freeing p->ret %p\n", (void *) p->ret));
	free_tree(p->ret, p, FR_RET);

	free(p->lh.expr);
    }
}

#define LS_DEBUG 0

static void real_reset_uvars (parser *p)
{
    if (p->err) {
	return;
    }

#if LS_DEBUG
    fprintf(stderr, "\nreal_reset_uvars (%s '%s') *\n",
	    getsymb(p->targ), p->lh.name);
#endif

    clear_uvnodes(p->tree);

    if (p->lhtree != NULL) {
	clear_uvnodes(p->lhtree);
    }

    p->lh.uv = NULL;
    p->lh.vnum = 0;
}

void genr_reset_uvars (parser *p)
{
    real_reset_uvars(p);
}

static void maybe_set_return_flags (parser *p)
{
    NODE *t = p->tree;

    if (t != NULL && t->t == UFUN) {
	p->flags |= P_UFRET;
    }
}

static int decl_check (parser *p, int flags)
{
    if (flags & P_COMPILE) {
	p->err = E_PARSE;
	gretl_errmsg_sprintf("%s:\n> '%s'",
			     _("Bare declarations are not allowed here"),
			     p->input);
    }

    return p->err;
}

static void autoreg_error (parser *p, int t)
{
    fprintf(stderr, "*** autoreg error at obs t = %d (t1 = %d):\n",
	    t, p->dset->t1);

    if (p->ret != NULL && p->ret->t != SERIES) {
	fprintf(stderr, " ret type != SERIES (=%d), p->err = %d\n",
		p->ret->t, p->err);
    } else if (p->ret == NULL) {
	fprintf(stderr, " ret = NULL, p->err = %d\n", p->err);
    }

    fprintf(stderr, " input = '%s'\n", p->input);

    if (!p->err) {
	p->err = E_DATA;
    }
}

int realgen (const char *s, parser *p, DATASET *dset, PRN *prn,
	     int flags, int targtype)
{
#if LHDEBUG || EDEBUG || AUX_NODES_DEBUG
    fprintf(stderr, "\n*** realgen: task = %s\n", (flags & P_COMPILE)?
	    "compile" : (flags & P_EXEC)? "exec" : "normal");
    if (s != NULL) {
	fprintf(stderr, "targ=%d (%s), input='%s'\n", targtype,
		(targtype < PUNCT_MAX)? gretl_type_get_name(targtype) :
		getsymb(targtype), s);
    }
#endif

    if (flags & P_EXEC) {
#if EDEBUG
	fprintf(stderr, "*** printing p->tree (before reinit)\n");
	print_tree(p->tree, p, 0);
#endif
	parser_reinit(p, dset, prn);
	if (p->err) {
	    fprintf(stderr, "error in parser_reinit\n");
	    return p->err;
	} else if (p->op == INC || p->op == DEC) {
	    /* more or less a no-op: the work is done by
	       save_generated_var()
	    */
	    return p->err;
	} else {
	    goto starteval;
	}
    } else {
	parser_init(p, s, dset, prn, flags, targtype);
	if (p->err) {
	    if (gretl_function_depth() == 0) {
		errmsg(p->err, prn);
	    }
	    return p->err;
	}
    }

#if EDEBUG
    fprintf(stderr, "after parser (re-)init, p->err = %d (decl? %s)\n",
	    p->err, (p->flags & P_DECL)? "yes" : "no");
#endif

    if (p->flags & P_DECL) {
	/* check validity of declaration(s) */
	decl_check(p, flags);
	return p->err;
    }

    if (p->op == INC || p->op == DEC) {
	/* implemented via save_generated_var() */
	return p->err;
    }

    /* fire up the lexer */
    lex(p);
    if (p->err) {
#if EDEBUG
	fprintf(stderr, "realgen %p ('%s'): exiting on lex() error %d\n",
		(void *) p, s, p->err);
#endif
	return p->err;
    }

    /* build the syntax tree */
    p->tree = expr(p);
    if (p->err) {
	return p->err;
    }

#if EDEBUG
    fprintf(stderr, "realgen: p->tree at %p, type %d (%s)\n", (void *) p->tree,
	    p->tree->t, getsymb(p->tree->t));
    if (p->ch == '\0') {
	fprintf(stderr, " p->ch = NUL, p->sym = %d\n", p->sym);
    } else {
	fprintf(stderr, " p->ch = '%c', p->sym = %d\n", p->ch, p->sym);
    }
#endif

    if (p->sym != EOT || p->ch != 0) {
	int c = p->ch;

	if (c == ' ') {
	    c = 0;
	} else if (c != 0) {
	    parser_ungetc(p);
	    c = p->ch;
	}
	context_error(c, p, "realgen");
	return p->err;
    }

    if (flags & P_NOEXEC) {
	/* we're done at this point */
	return p->err;
    }

    if (!p->err) {
	/* set P_UFRET here if relevant */
	maybe_set_return_flags(p);
    }

 starteval:

#if EDEBUG
    if (flags & P_EXEC) {
	fprintf(stderr, "*** printing p->tree (about to start eval)\n");
	print_tree(p->tree, p, 0);
    }
#endif

    if (autoreg(p)) {
	/* e.g. y = b*y(-1) : evaluate dynamically */
	int t;

	for (t=p->dset->t1; t<=p->dset->t2 && !p->err; t++) {
	    const double *x;

	    /* initialize for this observation */
	    p->obs = t;
#if EDEBUG
	    fprintf(stderr, "\n*** autoreg: p->obs = %d\n", p->obs);
#endif
	    p->ret = eval(p->tree, p);
	    if (p->ret != NULL && p->ret->t == SERIES) {
		x = p->ret->v.xvec;
		if (!na(x[t])) {
#if EDEBUG
		    fprintf(stderr, "writing xvec[%d] = %g into Z[%d][%d]\n",
			    t, x[t], p->lh.vnum, t);
#endif
		    p->dset->Z[p->lh.vnum][t] = x[t];
		}
	    } else {
		autoreg_error(p, t);
	    }
	    if (t == p->dset->t1) {
		p->flags &= ~P_START;
	    }
	}
    } else {
	/* standard non-dynamic evaluation */
	p->ret = eval(p->tree, p);
    }

    if (p->flags & P_EXEC) {
	p->callcount += 1;
    }

#if EDEBUG
    fprintf(stderr, "realgen: post-eval, err = %d\n", p->err);
# if EDEBUG > 1
    printnode(p->ret, p, 0);
    pputc(prn, '\n');
# endif
#endif

    gen_check_errvals(p);

    return p->err;
}
