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
#include "version.h"
#include "monte_carlo.h"
#include "gretl_func.h"
#include "libset.h"
#include "usermat.h"
#include "gretl_xml.h"
#include "cmd_private.h"
#include "gretl_string_table.h"
#include "gretl_array.h"
#include "gretl_typemap.h"
#include "gretl_zip.h"
#include "uservar.h"
#include "flow_control.h"
#include "kalman.h"
#include "system.h"
#include "genr_optim.h"

#include <errno.h>
#include <glib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define LSDEBUG 0

#define FNPARSE_DEBUG 0 /* debug parsing of function code */
#define EXEC_DEBUG 0    /* debugging of function execution */
#define UDEBUG 0        /* debug handling of args */
#define PKG_DEBUG 0     /* debug handling of function packages */
#define FN_DEBUG 0      /* miscellaneous debugging */
#define DDEBUG 0        /* debug the debugger */

#define INT_USE_XLIST (-999)

typedef struct fn_param_ fn_param;
typedef struct fn_arg_ fn_arg;
typedef struct fn_line_ fn_line;
typedef struct obsinfo_ obsinfo;
typedef struct fncall_ fncall;

/* structure representing a parameter of a user-defined function */

struct fn_param_ {
    char *name;     /* the name of the parameter */
    char type;      /* its type */
    char *descrip;  /* its description */
    char **labels;  /* value labels, if applicable */
    int nlabels;    /* number of value labels */
    char flags;     /* additional information (e.g. "const" flag) */
    double deflt;   /* default value */
    double min;     /* minimum value (scalar parameters only) */
    double max;     /* maximum value (scalar parameters only) */
    double step;    /* step increment (scalars only) */
};

/* structure representing a line of a user-defined function */

struct fn_line_ {
    int idx;        /* 1-based line index (allowing for blanks) */
    char *s;        /* text of command line */
    LOOPSET *loop;  /* attached "compiled" loop */
    int next_idx;   /* line index to skip to after loop */
    int ignore;     /* flag for comment lines */
};

#define UNSET_VALUE (-1.0e200)
#define default_unset(p) (p->deflt == UNSET_VALUE)

/* structure representing sample information at start of
   a function call */

struct obsinfo_ {
    int structure;      /* time-series, etc. */
    int pd;             /* data frequency */
    int t1, t2;         /* starting and ending observations */
    char changed;       /* sample has been changed within the function call? */
    char stobs[OBSLEN]; /* string representation of starting obs */
};

/* structure representing a call to a user-defined function */

struct fncall_ {
    ufunc *fun;    /* the function called */
    int argc;      /* argument count */
    fn_arg *args;  /* argument array */
    int *ptrvars;  /* list of pointer arguments */
    int *listvars; /* list of series included in a list argument */
    char *retname; /* name of return value (or dummy string) */
    int recursing; /* indicator for recursive call */
    obsinfo obs;   /* sample info */
};

/* structure representing a user-defined function */

struct ufunc_ {
    char name[FN_NAMELEN]; /* identifier */
    fnpkg *pkg;            /* pointer to parent package, or NULL */
    int pkg_role;          /* printer, plotter, etc. */
    UfunAttrs flags;       /* private, plugin, etc. */
    int line_idx;          /* current line index (compiling) */
    int n_lines;           /* number of lines of code */
    fn_line *lines;        /* array of lines of code */
    int n_params;          /* number of parameters */
    fn_param *params;      /* parameter info array */
    int argc;              /* argument count (call-dependent) */
    fn_arg *args;          /* argument array */
    int rettype;           /* return type (if any) */
    int debug;             /* are we debugging this function? */
};

/* structure representing a function package */

struct fnpkg_ {
    int ID;           /* unique ID (time stamp) */
    char name[FN_NAMELEN]; /* package name */
    char *fname;      /* filename */
    char *author;     /* author's name */
    char *email;      /* author's email address */
    char *version;    /* package version string */
    char *date;       /* last revision date */
    char *descrip;    /* package description */
    char *help;       /* package help text */
    char *gui_help;   /* GUI-specific help (optional) */
    char *sample;     /* sample caller script */
    char *help_fname;     /* filename: package help text */
    char *gui_help_fname; /* filename: GUI-specific help text */
    char *sample_fname;   /* filename: sample caller script */
    char *tags;       /* tag string(s) */
    char *label;      /* for use in GUI menus */
    char *mpath;      /* menu path in GUI */
    int minver;       /* minimum required gretl version */
    int uses_subdir;  /* lives in subdirectory (0/1) */
    FuncDataReq dreq; /* data requirement */
    int modelreq;     /* required model type, if applicable */
    ufunc **pub;      /* pointers to public interfaces */
    ufunc **priv;     /* pointers to private functions */
    int n_pub;        /* number of public functions */
    int n_priv;       /* number of private functions */
    char **datafiles; /* names of packaged data files */
    int n_files;      /* number of data files */
};

/* acceptable types for parameters of user-defined functions */

#define ok_function_arg_type(t) (t == GRETL_TYPE_BOOL ||		\
				 t == GRETL_TYPE_INT ||			\
				 t == GRETL_TYPE_OBS ||			\
				 t == GRETL_TYPE_DOUBLE ||		\
				 t == GRETL_TYPE_SERIES ||		\
				 t == GRETL_TYPE_LIST ||		\
				 t == GRETL_TYPE_MATRIX ||		\
				 t == GRETL_TYPE_STRING ||		\
				 t == GRETL_TYPE_BUNDLE ||		\
				 t == GRETL_TYPE_SCALAR_REF ||		\
				 t == GRETL_TYPE_SERIES_REF ||		\
				 t == GRETL_TYPE_MATRIX_REF ||		\
				 t == GRETL_TYPE_BUNDLE_REF ||		\
				 t == GRETL_TYPE_STRINGS ||		\
				 t == GRETL_TYPE_MATRICES ||		\
				 t == GRETL_TYPE_BUNDLES||		\
				 t == GRETL_TYPE_LISTS ||		\
				 t == GRETL_TYPE_STRINGS_REF ||		\
				 t == GRETL_TYPE_MATRICES_REF ||	\
				 t == GRETL_TYPE_BUNDLES_REF)

enum {
    ARG_OPTIONAL = 1 << 0,
    ARG_CONST    = 1 << 1
};

/* structure representing an argument to a user-defined function */

struct fn_arg_ {
    char type;             /* argument type */
    char flags;            /* ARG_OPTIONAL, ARG_CONST as appropriate */
    const char *name;      /* name as function param */
    char upname[VNAMELEN]; /* name of supplied arg at caller level */
    union {
	int idnum;        /* named series arg (series ID) */
	double x;         /* scalar arg */
	double *px;       /* anonymous series arg */
	gretl_matrix *m;  /* matrix arg */
	char *str;        /* string arg */
	int *list;        /* list arg */
	gretl_bundle *b;  /* anonymous bundle pointer */
	gretl_array *a;   /* array argument */
    } val;
};

static int n_ufuns;         /* number of user-defined functions in memory */
static ufunc **ufuns;       /* array of pointers to user-defined functions */
static ufunc *current_fdef; /* pointer to function currently being defined */
static GList *callstack;    /* stack of function calls */
static int n_pkgs;          /* number of loaded function packages */
static fnpkg **pkgs;        /* array of pointers to loaded packages */
static fnpkg *current_pkg;  /* pointer to package currently being edited */

static int function_package_record (fnpkg *pkg);
static void function_package_free (fnpkg *pkg);
static void free_args_array (fn_arg *args, int n);

/* record of state, and communication of state with outside world */

static int compiling;    /* boolean: are we compiling a function currently? */
static int fn_executing; /* depth of function call stack */
static int compiling_python;

#define function_is_private(f)   (f->flags & UFUN_PRIVATE)
#define function_is_plugin(f)    (f->flags & UFUN_PLUGIN)
#define function_is_noprint(f)   (f->flags & UFUN_NOPRINT)
#define function_is_menu_only(f) (f->flags & UFUN_MENU_ONLY)

struct flag_and_key {
    int flag;
    const char *key;
};

static struct flag_and_key pkg_lookups[] = {
    { UFUN_BUNDLE_PRINT, BUNDLE_PRINT },
    { UFUN_BUNDLE_PLOT,  BUNDLE_PLOT },
    { UFUN_BUNDLE_TEST,  BUNDLE_TEST },
    { UFUN_BUNDLE_FCAST, BUNDLE_FCAST },
    { UFUN_BUNDLE_EXTRA, BUNDLE_EXTRA },
    { UFUN_GUI_MAIN,     GUI_MAIN },
    { UFUN_GUI_PRECHECK, GUI_PRECHECK },
    { -1,                NULL }
};

#define pkg_aux_role(r) (r == UFUN_BUNDLE_PRINT || \
			 r == UFUN_BUNDLE_PLOT ||  \
			 r == UFUN_BUNDLE_TEST ||  \
			 r == UFUN_BUNDLE_FCAST || \
			 r == UFUN_BUNDLE_EXTRA)

static int pkg_key_get_role (const char *key)
{
    int i;

    if (key != NULL && *key != '\0') {
	for (i=0; pkg_lookups[i].flag > 0; i++) {
	    if (!strcmp(key, pkg_lookups[i].key)) {
		return pkg_lookups[i].flag;
	    }
	}
    }
    
    return UFUN_ROLE_NONE;
}

const char *package_role_get_key (int flag)
{
    int i;

    for (i=0; pkg_lookups[i].flag > 0; i++) {
	if (flag == pkg_lookups[i].flag) {
	    return pkg_lookups[i].key;
	}
    }
    
    return NULL;
}

static void set_function_private (ufunc *u, gboolean s)
{
    if (s) {
	u->flags |= UFUN_PRIVATE;
    } else {
	u->flags &= ~UFUN_PRIVATE;
    }
}

int gretl_compiling_function (void)
{
    return compiling;
}

int gretl_compiling_python (const char *line)
{
    if (compiling_python) {
	char s1[4], s2[8];

	if (sscanf(line, "%3s %7s", s1, s2) == 2 && 
	    !strcmp(s1, "end") && !strcmp(s2, "foreign")) {
	    compiling_python = 0;
	}
    }

    return compiling_python;
}

static void set_compiling_on (void)
{
    compiling = 1;
}

static void set_compiling_off (void)
{
    compiling = compiling_python = 0;
}

int gretl_function_depth (void)
{
    return fn_executing;
}

static void adjust_array_arg_type (fn_arg *arg)
{
    GretlType t = gretl_array_get_type(arg->val.a);

    if (arg->type == GRETL_TYPE_ARRAY_REF) {
	arg->type = gretl_type_get_ref_type(t);
    } else {
	arg->type = t;
    }
}

static int fn_arg_set_data (fn_arg *arg, const char *name,
			    GretlType type, void *p)
{
    int err = 0;
    
    arg->type = type;
    arg->flags = 0;
    arg->name = NULL;
    arg->upname[0] = '\0';

    if (name != NULL) {
	strcpy(arg->upname, name);
    }
    
    if (type == GRETL_TYPE_NONE) {
	arg->val.x = 0;
    } else if (type == GRETL_TYPE_DOUBLE || 
	       type == GRETL_TYPE_SCALAR_REF) {
	arg->val.x = *(double *) p;
    } else if (type == GRETL_TYPE_INT || 
	       type == GRETL_TYPE_OBS) {
	arg->val.x = *(int *) p;
    } else if (type == GRETL_TYPE_SERIES) {
	arg->val.px = (double *) p;
    } else if (type == GRETL_TYPE_MATRIX || 
	       type == GRETL_TYPE_MATRIX_REF) {
	arg->val.m = (gretl_matrix *) p;
    } else if (type == GRETL_TYPE_STRING) {
	arg->val.str = (char *) p;
    } else if (type == GRETL_TYPE_LIST) {
	arg->val.list = (int *) p;
    } else if (type == GRETL_TYPE_SERIES_REF ||
	       type == GRETL_TYPE_USERIES) {
	arg->val.idnum = *(int *) p;
    } else if (type == GRETL_TYPE_BUNDLE || 
	       type == GRETL_TYPE_BUNDLE_REF) {
	arg->val.b = (gretl_bundle *) p;
    } else if (type == GRETL_TYPE_ARRAY ||
	       type == GRETL_TYPE_ARRAY_REF) {
	arg->val.a = (gretl_array *) p;
	adjust_array_arg_type(arg);
    } else {
	err = E_TYPES;
    }

    return err;
}

static int ufunc_add_args_array (ufunc *u)
{
    int i, err = 0;
    
    u->args = malloc(u->n_params * sizeof *u->args);

    if (u->args == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<u->n_params; i++) {
	    u->args[i].type = 0;
	    u->args[i].flags = 0;
	    u->args[i].name = NULL;
	    u->args[i].upname[0] = '\0';
	}
    }

    return err;
}

/**
 * push_function_arg:
 * @fun: pointer to function.
 * @name: name of variable (or NULL for anonymous)
 * @type: type of argument to add.
 * @value: pointer to value to add.
 *
 * Writes a new argument of the specified type and value into the
 * argument array of @fun.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int push_function_arg (ufunc *fun, const char *name, GretlType type, 
		       void *value)
{
    int err = 0;

    if (fun == NULL) {
	err = E_DATA;
    } else if (fun->argc >= fun->n_params) {
	fprintf(stderr, "function %s has %d parameters but argc = %d\n",
		fun->name, fun->n_params, fun->argc);
	err = E_DATA;
    } else if (fun->args == NULL) {
	err = ufunc_add_args_array(fun);
    }

    if (!err) {
	err = fn_arg_set_data(&fun->args[fun->argc], name, type, value);
	fun->argc += 1;
    }

    return err;
}

void function_clear_args (ufunc *fun)
{
    if (fun != NULL && fun->args != NULL) {
	int i;

	for (i=0; i<fun->argc; i++) {
	    fun->args[i].type = 0;
	    fun->args[i].flags = 0;
	    fun->args[i].name = NULL;
	    fun->args[i].upname[0] = '\0';
	}

	fun->argc = 0;
    }
}

static fncall *fncall_new (ufunc *fun)
{
    fncall *call = malloc(sizeof *call);

    if (call != NULL) {
	call->fun = fun;
	call->ptrvars = NULL;
	call->listvars = NULL;
	call->retname = NULL;
	/* Allow for recursion: transfer the arg list
	   from @fun to @call so that it doesn't
	   overflow if @fun is called again. 
	*/
	call->argc = fun->argc;
	call->args = fun->args;
	call->recursing = 0;
	fun->args = NULL;
	fun->argc = 0;
    }

    return call;
}

static void fncall_free (fncall *call)
{
    if (call != NULL) {
	if (call->fun != NULL && call->args != NULL) {
	    if (call->recursing) {
		/* free (duplicate) args array */
		free_args_array(call->args, call->argc);
	    } else {
		/* hand args array back to @fun for
		   possible reuse */
		call->fun->args = call->args;
		function_clear_args(call->fun);
	    }
	}
	free(call->ptrvars);
	free(call->listvars);
	free(call->retname);
	free(call);
    }
}

static fnpkg *function_package_alloc (const char *fname)
{
    fnpkg *pkg = malloc(sizeof *pkg);

    if (pkg == NULL) {
	return NULL;
    }

    pkg->fname = gretl_strdup(fname);
    if (pkg->fname == NULL) {
	free(pkg);
	return NULL;
    }

#if PKG_DEBUG
    fprintf(stderr, "function_package_alloc: fname='%s'\n", fname);
#endif    

    pkg->ID = (int) time(NULL);
    pkg->name[0] = '\0';
    pkg->author = NULL;
    pkg->email = NULL;
    pkg->version = NULL;
    pkg->date = NULL;
    pkg->descrip = NULL;
    pkg->help = NULL;
    pkg->gui_help = NULL;
    pkg->sample = NULL;
    pkg->help_fname = NULL;
    pkg->gui_help_fname = NULL;
    pkg->sample_fname = NULL;
    pkg->tags = NULL;
    pkg->label = NULL;
    pkg->mpath = NULL;
    pkg->dreq = 0;
    pkg->modelreq = 0;
    pkg->minver = 0;
    pkg->uses_subdir = 0;
    
    pkg->pub = pkg->priv = NULL;
    pkg->n_pub = pkg->n_priv = 0;
    pkg->datafiles = NULL;
    pkg->n_files = 0;

    return pkg;
}

/* For function call @call, set the 'listarg' status of any named
   series provided to the called function via list arguments.  This is
   required for correct handling of namespaces.  
*/

static void set_listargs_from_call (fncall *call, DATASET *dset)
{
    int i, vi;

    if (dset == NULL) {
	return;
    }

    for (i=1; i<dset->v; i++) {
	series_unset_flag(dset, i, VAR_LISTARG);
    }

    if (call != NULL && call->listvars != NULL) {
	for (i=1; i<=call->listvars[0]; i++) {
	    vi = call->listvars[i];
#if UDEBUG
	    fprintf(stderr, "setting listarg status on var %d (%s)\n",
		    vi, dset->varname[vi]);
#endif
	    series_set_flag(dset, vi, VAR_LISTARG);
	}
    }
}

static void set_executing_off (fncall *call, DATASET *dset, PRN *prn)
{
    int dbg = gretl_debugging_on();
    fncall *popcall = NULL;

    destroy_option_params_at_level(fn_executing);
    fn_executing--;

    callstack = g_list_remove(callstack, call);

#if EXEC_DEBUG
    fprintf(stderr, "set_executing_off: removing call to %s, depth now %d\n",
	    call->fun->name, g_list_length(callstack));
#endif

    if (dbg) {
	pputs(prn, "*** ");
	bufspace(gretl_function_depth(), prn);
	pprintf(prn, "exiting function %s, ", call->fun->name);
    }

    fncall_free(call);

    if (fn_executing > 0) {
	GList *tmp = g_list_last(callstack);
	
	popcall = tmp->data;
    } else {
	g_list_free(callstack);
	callstack = NULL;
	gretl_insert_builtin_string("pkgdir", NULL);
    }

    if (dset != NULL) {
	set_listargs_from_call(popcall, dset);
    }

    if (popcall == NULL) {
	/* returning to main */
	switch_uservar_hash(0);
    }

    if (dbg) {
	if (popcall != NULL) {
	    pprintf(prn, "returning to %s\n", popcall->fun->name);
	} else {
	    pputs(prn, "returning to main\n");
	}
    }
}

/**
 * n_user_functions:
 *
 * Returns: the number of hansl functions currently loaded in memory.
 */

int n_user_functions (void)
{
    return n_ufuns;
}

/**
 * n_free_functions:
 *
 * Returns: the number of functions loaded in memory
 * that are not currently attached to any function package,
 * and are therefore available for packaging.
 */

int n_free_functions (void)
{
    int i, n = 0;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == NULL) {
	    n++;
	}
    }

    return n;
}

/**
 * get_user_function_by_index:
 * @idx: index number.
 *
 * Returns: pointer to the user-function that currently
 * occupies (0-based) slot @idx in the array of loaded
 * functions, or %NULL if @idx is out of bounds.
 */

const ufunc *get_user_function_by_index (int idx)
{
    return (idx < 0 || idx >= n_ufuns)? NULL : ufuns[idx];
}

/**
 * fn_n_params:
 * @fun: pointer to user-function.
 *
 * Returns: the number of parameters associated with @fun.
 */

int fn_n_params (const ufunc *fun)
{
    return fun->n_params;
}

/**
 * fn_param_type:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the type of parameter @i of function
 * @fun.
 */

int fn_param_type (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? 0 :
	fun->params[i].type;
}

/**
 * fn_param_name:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the name of parameter @i of function
 * @fun.
 */

const char *fn_param_name (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NULL :
	fun->params[i].name;
}

/**
 * fn_param_descrip:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the description of parameter @i of function
 * @fun (if any), otherwise %NULL.
 */

const char *fn_param_descrip (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NULL :
	fun->params[i].descrip;
}

/**
 * fn_param_value_labels:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * @n: location to receive number of labels.
 * 
 * Returns: the value-labels associated with parameter @i 
 * of function @fun (if any), otherwise %NULL.
 */

const char **fn_param_value_labels (const ufunc *fun, int i, 
				    int *n)
{
    if (i >= 0 && i < fun->n_params) {
	*n = fun->params[i].nlabels;
	return (const char **) fun->params[i].labels;
    } else {
	*n = 0;
	return NULL;
    }
}

/**
 * fn_param_has_default:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: 1 if the (scalar) parameter @i of function @fun 
 * (if any) has a default value set, otherwise 0.
 */

int fn_param_has_default (const ufunc *fun, int i)
{
    if (i < 0 || i >= fun->n_params) {
	return 0;
    } else {
	fn_param *fp = &fun->params[i];

	return !default_unset(fp);
    }
}   

/**
 * fn_param_default:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the default value of (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

double fn_param_default (const ufunc *fun, int i)
{
    if (i < 0 || i >= fun->n_params) {
	return NADBL;
    } else {
	fn_param *fp = &fun->params[i];

	return default_unset(fp)? NADBL : fp->deflt;
    }
}    

/**
 * fn_param_minval:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the minimum value of (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

double fn_param_minval (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NADBL :
	fun->params[i].min;
}

/**
 * fn_param_maxval:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the maximum value of (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

double fn_param_maxval (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NADBL :
	fun->params[i].max;
}

/**
 * fn_param_step:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the step value for (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

double fn_param_step (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NADBL :
	fun->params[i].step;
}

static int arg_may_be_optional (GretlType t)
{
    return gretl_ref_type(t) ||
	t == GRETL_TYPE_MATRIX ||
	t == GRETL_TYPE_LIST ||
	t == GRETL_TYPE_STRING;	
}

/**
 * fn_param_optional:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: 1 if parameter @i of function @fun is optional,
 * otherwise 0.
 */

int fn_param_optional (const ufunc *fun, int i)
{
    if (i < 0 || i >= fun->n_params) {
	return 0;
    }

    return arg_may_be_optional(fun->params[i].type) &&
	(fun->params[i].flags & ARG_OPTIONAL);
}

/**
 * fn_param_uses_xlist:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: 1 if parameter @i of function @fun is
 * designed to select an integer based on a gretl
 * model's list of regressors, otherwise 0.
 */

int fn_param_uses_xlist (const ufunc *fun, int i)
{
    if (i < 0 || i >= fun->n_params) {
	return 0;
    }

    return (fun->params[i].type == GRETL_TYPE_INT &&
	    fun->params[i].deflt == INT_USE_XLIST);
}

/**
 * user_func_get_return_type:
 * @fun: pointer to user-function.
 * 
 * Returns: the return type of function @fun.
 */

int user_func_get_return_type (const ufunc *fun)
{
    if (fun == NULL) {
	return GRETL_TYPE_NONE;
    } else {
	return fun->rettype;
    }
}

/**
 * user_func_is_noprint:
 * @fun: pointer to user-function.
 * 
 * Returns: 1 if the function is not designed to print anything.
 */

int user_func_is_noprint (const ufunc *fun)
{
    if (fun == NULL) {
	return 0;
    } else {
	return function_is_noprint(fun);
    }
}

/**
 * user_func_is_menu_only:
 * @fun: pointer to user-function.
 * 
 * Returns: 1 if the function is not designed to be called
 * other than via a GUI menu.
 */

int user_func_is_menu_only (const ufunc *fun)
{
    if (fun == NULL) {
	return 0;
    } else {
	return function_is_menu_only(fun);
    }
}

/**
 * user_function_name_by_index:
 * @i: the position of a user-function in the array of
 * loaded functions.
 * 
 * Returns: the name of the function, or %NULL if
 * @i is out of bounds.
 */

const char *user_function_name_by_index (int i)
{
    if (i >= 0 && i < n_ufuns) {
	return ufuns[i]->name;
    } else {
	return NULL;
    }
}

/**
 * user_function_index_by_name:
 * @name: function name.
 * @pkg: reference function package, or NULL.
 *
 * Looks up the 0-based index of the named function
 * in the current array of user-functions. If @pkg is 
 * non-NULL the search for @name is confined to functions
 * that belong to @pkg; otherwise the index of the first
 * match is returned.
 * 
 * Returns: 0-based index or -1 on failure.
 */

int user_function_index_by_name (const char *name,
				 fnpkg *pkg)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if ((pkg == NULL || ufuns[i]->pkg == pkg) && 
	    !strcmp(name, ufuns[i]->name)) {
	    return i;
	}
    }

    return -1;
}	

static int fname_idx;

void function_names_init (void)
{
    fname_idx = 0;
}

/* Apparatus used in the GUI selector for composing a new function
   package, or editing an existing one.  In the first case we want a
   list of names of currently unpackaged functions (and the @pkg
   argument should be NULL); in the second the list should include
   both unpackaged functions and those belonging to the package
   in question (specified via the @pkg arg).

   The caller should first invoke function_names_init(), then keep
   calling next_available_function_name() until it returns %NULL. The
   pointer argument @idxp provides a means to grab the "index number"
   (position in the current functions array) corresponding to the
   returned function name.
*/

const char *next_available_function_name (fnpkg *pkg, int *idxp)
{
    const char *ret = NULL;
    ufunc *fun;

    if (n_ufuns == 0) {
	fname_idx = 0;
	return NULL;
    }

    while (fname_idx < n_ufuns) {
	fun = ufuns[fname_idx++];
	if (fun->pkg == NULL || fun->pkg == pkg) {
	    ret = fun->name;
	    *idxp = fname_idx - 1;
	    break;
	}
    }

    return ret;
}

/* end GUI function-name apparatus */

static fncall *current_function_call (void)
{
    if (callstack != NULL) {
	GList *tmp = g_list_last(callstack);

	return tmp->data;
    } else {
	return NULL;
    }
}

static ufunc *currently_called_function (void)
{
    fncall *call = current_function_call();

    return (call != NULL)? call->fun : NULL;
}

void current_function_info (char const **funcname,
			    char const **pkgname)
{
    ufunc *u = currently_called_function();

    if (u != NULL) {
	*funcname = u->name;
	if (pkgname != NULL && u->pkg != NULL) {
	    *pkgname = u->pkg->name;
	}
    }
}

/* see if a function is currently employed in the call stack */

static int function_in_use (ufunc *fun)
{
    GList *tmp = callstack;
    fncall *call;

    while (tmp != NULL) {
	call = tmp->data;
	if (fun == call->fun) {
	    return 1;
	}
	tmp = tmp->next;
    }

    return 0;
}

int gretl_function_recursing (void)
{
    if (callstack != NULL) {
	GList *tmp = g_list_last(callstack);
	fncall *call = tmp->data;

	return call->recursing;
    } else {
	return 0;
    }
}

/**
 * get_ufunc_by_name:
 * @name: name to test.
 *
 * Returns: pointer to a user-function, if there exists a
 * function of the given name and it is accessible in
 * context (i.e. it's not private to a package other than
 * the one that's currently active, if any), otherwise
 * %NULL.
 */

ufunc *get_user_function_by_name (const char *name)
{
    fnpkg *pkg = current_pkg;
    ufunc *fun = NULL;
    int i;

    if (n_ufuns == 0) {
	return NULL;
    }

    if (pkg == NULL) {
	fun = currently_called_function();
	if (fun != NULL) {
	    pkg = fun->pkg;
	    fun = NULL;
	}
    }

    /* First pass: if there's no active function package, match any
       non-private function, but if there is an active package match
       only its member functions.  Try to optimize by putting the
       cheapest comparisons first.
    */

    for (i=0; i<n_ufuns; i++) {
	if ((pkg == NULL && !function_is_private(ufuns[i])) || ufuns[i]->pkg == pkg) {
	    if (!strcmp(name, ufuns[i]->name)) { 
		fun = ufuns[i];
		break;
	    }
	}
    }

    /* Second pass, if the first didn't work: match unpackaged
       functions or functions from other packages, so long as they're
       not private.
    */

    if (fun == NULL && pkg != NULL) {
	for (i=0; i<n_ufuns; i++) {
	    if (!function_is_private(ufuns[i]) && !strcmp(name, ufuns[i]->name)) {
		fun = ufuns[i];
		break;
	    }
	}
    }	

#if FN_DEBUG > 1
    if (fun != NULL) {
	fprintf(stderr, "get_user_function_by_name: name = '%s' (n_ufuns = %d);"
		" found match\n", name, n_ufuns);
    }
#endif

    return fun;
}

/**
 * get_function_from_package:
 * @funname: name if fuction to retrieve.
 * @pkg: function package.
 *
 * Returns: pointer to a user-function, if there exists a
 * function of the given @funname that is associated with
 * function package @pkg, otherwise NULL.  This is used
 * in the gretl function package editor.
 */

ufunc *get_function_from_package (const char *funname, fnpkg *pkg)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == pkg && 
	    !strcmp(funname, ufuns[i]->name)) {
	    return ufuns[i];
	}
    }

    return NULL;
}

/* allocate and initialize a new array of @n parameters */

static fn_param *allocate_params (int n)
{
    fn_param *params;
    int i;

    params = malloc(n * sizeof *params);
    if (params == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	params[i].name = NULL;
	params[i].type = 0;
	params[i].descrip = NULL;
	params[i].labels = NULL;
	params[i].nlabels = 0;
	params[i].flags = 0;
	params[i].deflt = UNSET_VALUE;
	params[i].min = NADBL;
	params[i].max = NADBL;
	params[i].step = NADBL;
    }

    return params;
}

static ufunc *ufunc_new (void)
{
    ufunc *fun = malloc(sizeof *fun);

    if (fun == NULL) {
	return NULL;
    }

    fun->name[0] = '\0';
    fun->pkg = NULL;
    fun->flags = 0;
    fun->pkg_role = 0;

    fun->n_lines = 0;
    fun->line_idx = 1;
    fun->lines = NULL;

    fun->n_params = 0;
    fun->params = NULL;
    fun->argc = 0;
    fun->args = NULL;

    fun->rettype = GRETL_TYPE_NONE;

    fun->debug = 0;

    return fun;
}

static void free_lines_array (fn_line *lines, int n)
{
    int i;

    if (lines == NULL) return;

    for (i=0; i<n; i++) {
	free(lines[i].s);
	if (lines[i].loop != NULL) {
	    gretl_loop_destroy(lines[i].loop);
	}
    }
    
    free(lines);
}

static void free_params_array (fn_param *params, int n)
{
    int i;

    if (params == NULL) return;

    for (i=0; i<n; i++) {
	free(params[i].name);
	free(params[i].descrip);
	strings_array_free(params[i].labels, params[i].nlabels);
    }
    free(params);
}

static void free_args_array (fn_arg *args, int n)
{
    if (args != NULL) {
	free(args);
    }
}

static void clear_ufunc_data (ufunc *fun)
{
    free_lines_array(fun->lines, fun->n_lines);
    free_params_array(fun->params, fun->n_params);
    free_args_array(fun->args, fun->n_params);
    
    fun->lines = NULL;
    fun->params = NULL;
    fun->args = NULL;

    fun->n_lines = 0;
    fun->line_idx = 1;
    fun->n_params = 0;
    fun->argc = 0;

    fun->rettype = GRETL_TYPE_NONE;
}

static void ufunc_free (ufunc *fun)
{
    free_lines_array(fun->lines, fun->n_lines);
    free_params_array(fun->params, fun->n_params);
    free_args_array(fun->args, fun->n_params);

    free(fun);
}

static int add_allocated_ufunc (ufunc *fun)
{
    int nf = n_ufuns;
    ufunc **myfuns;

    myfuns = realloc(ufuns, (nf + 1) * sizeof *myfuns);

    if (myfuns == NULL) {
	return E_ALLOC;
    }

    ufuns = myfuns;
    ufuns[nf] = fun;
    n_ufuns++;

#if PKG_DEBUG
    fprintf(stderr, "add_allocated_ufunc: name '%s', n_ufuns = %d\n",
	    fun->name, n_ufuns);
#endif

    return 0;
}

static ufunc *add_ufunc (const char *fname)
{
    ufunc *fun = ufunc_new();

#if FN_DEBUG
    fprintf(stderr, "add_ufunc: '%s'\n", fname);
#endif

    if (fun != NULL) {
	strncat(fun->name, fname, FN_NAMELEN - 1);
	if (add_allocated_ufunc(fun)) {
	    ufunc_free(fun);
	    fun = NULL;
	}
    }

    return fun;
}

static int no_scalar_default (fn_param *fp)
{
    int ret = 0;

    if (default_unset(fp)) {
	ret = 1;
    } else if (fp->type != GRETL_TYPE_DOUBLE && na(fp->deflt)) {
	ret = 1;
    } 

    return ret;
}

/* handling of XML function packages */

enum {
    FUNCS_INFO,
    FUNCS_LOAD,
    FUNCS_CODE,
    FUNCS_SAMPLE
};

static const char *arg_type_xml_string (int t)
{
    if (t == GRETL_TYPE_SCALAR_REF) {
	return "scalarref";
    } else if (t == GRETL_TYPE_SERIES_REF) {
	return "seriesref";
    } else if (t == GRETL_TYPE_MATRIX_REF) {
	return "matrixref";
    } else if (t == GRETL_TYPE_BUNDLE_REF) {
	return "bundleref";
    } else if (t == GRETL_TYPE_STRINGS_REF) {
	return "stringsref";
    } else if (t == GRETL_TYPE_MATRICES_REF) {
	return "matricesref";
    } else if (t == GRETL_TYPE_BUNDLES_REF) {
	return "bundlesref";
    } else if (t == GRETL_TYPE_LISTS_REF) {
	return "listsref";
    } else {
	return gretl_type_get_name(t);
    }
}

static GretlType return_type_from_string (const char *s,
					  int *err)
{
    GretlType t;

    if (!strcmp(s, "void")) {
	/* not OK as arg type, but OK as return */
	t = GRETL_TYPE_VOID;
    } else {
	t = gretl_type_from_string(s);
    }

    if (!ok_function_return_type(t)) {
	if (*s == '\0') {
	    gretl_errmsg_sprintf(_("Missing function return type"));
	} else if (t == GRETL_TYPE_NONE) {
	    gretl_errmsg_sprintf(_("Expected a function return type, found '%s'"),
				 s);
	} else {
	    gretl_errmsg_sprintf(_("%s: invalid return type for function"),
				 s);
	}
	*err = E_TYPES;
    }

    return t;
}

static GretlType param_field_to_type (const char *s, 
				      const char *funname,
				      int *err)
{
    GretlType t = gretl_type_from_string(s);

    if (!ok_function_arg_type(t)) {
	gretl_errmsg_sprintf("function %s: invalid parameter type '%s'", 
			     funname, s);
	*err = E_INVARG;
    }

    return t;
}    

/* read the parameter info for a function from XML file */

static int func_read_params (xmlNodePtr node, xmlDocPtr doc,
			     ufunc *fun)
{
    xmlNodePtr cur;
    char *field;
    int n, err = 0;

    if (!gretl_xml_get_prop_as_int(node, "count", &n) || n < 0) {
	fprintf(stderr, "Couldn't read param count\n");
	return E_DATA;
    }  

    if (n == 0) {
	return 0;
    }

    fun->params = allocate_params(n);
    if (fun->params == NULL) {
	return E_ALLOC;
    } 

    fun->n_params = n;

    gretl_push_c_numeric_locale();

    cur = node->xmlChildrenNode;
    n = 0;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "param")) {
	    fn_param *param = &fun->params[n++];

	    if (gretl_xml_get_prop_as_string(cur, "name", &field)) {
		param->name = field;
	    } else {
		err = E_DATA;
		break;
	    }
	    if (gretl_xml_get_prop_as_string(cur, "type", &field)) {
		param->type = param_field_to_type(field, fun->name, &err);
		free(field);
		if (gretl_scalar_type(param->type)) {
		    double x;

		    if (gretl_xml_get_prop_as_double(cur, "default", &x)) {
			param->deflt = x;
		    } else {
			param->deflt = UNSET_VALUE;
		    } 
		    if (param->type != GRETL_TYPE_BOOL) {
			gretl_xml_get_prop_as_double(cur, "min", &param->min);
			gretl_xml_get_prop_as_double(cur, "max", &param->max);
			gretl_xml_get_prop_as_double(cur, "step", &param->step);
		    }
		}
		if (gretl_xml_get_prop_as_bool(cur, "optional")) {
		    param->flags |= ARG_OPTIONAL;
		}
		if (gretl_xml_get_prop_as_bool(cur, "const")) {
		    param->flags |= ARG_CONST;
		}
	    } else {
		err = E_DATA;
		break;
	    }
	    gretl_xml_child_get_string(cur, doc, "description", 
				       &param->descrip);
	    gretl_xml_child_get_strings_array(cur, doc, "labels",
					      &param->labels,
					      &param->nlabels);
	}	    
	cur = cur->next;
    }

    gretl_pop_c_numeric_locale();

    if (!err && n != fun->n_params) {
	err = E_DATA;
    }

    return err;
}

static int push_function_line (ufunc *fun, const char *s)
{
    int err = 0;

    if (string_is_blank(s)) {
	fun->line_idx += 1;
    } else {
	fn_line *lines;
	int n = fun->n_lines + 1;

	lines = realloc(fun->lines, n * sizeof *lines);

	if (lines == NULL) {
	    err = E_ALLOC;
	} else {
	    int i = n - 1;
	    
	    fun->lines = lines;
	    lines[i].idx = fun->line_idx;
	    lines[i].s = gretl_strdup(s);
	    lines[i].loop = NULL;
	    lines[i].next_idx = -1;
	    lines[i].ignore = 0;
	    if (lines[i].s == NULL) {
		err = E_ALLOC;
	    } else {
		fun->n_lines = n;
		fun->line_idx += 1;
	    }
	}
    }

    return err;
}

/* read the actual code lines from the XML representation of a 
   function */

static int func_read_code (xmlNodePtr node, xmlDocPtr doc, ufunc *fun)
{
    char line[MAXLINE];
    char *buf, *s;
    int err = 0;

    buf = (char *) xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
    if (buf == NULL) {
	return 1;
    }

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf) && !err) {
	s = line;
	while (isspace(*s)) s++;
	tailstrip(s);
	err = push_function_line(fun, s);
    }

    bufgets_finalize(buf);

    free(buf);

    return err;
}

static void print_opt_flags (fn_param *param, PRN *prn)
{
    if (param->flags & ARG_OPTIONAL) {
	pputs(prn, "[null]");
    }
}

static void print_param_description (fn_param *param, PRN *prn)
{
    if (param->descrip != NULL && *param->descrip != '\0') {
	pprintf(prn, " \"%s\"", param->descrip);
    }
}

static void print_param_labels (fn_param *param, PRN *prn)
{
    int i;

    pputs(prn, " {");

    for (i=0; i<param->nlabels; i++) {
	pprintf(prn, "\"%s\"", param->labels[i]);
	if (i < param->nlabels - 1) {
	    pputs(prn, ", ");
	}
    }

    pputc(prn, '}');
}

static void print_min_max_deflt (fn_param *param, PRN *prn)
{
    if (na(param->min) && na(param->max) && default_unset(param)) {
	return; /* no-op */
    } else if (na(param->min) && na(param->max)) {
	/* got a default value only? */
	if (!default_unset(param)) {
	    if (na(param->deflt)) {
		pputs(prn, "[NA]");
	    } else {
		pprintf(prn, "[%g]", param->deflt);
	    }
	}
	return;
    }

    pputc(prn, '[');

    /* minimum */
    if (!na(param->min)) pprintf(prn, "%g", param->min);
    pputc(prn, ':');

    /* maximum */
    if (!na(param->max)) pprintf(prn, "%g", param->max);
    pputc(prn, ':');

    /* default */
    if (!default_unset(param)) {
	if (na(param->deflt)) {
	    pputs(prn, "NA");
	} else {
	    pprintf(prn, "%g", param->deflt);
	}
    }

    if (!na(param->step)) {
	/* step */
	pputc(prn, ':');
	pprintf(prn, "%g", param->step);
    }

    pputc(prn, ']');
}

/* free @fun and also remove it from the list of loaded
   functions */

static void ufunc_unload (ufunc *fun)
{
    int i, j, found = 0;

    if (n_ufuns == 0 || fun == NULL) {
	/* "can't happen" */
	return;
    }

    /* remove this function from the array of loaded functions */

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i] == fun) {
	    for (j=i; j<n_ufuns-1; j++) {
		ufuns[j] = ufuns[j+1];
	    }
	    found = 1;
	    break;
	}
    }

    ufunc_free(fun);

    if (found) {
	n_ufuns--;
    }
}

/* remove a package from the current listing and free it; 
   if @full is non-zero, also unload any member functions
*/

static void real_function_package_unload (fnpkg *pkg, int full)
{
    int i, j, found = 0;

#if PKG_DEBUG
    fprintf(stderr, "function_package_unload: unloading '%s'\n", pkg->name);
#endif

    /* unload children? */
    if (full) {
	for (i=0; i<n_ufuns; i++) {
	    if (ufuns[i]->pkg == pkg) {
		ufunc_unload(ufuns[i]);
	    }
	}  
    }  

    /* remove this package from the array of loaded packages */
    for (i=0; i<n_pkgs; i++) {
	if (pkgs[i] == pkg) {
	    for (j=i; j<n_pkgs-1; j++) {
		pkgs[j] = pkgs[j+1];
	    }
	    found = 1;
	    break;
	}
    }

    /* free the package itself */
    function_package_free(pkg);

    if (found) {
	n_pkgs--;
    }
}

static void function_package_unload_full (fnpkg *pkg)
{
    real_function_package_unload(pkg, 1);
}

/* Append a pointer to @fun to the array of child-pointers in @pkg: we
   do this when reading the definition of a packaged function from an
   XML file.  Note that this action does not add the functions to the
   array of loaded functions -- that's done separately, if we're
   loading the package 'for real'.
*/

static int attach_ufunc_to_package (ufunc *fun, fnpkg *pkg)
{
    ufunc **uf;
    int n, err = 0;

    if (function_is_private(fun)) {
	n = pkg->n_priv;
	uf = realloc(pkg->priv, (n + 1) * sizeof *uf);
	if (uf == NULL) {
	    err = E_ALLOC;
	} else {
	    pkg->priv = uf;
	    pkg->priv[n] = fun;
	    pkg->n_priv += 1;
	}
    } else {
	n = pkg->n_pub;
	uf = realloc(pkg->pub, (n + 1) * sizeof *uf);
	if (uf == NULL) {
	    err = E_ALLOC;
	} else {
	    pkg->pub = uf;
	    pkg->pub[n] = fun;
	    pkg->n_pub += 1;
	}	
    } 

#if PKG_DEBUG
    fprintf(stderr, "attach_ufunc_to_package: package = %s, "
	    "private = %d, err = %d\n", pkg->name, 
	    function_is_private(fun), err);
#endif

    return err;
}

static void maybe_set_menu_only (ufunc *fun, fnpkg *pkg)
{
    if (pkg->mpath != NULL && strstr(pkg->mpath, "MODELWIN")) {
	if (fun->pkg_role == UFUN_GUI_MAIN) {
	    fun->flags |= UFUN_MENU_ONLY;
	}
    }
}

/* read a single user-function definition from XML file: if the
   function is a child of a package, the @pkg argument will
   be non-NULL
*/

static int read_ufunc_from_xml (xmlNodePtr node, xmlDocPtr doc, fnpkg *pkg)
{
    ufunc *fun = ufunc_new();
    xmlNodePtr cur;
    char *tmp;
    int err = 0;

    if (fun == NULL) {
	return E_ALLOC;
    }

    if (!gretl_xml_get_prop_as_string(node, "name", &tmp)) {
	ufunc_free(fun);
	return E_DATA;
    }

    strncat(fun->name, tmp, FN_NAMELEN - 1);
    free(tmp);

    if (pkg != NULL) {
	fun->pkg = pkg;
    }

    if (gretl_xml_get_prop_as_string(node, "type", &tmp)) {
	fun->rettype = return_type_from_string(tmp, &err);
	free(tmp);
    } else {
	fun->rettype = GRETL_TYPE_VOID;
    }

    if (err) {
	return err;
    }

    if (gretl_xml_get_prop_as_bool(node, "private")) {
	fun->flags |= UFUN_PRIVATE;
    }

    if (gretl_xml_get_prop_as_bool(node, "plugin-wrapper")) {
	fun->flags |= UFUN_PLUGIN;
    }

    if (gretl_xml_get_prop_as_bool(node, "no-print")) {
	fun->flags |= UFUN_NOPRINT;
    }

    if (gretl_xml_get_prop_as_bool(node, "menu-only")) {
	fun->flags |= UFUN_MENU_ONLY;
    }

    if (gretl_xml_get_prop_as_string(node, "pkg-role", &tmp)) {
	fun->pkg_role = pkg_key_get_role(tmp);
	free(tmp);
    }

    if (!(fun->flags & UFUN_MENU_ONLY) && pkg != NULL) {
	maybe_set_menu_only(fun, pkg);
    }

#if PKG_DEBUG
    fprintf(stderr, "read_ufunc_from_xml: name '%s', type %d\n"
	    " private = %d, plugin = %d\n", fun->name, fun->rettype, 
	    function_is_private(fun), function_is_plugin(fun));
#endif

    if (pkg == NULL && (function_is_private(fun) || function_is_plugin(fun))) {
	fprintf(stderr, "unpackaged function: can't be private or plugin\n");
	ufunc_free(fun);
	return E_DATA;
    }	

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "help")) {
	    /* backward compatibility: help used to be attached to
	       a package's public interface */
	    if (pkg->help != NULL) {
		free(pkg->help);
	    }
	    pkg->help = gretl_xml_get_string(cur, doc);
	} else if (!xmlStrcmp(cur->name, (XUC) "params")) {
	    err = func_read_params(cur, doc, fun);
	    if (err) {
		fprintf(stderr, "%s: error parsing function parameters\n",
			fun->name);
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "return")) {
	    gretl_errmsg_set("Old-style function definitions no longer supported");
	    err = E_DATA;
	} else if (!xmlStrcmp(cur->name, (XUC) "code")) {
	    err = func_read_code(cur, doc, fun);
	}
	cur = cur->next;
    }

    if (!err) {
	if (pkg != NULL) {
	    /* function belongs to a package */
	    err = attach_ufunc_to_package(fun, pkg);
	} else {
	    /* reading a standalone function from session file */
	    err = add_allocated_ufunc(fun);
	}
    }

    if (err) {
	ufunc_free(fun);
    } 

#if PKG_DEBUG
    fprintf(stderr, "read_ufunc_from_xml: returning %d\n", err);
#endif

    return err;
}

static int wordmatch (const char *s, const char *test)
{
    int n = strlen(test);
    
    return (!strncmp(s, test, n) && (s[n] == '\0' || isspace(s[n])));
}

void adjust_indent (const char *s, int *this_indent, int *next_indent)
{
    int ti = *next_indent;
    int ni = *next_indent;

    if (*s == '\0') {
	*this_indent = *next_indent;
	return;
    }

    if (!strncmp(s, "catch ", 6)) {
	s += 6;
	s += strspn(s, " ");
    }

    if (wordmatch(s, "loop")) {
	ni++;
    } else if (wordmatch(s, "if")) {
	ni++;
    } else if (wordmatch(s, "nls")) {
	ni++;
    } else if (wordmatch(s, "mle")) {
	ni++;
    } else if (wordmatch(s, "gmm")) {
	ni++;
    } else if (wordmatch(s, "mpi")) {
	ni++;
    } else if (wordmatch(s, "plot")) {
	ni++;
    } else if (wordmatch(s, "function")) {
	ni++;
    } else if (wordmatch(s, "restrict")) {
	ni++;
    } else if (wordmatch(s, "system")) {
	ni++;
    } else if (wordmatch(s, "foreign")) {
	ni++;
    } else if (wordmatch(s, "kalman")) {
	ni++;
    } else if (wordmatch(s, "end") ||
	       wordmatch(s, "endif") ||
	       wordmatch(s, "endloop")) {
	ti--;
	ni--;
    } else if (wordmatch(s, "else") ||
	       wordmatch(s, "elif")) {
	ni = ti;
	ti--;
    } 

    *this_indent = ti;
    *next_indent = ni;
}

/* ensure use of canonical forms "endif", "endloop" */

static void maybe_correct_line (char *line)
{
    char *p = strstr(line, "end if");

    if (p == NULL) {
	p = strstr(line, "end loop");
    }

    if (p != NULL && (p == line || *(p-1) == ' ')) {
	shift_string_left(p + 3, 1);
    }
}

#define parm_has_children(p) (p->descrip != NULL || p->nlabels > 0)

/* write out a single user-defined function as XML, according to
   gretlfunc.dtd */

static int write_function_xml (ufunc *fun, FILE *fp)
{
    int rtype = fun->rettype;
    int this_indent = 0;
    int next_indent = 0;
    int i, j;

    if (rtype == GRETL_TYPE_NONE) {
	rtype = GRETL_TYPE_VOID;
    }

    fprintf(fp, "<gretl-function name=\"%s\" type=\"%s\"", 
	    fun->name, gretl_type_get_name(rtype));

    if (function_is_private(fun)) {
	fputs(" private=\"1\"", fp);
    }    
    if (function_is_plugin(fun)) {
	fputs(" plugin-wrapper=\"1\"", fp);
    }
    if (function_is_noprint(fun)) {
	fputs(" no-print=\"1\"", fp);
    }
    if (function_is_menu_only(fun)) {
	fputs(" menu-only=\"1\"", fp);
    }    

    if (fun->pkg_role) {
	fprintf(fp, " pkg-role=\"%s\"", package_role_get_key(fun->pkg_role));
    }

    fputs(">\n", fp);

    if (fun->n_params > 0) {

	gretl_push_c_numeric_locale();

	fprintf(fp, " <params count=\"%d\">\n", fun->n_params);
	for (i=0; i<fun->n_params; i++) {
	    fn_param *param = &fun->params[i];

	    fprintf(fp, "  <param name=\"%s\" type=\"%s\"",
		    param->name, arg_type_xml_string(param->type));
	    if (!na(param->min)) {
		fprintf(fp, " min=\"%g\"", param->min);
	    }
	    if (!na(param->max)) {
		fprintf(fp, " max=\"%g\"", param->max);
	    }
	    if (!default_unset(param)) {
		if (na(param->deflt)) {
		    fputs(" default=\"NA\"", fp);
		} else {
		    fprintf(fp, " default=\"%g\"", param->deflt);
		}
	    }
	    if (!na(param->step)) {
		fprintf(fp, " step=\"%g\"", param->step);
	    }
	    if (param->flags & ARG_OPTIONAL) {
		fputs(" optional=\"true\"", fp);
	    }
	    if (param->flags & ARG_CONST) {
		fputs(" const=\"true\"", fp);
	    }
	    if (parm_has_children(param)) {
		fputs(">\n", fp); /* terminate opening tag */
		if (param->descrip != NULL) {
		    gretl_xml_put_tagged_string("description", 
						param->descrip,
						fp);
		} 
		if (param->nlabels > 0) {
		    gretl_xml_put_strings_array_quoted("labels", 
						       (const char **) param->labels, 
						       param->nlabels, fp);
		}
		fputs("  </param>\n", fp);
	    } else {
		fputs("/>\n", fp); /* terminate opening tag */
	    }
	}
	fputs(" </params>\n", fp);

	gretl_pop_c_numeric_locale();
    }

    fputs("<code>", fp);

    for (i=0; i<fun->n_lines; i++) {
	adjust_indent(fun->lines[i].s, &this_indent, &next_indent);
	for (j=0; j<this_indent; j++) {
	    fputs("  ", fp);
	}
	maybe_correct_line(fun->lines[i].s);
	gretl_xml_put_string(fun->lines[i].s, fp);
	fputc('\n', fp);
    }

    fputs("</code>\n", fp);

    fputs("</gretl-function>\n", fp);
    
    return 0;
}

/* script-style output */

static void print_function_start (ufunc *fun, PRN *prn)
{
    const char *s;
    int i, pos = 0;

    if (fun->rettype == GRETL_TYPE_NONE) {
	pos += pprintf(prn, "function void %s ", fun->name);
    } else {
	const char *typestr = gretl_type_get_name(fun->rettype);

	pos += pprintf(prn, "function %s %s ", typestr, fun->name);
    }

    gretl_push_c_numeric_locale();

    if (fun->n_params == 0) {
	pputs(prn, "(void)");
    } else {
	pos += pputc(prn, '(');
    }

    for (i=0; i<fun->n_params; i++) {
	fn_param *fp = &fun->params[i];

	if (fp->flags & ARG_CONST) {
	    pputs(prn, "const ");
	}
	s = gretl_type_get_name(fp->type);
	if (s[strlen(s) - 1] == '*') {
	    pprintf(prn, "%s%s", s, fp->name);
	} else {
	    pprintf(prn, "%s %s", s, fp->name);
	}
	if (fp->type == GRETL_TYPE_BOOL) {
	    if (!default_unset(fp) && !na(fp->deflt)) {
		pprintf(prn, "[%g]", fp->deflt);
	    }
	} else if (gretl_scalar_type(fp->type)) {
	    print_min_max_deflt(fp, prn);
	} else if (arg_may_be_optional(fp->type)) {
	    print_opt_flags(&fun->params[i], prn);
	}
	print_param_description(fp, prn);
	if (fp->nlabels > 0) {
	    print_param_labels(fp, prn);
	}
	if (i == fun->n_params - 1) {
	    pputc(prn, ')');
	} else {
	    pputs(prn, ",\n");
	    bufspace(pos, prn);
	}
    }

    pputc(prn, '\n');

    gretl_pop_c_numeric_locale();
}

/**
 * gretl_function_print_code:
 * @u: pointer to user-function.
 * @tabwidth: number of spaces per "tab" (logical indent).
 * @prn: printing struct.
 *
 * Prints out function @fun to @prn, script-style.
 *
 * Returns: 0 on success, non-zero if @fun is %NULL.
 */

int gretl_function_print_code (ufunc *u, int tabwidth, PRN *prn)
{
    int this_indent = 0;
    int next_indent = 0;
    int i, j;

    if (u == NULL) {
	return E_DATA;
    }

    if (tabwidth == 0) {
	tabwidth = 2;
    }
   
    print_function_start(u, prn);

    for (i=0; i<u->n_lines; i++) {
	adjust_indent(u->lines[i].s, &this_indent, &next_indent);
	for (j=0; j<=this_indent; j++) {
	    bufspace(tabwidth, prn);
	}
	pputs(prn, u->lines[i].s);
	if (i < u->n_lines - 1) {
	    if (u->lines[i+1].idx > u->lines[i].idx + 1) {
		pputc(prn, '\n');
	    }
	}
	pputc(prn, '\n');
    }

    pputs(prn, "end function\n");

    return 0;
}

char **gretl_function_retrieve_code (ufunc *u, int *nlines)
{
    char **S = NULL;
    int i, j = 0;

    for (i=0; i<u->n_lines; i++) {
	if (!u->lines[i].ignore) {
	    j++;
	}
    }

    if (j > 0) {
	S = strings_array_new(j);
    }

    if (S != NULL) {
	*nlines = j;
	j = 0;
	for (i=0; i<u->n_lines; i++) {
	    if (!u->lines[i].ignore) {
		S[j++] = u->lines[i].s;
	    }
	}	
    }

    return S;
}

/* construct a name for @pkg based on its filename member:
   take the basename and knock off ".gfn"
*/

static void name_package_from_filename (fnpkg *pkg)
{
    char *p = strrchr(pkg->fname, SLASH);
    int n;

    if (p != NULL) {
	p++;
    } else {
	p = pkg->fname;
    }

    n = strlen(p);
    if (has_suffix(p, ".gfn")) {
	n -= 4;
    }

    if (n > FN_NAMELEN - 1) {
	n = FN_NAMELEN - 1;
    }

    *pkg->name = '\0';
    strncat(pkg->name, p, n);

#if PKG_DEBUG
    fprintf(stderr, "filename '%s' -> pkgname '%s'\n",
	    pkg->fname, pkg->name);
#endif
}

/* name lookup for functions to be connected to a package,
   allowing for the possibility that they're already
   connected */

static ufunc *get_uf_array_member (const char *name, fnpkg *pkg)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == pkg || ufuns[i]->pkg == NULL) {
	    if (!strcmp(name, ufuns[i]->name)) {
		return ufuns[i];
	    }
	}
    }   

    return NULL;
}

static void check_special_comments (ufunc *fun)
{
    int i;

    for (i=0; i<fun->n_lines; i++) {
	if (strstr(fun->lines[i].s, "## plugin-wrapper ##")) {
	    fun->flags |= UFUN_PLUGIN;
	} else if (strstr(fun->lines[i].s, "## no-print ##")) {
	    fun->flags |= UFUN_NOPRINT;
	} else if (strstr(fun->lines[i].s, "## menu-only ##")) {
	    fun->flags |= UFUN_MENU_ONLY;
	}
    }
}

/* before revising the function members of an existing
   package, detach all its current functions
*/

static void package_disconnect_funcs (fnpkg *pkg)
{
    int i;

    if (pkg->pub != NULL) {
	for (i=0; i<pkg->n_pub; i++) {
	    pkg->pub[i]->pkg = NULL;
	}
	free(pkg->pub);
	pkg->pub = NULL;
	pkg->n_pub = 0;
    }

    if (pkg->priv != NULL) {
	for (i=0; i<pkg->n_priv; i++) {
	    pkg->priv[i]->pkg = NULL;
	    set_function_private(pkg->priv[i], FALSE);
	}
	free(pkg->priv);
	pkg->priv = NULL;
	pkg->n_priv = 0;
    }
}

/* Given an array of @n function names in @names, set up a
   corresponding array of pointers to the named functions in @pkg,
   Flag an error if any of the function names are bad, or if
   allocation fails.
*/

static int set_uf_array_from_names (fnpkg *pkg, char **names, 
				    int n, int priv)
{
    ufunc **uf = NULL;
    ufunc *fun;
    int i;

    /* check the supplied names */
    for (i=0; i<n; i++) {
	if (get_uf_array_member(names[i], NULL) == NULL) {
	    fprintf(stderr, "%s: function not found!\n", names[i]);
	    return E_DATA;
	}
    }

    /* allocate storage */
    if (n > 0) {
	uf = malloc(n * sizeof *uf);
	if (uf == NULL) {
	    return E_ALLOC;
	}
    }

    /* connect the specified functions */
    for (i=0; i<n; i++) {
	fun = get_uf_array_member(names[i], NULL);
	fun->pkg = pkg;
	set_function_private(fun, priv);
	check_special_comments(fun);
	uf[i] = fun;
    }

    /* set up the package info */
    if (priv) {
	pkg->priv = uf;
	pkg->n_priv = n;
    } else {
	pkg->pub = uf;
	pkg->n_pub = n;
    }

    return 0;
}

/**
 * function_package_connect_funcs:
 * @pkg: function package.
 * @pubnames: array of names of public functions.
 * @n_pub: number of strings in @pubnames.
 * @privnames: array of names of private functions (or %NULL).
 * @n_priv: number of strings in @privnames (may be 0).
 *
 * Looks up the functions named in @pubnames and @privnames
 * and adds pointers to these functions to @pkg, hence marking
 * the functions as belonging to @pkg.
 *
 * Returns: 0 on success, non-zero on error.
 */

int function_package_connect_funcs (fnpkg *pkg, 
				    char **pubnames, int n_pub,
				    char **privnames, int n_priv) 
{
    int err;

    /* clear the decks first */
    package_disconnect_funcs(pkg);  

    err = set_uf_array_from_names(pkg, pubnames, n_pub, 0);

    if (!err) {
	err = set_uf_array_from_names(pkg, privnames, n_priv, 1);
    }

    return err;
}

/**
 * function_package_new:
 * @fname: filename for package.
 * @pubnames: array of names of public functions.
 * @n_pub: number of strings in @pubnames.
 * @privnames: array of names of private functions (or %NULL).
 * @n_priv: number of strings in @privnames (may be 0).
 * @err: location to receive error code.
 *
 * Allocates a new package with filename-member @fname, including
 * the public and private functions named in @pubnames and
 * @privnames.  Note that this function does not cause the
 * package to be written to file; for that, see
 * write_function_package().
 *
 * Returns: pointer to package on success, %NULL on error.
 */

fnpkg *function_package_new (const char *fname, 
			     char **pubnames, int n_pub,
			     char **privnames, int n_priv, 
			     int *err)
{
    fnpkg *pkg = NULL;

    if (n_pub <= 0) {
	/* we need at least one public function */
	*err = E_DATA;
    } else {
	pkg = function_package_alloc(fname);
	if (pkg == NULL) {
	    *err = E_ALLOC;
	}
    } 

    if (*err) {
	return NULL;
    }

    name_package_from_filename(pkg);

    *err = function_package_connect_funcs(pkg, pubnames, n_pub,
					  privnames, n_priv);

    if (!*err) {
	*err = function_package_record(pkg);
    }

    if (*err) {
	/* note: this does not free the packaged functions, if any */
	function_package_free(pkg);
	pkg = NULL;
    }

    return pkg;
}

static char *trim_text (char *s)
{
    while (isspace(*s)) s++;
    return tailstrip(s);
}

static int package_write_translatable_strings (fnpkg *pkg, PRN *prn)
{
    FILE *fp;
    gchar *trname;
    char **S = NULL;
    int i, n = 0;

    trname = g_strdup_printf("%s-i18n.c", pkg->name);

    fp = gretl_fopen(trname, "w");
    if (fp == NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s"), trname);
	g_free(trname);
	return E_FOPEN;
    }

    if (pkg->pub != NULL) {
	int j, k;

	for (i=0; i<pkg->n_pub; i++) {
	    ufunc *fun = pkg->pub[i];

	    for (j=0; j<fun->n_params; j++) {
		fn_param *param = &fun->params[j];

		if (param->descrip != NULL) {
		    strings_array_add(&S, &n, param->descrip);
		}
		for (k=0; k<param->nlabels; k++) {
		    strings_array_add(&S, &n, param->labels[k]);
		}
	    }
	}
    }

    if (pkg->label != NULL || S != NULL) {
	fprintf(fp, "const char *%s_translations[] = {\n", pkg->name);
	if (pkg->label != NULL) {
	    fprintf(fp, "    N_(\"%s\"),\n", pkg->label);
	}	
	if (S != NULL) {
	    strings_array_sort(&S, &n, OPT_U);
	    for (i=0; i<n; i++) {
		fprintf(fp, "    N_(\"%s\")", S[i]);
		if (i < n-1) {
		    fputc(',', fp);
		}
		fputc('\n', fp);
	    }
	    strings_array_free(S, n);
	}
	fputs("};\n", fp);
    }

    fclose(fp);

    pprintf(prn, "Wrote translations file %s\n", trname);
    g_free(trname);

    return 0;
}

static int package_write_index (fnpkg *pkg, PRN *prn)
{
    gchar *idxname;
    FILE *fp;

    idxname = g_strdup_printf("%s.xml", pkg->name);

    fp = gretl_fopen(idxname, "w");
    if (fp == NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s"), idxname);
	g_free(idxname);
	return E_FOPEN;
    }

    gretl_xml_header(fp);

    fputs("<gretl-addon", fp);

    fprintf(fp, " name=\"%s\"", pkg->name);

    if (pkg->dreq == FN_NEEDS_TS) {
	fprintf(fp, " %s=\"true\"", NEEDS_TS);
    } else if (pkg->dreq == FN_NEEDS_QM) {
	fprintf(fp, " %s=\"true\"", NEEDS_QM);
    } else if (pkg->dreq == FN_NEEDS_PANEL) {
	fprintf(fp, " %s=\"true\"", NEEDS_PANEL);
    } else if (pkg->dreq == FN_NODATA_OK) {
	fprintf(fp, " %s=\"true\"", NO_DATA_OK);
    }

    if (pkg->modelreq > 0) {
	fprintf(fp, " model-requirement=\"%s\"", 
		gretl_command_word(pkg->modelreq));
    }

    if (pkg->minver > 0) {
	char vstr[8];

	fprintf(fp, " minver=\"%s\"", gretl_version_string(vstr, pkg->minver));
    }

    if (pkg->uses_subdir) {
	fprintf(fp, " lives-in-subdir=\"true\"");
    }    

    fputs(">\n", fp);

    gretl_xml_put_tagged_string("author",  pkg->author, fp);
    gretl_xml_put_tagged_string("version", pkg->version, fp);
    gretl_xml_put_tagged_string("date",    pkg->date, fp);
    gretl_xml_put_tagged_string("description", pkg->descrip, fp);

    if (pkg->tags != NULL) {
	gretl_xml_put_tagged_string("tags", pkg->tags, fp);
    }    

    if (pkg->label != NULL) {
	gretl_xml_put_tagged_string("label", pkg->label, fp);
    }

    if (pkg->mpath != NULL) {
	gretl_xml_put_tagged_string("menu-attachment", pkg->mpath, fp);
    }

    fputs("</gretl-addon>\n", fp);

    fclose(fp);

    pprintf(prn, "Wrote index file %s\n", idxname);
    g_free(idxname);

    return 0;
}

/* Write out a function package as XML. If @fp is NULL, then we write
   this as a package file in its own right, using the filename given
   in the package (the 'standalone' case), otherwise we're writing it
   as a component of the functions listing in a gretl session file,
   which already has an associated open stream.  
*/

static int real_write_function_package (fnpkg *pkg, FILE *fp)
{
    int standalone = (fp == NULL);
    int i, err = 0;

    if (standalone) {
	fp = gretl_fopen(pkg->fname, "w");
	if (fp == NULL) {
	    gretl_errmsg_sprintf(_("Couldn't open %s"), pkg->fname);
	    return E_FOPEN;
	} else {
	    gretl_xml_header(fp);
	    fputs("<gretl-functions>\n", fp);
	}
    }

    fputs("<gretl-function-package", fp);

    name_package_from_filename(pkg); /* ? */

    fprintf(fp, " name=\"%s\"", pkg->name);
    fprintf(fp, " ID=\"%d\"", pkg->ID);

    if (pkg->dreq == FN_NEEDS_TS) {
	fprintf(fp, " %s=\"true\"", NEEDS_TS);
    } else if (pkg->dreq == FN_NEEDS_QM) {
	fprintf(fp, " %s=\"true\"", NEEDS_QM);
    } else if (pkg->dreq == FN_NEEDS_PANEL) {
	fprintf(fp, " %s=\"true\"", NEEDS_PANEL);
    } else if (pkg->dreq == FN_NODATA_OK) {
	fprintf(fp, " %s=\"true\"", NO_DATA_OK);
    }

    if (pkg->modelreq > 0) {
	fprintf(fp, " model-requirement=\"%s\"", 
		gretl_command_word(pkg->modelreq));
    }

    if (pkg->minver > 0) {
	char vstr[8];

	fprintf(fp, " minver=\"%s\"", gretl_version_string(vstr, pkg->minver));
    }

    if (pkg->uses_subdir) {
	fprintf(fp, " lives-in-subdir=\"true\"");
    }    

    fputs(">\n", fp);

    if (pkg->email != NULL && *pkg->email != '\0') {
	gretl_xml_put_tagged_string_plus("author", pkg->author,
					 "email", pkg->email,
					 fp);
    } else {
	gretl_xml_put_tagged_string("author", pkg->author, fp);
    }
    gretl_xml_put_tagged_string("version", pkg->version, fp);
    gretl_xml_put_tagged_string("date",    pkg->date, fp);
    gretl_xml_put_tagged_string("description", pkg->descrip, fp);

    if (pkg->tags != NULL) {
	gretl_xml_put_tagged_string("tags", pkg->tags, fp);
    }    

    if (pkg->label != NULL) {
	gretl_xml_put_tagged_string("label", pkg->label, fp);
    }

    if (pkg->mpath != NULL) {
	gretl_xml_put_tagged_string("menu-attachment", pkg->mpath, fp);
    }

    if (pkg->help != NULL) {
	if (pkg->help_fname != NULL) {
	    fprintf(fp, "<help filename=\"%s\">\n", pkg->help_fname);
	} else {
	    fputs("<help>\n", fp);
	}
	gretl_xml_put_string(trim_text(pkg->help), fp);
	fputs("\n</help>\n", fp);
    } 

    if (pkg->gui_help != NULL) {
	if (pkg->gui_help_fname != NULL) {
	    fprintf(fp, "<gui-help filename=\"%s\">\n",
		    pkg->gui_help_fname);
	} else {
	    fputs("<gui-help>\n", fp);
	}	
	gretl_xml_put_string(trim_text(pkg->gui_help), fp);
	fputs("\n</gui-help>\n", fp);
    }

    if (pkg->datafiles != NULL) {
	gretl_xml_put_strings_array("data-files",
				    (const char **) pkg->datafiles,
				    pkg->n_files, fp);
    }    

    if (pkg->pub != NULL) {
	for (i=0; i<pkg->n_pub; i++) {
	    write_function_xml(pkg->pub[i], fp);
	}
    }    

    if (pkg->priv != NULL) {
	for (i=0; i<pkg->n_priv; i++) {
	    write_function_xml(pkg->priv[i], fp);
	}
    }

    if (pkg->sample != NULL) {
	if (pkg->sample_fname != NULL) {
	    fprintf(fp, "<sample-script filename=\"%s\">\n",
		    pkg->sample_fname);
	} else {
	    fputs("<sample-script>\n", fp);
	}
	gretl_xml_put_string(trim_text(pkg->sample), fp);
	fputs("\n</sample-script>\n", fp);	
    }

    fputs("</gretl-function-package>\n", fp);

    if (standalone) {
	fputs("</gretl-functions>\n", fp);
	fclose(fp);
    }

    return err;
}

/**
 * function_package_write_file:
 * @pkg: function package.
 *
 * Write out @pkg as an XML file, using the filename
 * recorded in the package.
 *
 * Returns: 0 on success, non-zero on error.
 */

int function_package_write_file (fnpkg *pkg)
{
    return real_write_function_package(pkg, NULL);
}

/* below: apparatus for constructing and saving a gfn function
   package from the command line
*/

static int is_unclaimed (const char *s, char **S, int n)
{
    int i;

    for (i=0; i<n; i++) {
	if (strcmp(s, S[i]) == 0) {
	    return 0;
	}
    }

    return 1;
}

static gchar *pkg_aux_content (const char *fname, int *err)
{
    gchar *ret = NULL;
    GError *gerr = NULL;
    gsize len = 0;

    g_file_get_contents(fname, &ret, &len, &gerr);

    if (gerr != NULL) {
	gretl_errmsg_set(gerr->message);
	g_error_free(gerr);
	*err = E_FOPEN;
    }

    return ret;
}

static int pkg_set_dreq (fnpkg *pkg, const char *s)
{
    int err = 0;

    if (!strcmp(s, NEEDS_TS)) {
	pkg->dreq = FN_NEEDS_TS;
    } else if (!strcmp(s, NEEDS_QM)) {
	pkg->dreq = FN_NEEDS_QM;
    } else if (!strcmp(s, NEEDS_PANEL)) {
	pkg->dreq = FN_NEEDS_PANEL;
    } else if (!strcmp(s, NO_DATA_OK)) {
	pkg->dreq = FN_NODATA_OK;
    } else {
	err = E_PARSE;
    }

    return err;
}

static int pkg_set_modelreq (fnpkg *pkg, const char *s)
{
    int ci = gretl_command_number(s);

    if (ci > 0) {
	pkg->modelreq = ci;
	return 0;
    } else {
	return E_PARSE;
    }
}

static int pkg_set_datafiles (fnpkg *pkg, const char *s)
{
    int err = 0;

    if (string_is_blank(s)) {
	err = E_DATA;
    } else {
	int n = 0;
	
	pkg->datafiles = gretl_string_split(s, &n, NULL);
	if (pkg->datafiles == NULL) {
	    pkg->n_files = 0;
	    err = E_ALLOC;
	} else {
	    pkg->n_files = n;
	}
    }

    if (!err) {
	pkg->uses_subdir = 1;
    }

    return err;
}

static int pkg_boolean_from_string (const char *s)
{
    if (!strcmp(s, "true")) {
	return 1;
    } else {
	return 0;
    }
}

static int pkg_remove_role (fnpkg *pkg, int role)
{
    ufunc *u;
    int i;

    for (i=0; i<pkg->n_priv; i++) {
	u = pkg->priv[i];
	if (u->pkg_role == role) {
	    u->pkg_role = UFUN_ROLE_NONE;
	    return 0;
	}
    }
    for (i=0; i<pkg->n_pub; i++) {
	u = pkg->pub[i];
	if (u->pkg_role == role) {
	    u->pkg_role = UFUN_ROLE_NONE;
	    return 0;
	}
    }

    return E_DATA;
}

int function_set_package_role (const char *name, fnpkg *pkg,
			       const char *attr, PRN *prn)
{
    ufunc *u = NULL;
    int role = pkg_key_get_role(attr);
    int i, j, err = 0;

    if (name == NULL) {
	/* removing a role */
	pkg_remove_role(pkg, role);
	return 0;
    }

    if (role == UFUN_ROLE_NONE) {
	for (i=0; i<pkg->n_priv; i++) {
	    if (!strcmp(name, pkg->priv[i]->name)) {
		u = pkg->priv[i];
		u->pkg_role = role;
		return 0;
	    }
	}
	for (i=0; i<pkg->n_pub; i++) {
	    if (!strcmp(name, pkg->pub[i]->name)) {
		u = pkg->pub[i];
		u->pkg_role = role;
		return 0;
	    }
	}
	return E_DATA;
    }

    /* check that the function in question satisfies the
       requirements for its role, and if so, hook it up 
    */

    if (role == UFUN_GUI_PRECHECK) {
	/* the pre-checker must be a private function */
	for (i=0; i<pkg->n_priv; i++) {
	    if (!strcmp(name, pkg->priv[i]->name)) {
		u = pkg->priv[i];
		if (u->rettype != GRETL_TYPE_DOUBLE) {
		    pprintf(prn, "%s: must return a scalar\n", attr);
		    err = E_TYPES;
		} else if (u->n_params > 0) {
		    pprintf(prn, "%s: no parameters are allowed\n", attr);
		    err = E_TYPES;
		}
		if (!err) {
		    u->pkg_role = role;
		}		
		return err; /* found */
	    }
	}
	/* not found */
	pprintf(prn, "%s: %s: no such private function\n", attr, name);
	return E_DATA;
    }

    for (i=0; i<pkg->n_pub; i++) {
	/* all other special-role functions must be public */
	if (!strcmp(name, pkg->pub[i]->name)) {
	    u = pkg->pub[i];
	    if (role == UFUN_GUI_MAIN) {
#if 1
		; /* OK, type does not matter */
#else	
		if (u->rettype != GRETL_TYPE_BUNDLE && u->rettype != GRETL_TYPE_VOID) {
		    pprintf(prn, "%s: must return a bundle, or nothing\n", attr);
		    err = E_TYPES;
		}
#endif
	    } else {
		/* bundle-print, bundle-plot, etc. */
		if (u->n_params == 0) {
		    pprintf(prn, "%s: must take a %s argument\n", attr,
			    gretl_type_get_name(GRETL_TYPE_BUNDLE_REF));
		    err = E_TYPES;
		}
		for (j=0; j<u->n_params && !err; j++) {
		    if (j == 0 && u->params[j].type != GRETL_TYPE_BUNDLE_REF) {
			pprintf(prn, "%s: first param type must be %s\n",
				attr, gretl_type_get_name(GRETL_TYPE_BUNDLE_REF));
			err = E_TYPES;
		    } else if (j == 1 && u->params[j].type != GRETL_TYPE_INT) {
			pprintf(prn, "%s: second param type must be %s\n",
				attr, gretl_type_get_name(GRETL_TYPE_INT));
			err = E_TYPES;
		    } else if (j > 1 && !fn_param_optional(u, j) &&
			       na(fn_param_default(u, j))) {
			pprintf(prn, "%s: params beyond the second must be optional\n",
				attr);
			err = E_TYPES;
		    }
		}
	    }
	    if (!err) {
		u->pkg_role = role;
	    }
	    return err;
	}
    }

    pprintf(prn, "%s: %s: no such public function\n", attr, name);

    return E_DATA;
}

/* called from the GUI package editor to check whether
   a given function of name @name can be shown as a 
   candidate for a specified GUI-special @role
*/

int function_ok_for_package_role (const char *name,
				  int role)
{
    ufunc *u = NULL;
    int i, err = 0;

    if (name == NULL || role == UFUN_ROLE_NONE) {
	return 0;
    }
    
    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(name, ufuns[i]->name)) { 
	    u = ufuns[i];
	    break;
	}
    }

    if (u == NULL) {
	return 0;
    }

    if (role == UFUN_GUI_PRECHECK) {
	if (u->rettype != GRETL_TYPE_DOUBLE) {
	    err = E_TYPES;
	} else if (u->n_params > 0) {
	    err = E_TYPES;
	}
	return !err; /* found */
    }

    if (role == UFUN_GUI_MAIN) {
#if 1
	; /* OK, we don't mind what type it is */
#else
	if (u->rettype != GRETL_TYPE_BUNDLE && u->rettype != GRETL_TYPE_VOID) {
	    err = E_TYPES;
	}
#endif
    } else {
	/* bundle-print, bundle-plot, etc. */
	if (u->n_params == 0) {
	    err = E_TYPES;
	}
	for (i=0; i<u->n_params && !err; i++) {
	    if (i == 0 && u->params[i].type != GRETL_TYPE_BUNDLE_REF) {
		err = E_TYPES;
	    } else if (i == 1 && u->params[i].type != GRETL_TYPE_INT) {
		err = E_TYPES;
	    } else if (i > 1 && !fn_param_optional(u, i) &&
		       na(fn_param_default(u, i))) {
		err = E_TYPES;
	    }
	}
    }

    return !err;
}

static int pkg_set_funcs_attribute (fnpkg *pkg, const char *s,
				    int flag)
{
    char **S;
    int ns = 0, err = 0;

    S = gretl_string_split(s, &ns, NULL);
    if (ns == 0) {
	err = E_DATA;
    } else {
	int i, j, match;

	for (i=0; i<ns && !err; i++) {
	    match = 0;
	    for (j=0; j<pkg->n_pub; j++) {
		if (!strcmp(S[i], pkg->pub[j]->name)) {
		    pkg->pub[j]->flags |= flag;
		    match = 1;
		    break;
		}
	    }
	    if (!match) {
		err = E_DATA;
	    }
	}

	strings_array_free(S, ns);
    }

    return err;
}

static void function_package_set_auxfile (fnpkg *pkg,
					  const char *id,
					  const char *fname)
{
    gchar *test = NULL;
    
    /* Maybe set the source filename for an element
       of a function package read from file, but only
       if it's not the standard, default filename for
       the element in question.
    */

    if (!strcmp(id, "help-fname")) {
	test = g_strdup_printf("%s_help.txt", pkg->name);
    } else if (!strcmp(id, "gui-help-fname")) {
	test = g_strdup_printf("%s_gui_help.txt", pkg->name);
    } else if (!strcmp(id, "sample-fname")) {
	test = g_strdup_printf("%s_sample.inp", pkg->name);
    }

    if (test != NULL) {
	if (strcmp(fname, test)) {
	   function_package_set_properties(pkg, id, fname, NULL);
	}
	g_free(test);
    }
}

/* Having assembled and checked the function-listing for a new
   package, now retrieve the additional information from the
   spec file (named by @fname, opened as @fp).
*/

static int new_package_info_from_spec (fnpkg *pkg, const char *fname,
				       FILE *fp, gretlopt opt,
				       PRN *prn)
{
    const char *okstr, *failed;
    int quiet = (opt & OPT_Q);
    char *p, line[1024];
    gchar *tmp;
    int dfd = -1;
    int got = 0;
    int err = 0;

    if (opt & OPT_G) {
	/* GUI use */
	okstr = "<@ok>\n";
	failed = "<@fail>\n";
    } else {
	okstr = "OK\n";
	failed = "failed\n";
    }

    if (strrchr(fname, SLASH) != NULL) {
	char dirname[FILENAME_MAX];

	strcpy(dirname, fname);
	p = strrchr(dirname, SLASH);
	*p = '\0';
	dfd = gretl_open(".", O_RDONLY, 0);
	err = gretl_chdir(dirname);
    }

    while (fgets(line, sizeof line, fp) && !err) {
	if (*line == '#' || string_is_blank(line)) {
	    continue;
	}
	tailstrip(line);
	p = strchr(line, '=');
	if (p == NULL) {
	    continue;
	} else {
	    p++;
	    p += strspn(p, " ");
	    if (!strncmp(line, "author", 6)) {
		err = function_package_set_properties(pkg, "author", p, NULL);
		if (!err) got++;
	    } else if (!strncmp(line, "email", 5)) {
		err = function_package_set_properties(pkg, "email", p, NULL);
	    } else if (!strncmp(line, "version", 7)) {
		err = function_package_set_properties(pkg, "version", p, NULL);
		if (!err) got++;
	    } else if (!strncmp(line, "date", 4)) {
		err = function_package_set_properties(pkg, "date", p, NULL);
		if (!err) got++;
	    } else if (!strncmp(line, "description", 11)) {
		err = function_package_set_properties(pkg, "description", p, NULL);
		if (!err) got++;
	    } else if (!strncmp(line, "tags", 4)) {
		err = function_package_set_properties(pkg, "tags", p, NULL);
	    } else if (!strncmp(line, "label", 5)) {
		err = function_package_set_properties(pkg, "label", p, NULL);
	    } else if (!strncmp(line, "menu-attachment", 15)) {
		err = function_package_set_properties(pkg, "menu-attachment", p, NULL);
	    } else if (!strncmp(line, "help", 4)) {
		int pdfdoc = 0;
		
		if (has_suffix(p, ".pdf")) {
		    if (!quiet) {
			pprintf(prn, "Recording help reference %s\n", p);
		    }
		    tmp = g_strdup_printf("pdfdoc:%s", p);
		    pdfdoc = 1;
		} else {
		    if (!quiet) {
			pprintf(prn, "Looking for help text in %s... ", p);
		    }
		    tmp = pkg_aux_content(p, &err);
		    if (err) {
			pputs(prn, failed);
		    } else {
			pputs(prn, okstr);
		    }
		}
		if (!err) {
		    err = function_package_set_properties(pkg, "help", tmp, NULL);
		    if (!err) {
			got++;
			if (!pdfdoc) {
			    function_package_set_auxfile(pkg, "help-fname", p);
			}
		    }
		    g_free(tmp);
		}
	    } else if (!strncmp(line, "gui-help", 8)) {
		if (!quiet) {
		    pprintf(prn, "Looking for GUI help text in %s... ", p);
		}
		tmp = pkg_aux_content(p, &err);
		if (err) {
		    pputs(prn, failed);
		} else {
		    pputs(prn, okstr);
		    err = function_package_set_properties(pkg, "gui-help", tmp, NULL);
		    if (!err) {
			function_package_set_auxfile(pkg, "gui-help-fname", p);
		    }
		    g_free(tmp);
		}
	    } else if (!strncmp(line, "sample-script", 13)) {
		if (!quiet) {
		    pprintf(prn, "Looking for sample script in %s... ", p);
		}
		tmp = pkg_aux_content(p, &err);
		if (err) {
		    pputs(prn, failed);
		} else {
		    pputs(prn, okstr);
		    err = function_package_set_properties(pkg, "sample-script", tmp, NULL);
		    if (!err) {
			got++;
			function_package_set_auxfile(pkg, "sample-fname", p);
		    }
		    g_free(tmp);
		}
	    } else if (!strncmp(line, "data-files", 10)) {
		if (!quiet) {
		    pprintf(prn, "Recording data-file list: %s\n", p);
		}
		err = pkg_set_datafiles(pkg, p);
	    } else if (!strncmp(line, "data-requirement", 16)) {
		err = pkg_set_dreq(pkg, p);
	    } else if (!strncmp(line, "model-requirement", 17)) {
		err = pkg_set_modelreq(pkg, p);
	    } else if (!strncmp(line, "min-version", 11)) {
		pkg->minver = gretl_version_number(p);
		got++;
	    } else if (!strncmp(line, "lives-in-subdir", 15)) {
		pkg->uses_subdir = pkg_boolean_from_string(p);
	    } else if (!strncmp(line, "no-print", 8)) {
		err = pkg_set_funcs_attribute(pkg, p, UFUN_NOPRINT);
	    } else if (!strncmp(line, "menu-only", 9)) {
		err = pkg_set_funcs_attribute(pkg, p, UFUN_MENU_ONLY);
	    } else {
		const char *key;
		int i;

		for (i=0; pkg_lookups[i].key != NULL; i++) {
		    key = pkg_lookups[i].key;
		    if (!strncmp(line, key, strlen(key))) {
			err = function_set_package_role(p, pkg, key, prn);
			if (!err && !quiet) {
			    pprintf(prn, "%s function is %s, %s", key, p, okstr);
			}
			break;
		    }
		}
	    }
	}
    }

    if (dfd >= 0) {
	/* go back where we came from */
	gretl_fchdir(dfd);
    }

    if (!err && got < 7) {
	gretl_errmsg_set("Some required information was missing");
	err = E_DATA;
    }

    return err;
}

static fnpkg *new_pkg_from_spec_file (const char *gfnname, gretlopt opt,
				      PRN *prn, int *err)
{
    fnpkg *pkg = NULL;
    char fname[FILENAME_MAX];
    char line[4096], cont[1024];
    int quiet = (opt & OPT_Q);
    FILE *fp;

    if (!has_suffix(gfnname, ".gfn")) {
	gretl_errmsg_set("Output must have extension \".gfn\"");
	*err = E_DATA;
	return NULL;
    }

    switch_ext(fname, gfnname, "spec");
    fp = gretl_fopen(fname, "r");

    if (fp == NULL) {
	*err = E_FOPEN;
    } else {
	char **pubnames = NULL;
	char **privnames = NULL;
	int npub = 0, npriv = 0;
	ufunc *uf;
	int i;

	if (!quiet) {
	    pprintf(prn, "Found spec file '%s'\n", fname);
	}

	/* first pass: gather names of public functions */

	while (fgets(line, sizeof line, fp) && !*err) {
	    if (!strncmp(line, "public =", 8)) {
		while (ends_with_backslash(line)) {
		    gretl_charsub(line, '\\', '\0');
		    *cont = '\0';
		    if (fgets(cont, sizeof cont, fp)) {
			strcat(line, cont);
		    } else {
			*err = E_DATA;
		    }
		}
		if (!*err) {
		    tailstrip(line);
		    pubnames = gretl_string_split(line + 8, &npub, NULL);
		    if (npub == 0) {
			*err = E_DATA;
		    }
		}
	    } 
	}

	if (!*err) {
	    if (!quiet) {
		pprintf(prn, "number of public interfaces = %d\n", npub);
	    }
	    for (i=0; i<npub && !*err; i++) {
		uf = get_user_function_by_name(pubnames[i]);
		if (!quiet) {
		    pprintf(prn, " %s", pubnames[i]);
		}
		if (uf == NULL) {
		    if (!quiet) {
			pputs(prn, ": *** not found");
		    }
		    gretl_errmsg_sprintf("'%s': no such function", pubnames[i]);
		    *err = E_DATA;
		}
		if (!quiet) {
		    pputc(prn, '\n');
		}
	    }
	}

	/* note: any other functions currently loaded are assumed to be
	   private functions for this package */

	if (!*err) {
	    npriv = n_free_functions() - npub;
	    if (npriv < 0) {
		npriv = 0;
	    }
	}	

	if (!*err && npriv > 0) {
	    npriv = 0;
	    for (i=0; i<n_ufuns && !*err; i++) {
		if (ufuns[i]->pkg == NULL && is_unclaimed(ufuns[i]->name, pubnames, npub)) {
		    *err = strings_array_add(&privnames, &npriv, ufuns[i]->name);
		}
	    }
	}

	if (!quiet && !*err && npriv > 0) {
	    pprintf(prn, "number of private functions = %d\n", npriv);
	    for (i=0; i<npriv; i++) {
		pprintf(prn, " %s\n", privnames[i]);
	    }
	}

	if (!*err) {
	    pkg = function_package_new(gfnname, pubnames, npub,
				       privnames, npriv, err);
	}

	strings_array_free(pubnames, npub);
	strings_array_free(privnames, npriv);

	if (!*err) {
	    rewind(fp);
	    *err = new_package_info_from_spec(pkg, fname, fp, opt, prn);
	}

	if (*err && pkg != NULL) {
	    real_function_package_unload(pkg, 0);
	    pkg = NULL;
	}

	fclose(fp);
    }

    return pkg;
}

static int cli_validate_package_file (const char *fname,
				      gretlopt opt, PRN *prn)
{
    char dtdname[FILENAME_MAX];
    xmlDocPtr doc = NULL;
    xmlDtdPtr dtd = NULL;
    int err;

    err = gretl_xml_open_doc_root(fname, NULL, &doc, NULL);
    if (err) {
	pprintf(prn, "Couldn't parse %s\n", fname);
	return 1;
    }

    *dtdname = '\0';

    if (opt & OPT_D) {
	const char *dpath = get_optval_string(MAKEPKG, OPT_D);

	if (dpath != NULL && *dpath != '\0') {
	    strcat(dtdname, dpath);
	}
    } else {
	sprintf(dtdname, "%sfunctions%cgretlfunc.dtd", gretl_home(), SLASH);
    }

    if (*dtdname != '\0') {
	dtd = xmlParseDTD(NULL, (const xmlChar *) dtdname);
    }

    if (dtd == NULL) {
	pputs(prn, "Couldn't open DTD to check package\n");
    } else {
	const char *pkgname = path_last_element(fname);
	xmlValidCtxtPtr cvp = xmlNewValidCtxt();

	if (cvp == NULL) {
	    pputs(prn, "Couldn't get an XML validation context\n");
	    xmlFreeDtd(dtd);
	    xmlFreeDoc(doc);
	    return 0;
	}

	cvp->userData = (void *) prn;
	cvp->error    = (xmlValidityErrorFunc) pprintf;
	cvp->warning  = (xmlValidityWarningFunc) pprintf;

	pprintf(prn, "Checking against %s\n", dtdname);

	if (!xmlValidateDtd(cvp, doc, dtd)) {
	    err = 1;
	} else {
	    pprintf(prn, _("%s: validated against DTD OK"), pkgname);
	    pputc(prn, '\n');
	}

	xmlFreeValidCtxt(cvp);
	xmlFreeDtd(dtd);
    } 

    xmlFreeDoc(doc);

    return err;
}

static int get_gfn_info_for_zip (const char *fname,
				 int *pdfdoc,
				 char ***datafiles,
				 int *n_datafiles)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr sub;
    int found = 0;
    int err = 0;

    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (err) {
	return err;
    }

    node = node->xmlChildrenNode;
    while (node != NULL && found < 2) {
	if (!xmlStrcmp(node->name, (XUC) "gretl-function-package")) {
	    sub = node->xmlChildrenNode;
	    while (sub != NULL && found < 2) {
		if (!xmlStrcmp(sub->name, (XUC) "help")) {
		    char *s = NULL;
		    
		    gretl_xml_node_get_trimmed_string(sub, doc, &s);
		    *pdfdoc = has_suffix(s, ".pdf");
		    free(s);
		    found++;
		} else if (!xmlStrcmp(sub->name, (XUC) "data-files")) {
		    *datafiles =
			gretl_xml_get_strings_array(sub, doc, n_datafiles,
						    0, &err);
		    found++;
		} else if (!xmlStrcmp(sub->name, (XUC) "gretl-function")) {
		    /* we've overshot */
		    found = 2;
		}
		sub = sub->next;
	    }
	}
	node = node->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return err;
}

static int cli_build_zip_package (const char *fname,
				  const char *gfnname,
				  fnpkg *pkg,
				  PRN *prn)
{
    char **datafiles = NULL;
    int n_datafiles = 0;
    int pdfdoc = 0;
    int freeit = 0;
    int err = 0;

    if (pkg != NULL) {
	datafiles = pkg->datafiles;
	n_datafiles = pkg->n_files;
	pdfdoc = has_suffix(pkg->help, ".pdf");
    } else {
	/* grope in the gfn file */
	err = get_gfn_info_for_zip(gfnname,
				   &pdfdoc,
				   &datafiles,
				   &n_datafiles);
	freeit = 1;
    }
	
    if (!err) {
	err = package_make_zipfile(gfnname,
				   pdfdoc,
				   datafiles,
				   n_datafiles,
				   NULL,
				   fname,
				   OPT_NONE,
				   prn);
    }

    if (freeit) {
	strings_array_free(datafiles, n_datafiles);
    }

    return err;
}

/**
 * create_and_write_function_package:
 * @fname: filename for function package.
 * @opt: may include OPT_I to write a package-index entry,
 * OPT_T to write translatable strings.
 * @prn: printer struct for feedback.
 *
 * Create a package based on the functions currently loaded, and
 * write it out as an XML file. Responds to the makepkg command.
 *
 * Returns: 0 on success, non-zero on error.
 */

int create_and_write_function_package (const char *fname, 
				       gretlopt opt,
				       PRN *prn)
{
    char gfnname[FILENAME_MAX];
    fnpkg *pkg = NULL;
    int build_gfn = 1;
    int build_zip = 0;
    int err = 0;

    if (has_suffix(fname, ".zip")) {
	/* building a zip package */
	FILE *fp;

	switch_ext(gfnname, fname, "gfn");
	fp = gretl_fopen(gfnname, "r");
	if (fp != NULL) {
	    /* gfn is already made */
	    build_gfn = 0;
	    fclose(fp);
	}
	build_zip = 1;
    } else {
	/* just building a gfn file */
	strcpy(gfnname, fname);
    }

    if (build_gfn && n_free_functions() == 0) {
	gretl_errmsg_set(_("No functions are available for packaging at present."));
	err = E_DATA;
    } else if (build_gfn) {
	pkg = new_pkg_from_spec_file(gfnname, opt, prn, &err);
	if (pkg != NULL) {
	    err = function_package_write_file(pkg);
	    if (!err) {
		err = cli_validate_package_file(gfnname, opt, prn);
		/* should we delete @gfnname ? */
	    }
	    if (!err) {
		if (opt & OPT_T) {
		    package_write_translatable_strings(pkg, prn);
		}
		if (opt & OPT_I) {
		    package_write_index(pkg, prn);
		}		
	    }
	}
    }

    if (!err && build_zip) {
	err = cli_build_zip_package(fname, gfnname, pkg, prn);
    }

    return err;
}

/**
 * function_package_get_name:
 * @pkg: function package.
 *
 * Returns: the name of the package.
 */

const char *function_package_get_name (fnpkg *pkg)
{
    return (pkg != NULL)? pkg->name : NULL;
}

static int maybe_replace_string_var (char **svar, const char *src)
{
    if (src == NULL) {
	gretl_errmsg_set("string value is missing");
	return E_DATA;
    } else {
	free(*svar);
	*svar = gretl_strdup(src);
	return (*svar == NULL)? E_ALLOC : 0;
    }
}

/* unlike the case above, here we'll accept NULL for @src, to wipe
   out an existing string var */

static int maybe_replace_optional_string_var (char **svar, const char *src)
{
    if (src == NULL) {
	free(*svar);
	*svar = NULL;
	return 0;
    } else {
	free(*svar);
	*svar = gretl_strdup(src);
	return (*svar == NULL)? E_ALLOC : 0;
    }
}

/* Called (indirectly) from GUI function packager. Note that we're
   careful not to touch UFUN_PRIVATE or UFUN_PLUGIN on the uf->flags
   side, since these flags are not represented in @attrs.
*/

static void pkg_set_gui_attrs (fnpkg *pkg, const unsigned char *attrs)
{
    ufunc *uf;
    int i, r;

    for (r=1; r<UFUN_GUI_PRECHECK; r++) {
	uf = NULL;
	for (i=0; i<n_ufuns; i++) {
	    if (ufuns[i]->pkg == pkg && ufuns[i]->pkg_role == r) {
		uf = ufuns[i];
		break;
	    }
	}
	if (uf != NULL) {
	    if (attrs[r-1] & UFUN_NOPRINT) {
		uf->flags |= UFUN_NOPRINT;
	    } else {
		uf->flags &= ~UFUN_NOPRINT;
	    }
	    if (attrs[r-1] & UFUN_MENU_ONLY) {
		uf->flags |= UFUN_MENU_ONLY;
	    } else {
		uf->flags &= ~UFUN_MENU_ONLY;
	    }	    
	}
    }
}

/* Called (indirectly) from GUI function packager. Note that we scrub
   UFUN_PRIVATE and UFUN_PLUGIN from the user function flags passed to
   the packager: UFUN_PLUGIN is irrelevant and UFUN_PRIVATE is handled
   by a different mechanism.
*/

static void pkg_get_gui_attrs (fnpkg *pkg, unsigned char *attrs)
{
    ufunc *uf;
    int i, r;

    for (r=1; r<UFUN_GUI_PRECHECK; r++) {
	uf = NULL;
	for (i=0; i<n_ufuns; i++) {
	    if (ufuns[i]->pkg == pkg && ufuns[i]->pkg_role == r) {
		uf = ufuns[i];
		break;
	    }
	}
	if (uf == NULL) {
	    attrs[r-1] = 0;
	} else {
	    attrs[r-1] = uf->flags & ~(UFUN_PRIVATE | UFUN_PLUGIN);
	}
    }
}

static char **pkg_strvar_pointer (fnpkg *pkg, const char *key,
				  int *optional)
{
    *optional = 0;
    
    if (!strcmp(key, "fname")) {
	return &pkg->fname;
    } else if (!strcmp(key, "author")) {
	return &pkg->author;
    } else if (!strcmp(key, "email")) {
	return &pkg->email;
    } else if (!strcmp(key, "version")) {
	return &pkg->version;
    } else if (!strcmp(key, "date")) {
	return &pkg->date;
    } else if (!strcmp(key, "description")) {
	return &pkg->descrip;
    } else if (!strcmp(key, "help")) {
	return &pkg->help;
    } else if (!strcmp(key, "sample-script")) {
	return &pkg->sample;
    }

    *optional = 1;

    if (!strcmp(key, "tags")) {
	return &pkg->tags; /* FIXME should become non-optional */
    } else if (!strcmp(key, "label")) {
	return &pkg->label;
    } else if (!strcmp(key, "menu-attachment")) {
	return &pkg->mpath;
    } else if (!strcmp(key, "gui-help")) {
	return &pkg->gui_help;
    } else if (!strcmp(key, "help-fname")) {
	return &pkg->help_fname;
    } else if (!strcmp(key, "gui-help-fname")) {
	return &pkg->gui_help_fname;
    } else if (!strcmp(key, "sample-fname")) {
	return &pkg->sample_fname;
    }

    return NULL;
}

/* varargs function for setting the properties of a function
   package: the settings take the form of a NULL-terminated
   set of (key, value) pairs.
*/

int function_package_set_properties (fnpkg *pkg, ...)
{
    va_list ap;
    const char *key;
    char **sptr;
    int optional;
    int i, err = 0;

    va_start(ap, pkg);

    for (i=1; !err; i++) {
	key = va_arg(ap, const char *);
	if (key == NULL) {
	    break;
	}

	sptr = pkg_strvar_pointer(pkg, key, &optional);

	if (sptr != NULL) {
	    const char *sval = va_arg(ap, const char *);

	    if (optional) {
		err = maybe_replace_optional_string_var(sptr, sval);
	    } else {
		err = maybe_replace_string_var(sptr, sval);
	    }
	    
	    if (!err && !strcmp(key, "help")) {
		if (!strncmp(sval, "pdfdoc", 6) || has_suffix(sval, ".pdf")) {
		    pkg->uses_subdir = 1;
		}
	    }
	} else if (!strcmp(key, "gui-attrs")) {
	    const unsigned char *np = va_arg(ap, const unsigned char *);

	    pkg_set_gui_attrs(pkg, np);
	} else {
	    int ival = va_arg(ap, int);

	    if (!strcmp(key, "data-requirement")) {
		pkg->dreq = ival;
	    } else if (!strcmp(key, "model-requirement")) {
		pkg->modelreq = ival;
	    } else if (!strcmp(key, "min-version")) {
		pkg->minver = ival;
	    } else if (!strcmp(key, "lives-in-subdir")) {
		pkg->uses_subdir = (ival != 0);
	    }
	} 
    }

    va_end(ap);

    return err;
}

enum {
    PUBLIST,
    GUILIST,
    PRIVLIST
};

/* From a function package get a list of either its public or its
   private functions, in the form of a simple gretl list.  The
   elements of this list are just the positions of these functions in
   the current recorder array for user-functions.  This is fine if
   the information is to be used right away (as in the GUI function
   call dialog), but it should _not_ be assumed that the identifiers
   will remain valid indefinitely.
*/

static int *function_package_get_list (fnpkg *pkg, int code, int n)
{
    int *list = NULL;
    int j = 0;

    if (n > 0) {
 	list = gretl_list_new(n);
	if (list != NULL) {
	    int i, priv, menu_only;

	    for (i=0; i<n_ufuns; i++) {
		if (ufuns[i]->pkg == pkg) {
		    priv = function_is_private(ufuns[i]);
		    menu_only = function_is_menu_only(ufuns[i]);
		    if (code == PRIVLIST && priv) {
			list[++j] = i;
		    } else if (code == PUBLIST && !priv) {
			list[++j] = i;
		    } else if (code == GUILIST && !priv && !menu_only &&
			       !pkg_aux_role(ufuns[i]->pkg_role)) {
			/* in the GUI list of public functions, don't
			   display post-processing functions
			*/
			list[++j] = i;
		    }
		}
	    }
	}
    } 

    if (list != NULL) {
	if (j == 0) {
	    free(list);
	    list = NULL;
	} else {
	    list[0] = j;
	}
    }

    return list;
}

static char *pkg_get_special_func (fnpkg *pkg, int role)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == pkg && ufuns[i]->pkg_role == role) {
	    return g_strdup(ufuns[i]->name);
	}
    }	    

    return NULL;
}

static int pkg_get_special_func_id (fnpkg *pkg, int role)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == pkg && ufuns[i]->pkg_role == role) {
	    return i;
	}
    }	    

    return -1;
}

static void handle_optional_string (char **ps, const char *src)
{
    if (src == NULL) {
	*ps = NULL;
    } else {
	*ps = g_strdup(src);
    }
}

/* varargs function for retrieving the properties of a function
   package: the arguments after @pkg take the form of a
   NULL-terminated set of (key, pointer) pairs; values are written to
   the locations given by the pointers.
*/

int function_package_get_properties (fnpkg *pkg, ...)
{
    va_list ap;
    int npub = 0;
    int npriv = 0;
    int **plist;
    const char *key;
    void *ptr;
    char **ps;
    int *pi;
    int i, err = 0;

    g_return_val_if_fail(pkg != NULL, E_DATA);

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == pkg) {
	    if (function_is_private(ufuns[i])) {
		npriv++;
	    } else {
		npub++;
	    }
	}
    }

    va_start(ap, pkg);

    for (i=1; ; i++) {
	key = va_arg(ap, const char *);
	if (key == NULL) {
	    break;
	}
	ptr = va_arg(ap, void *);
	if (ptr == NULL) {
	    break;
	}
	if (!strcmp(key, "name")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->name);
	} else if (!strcmp(key, "author")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->author);
	} else if (!strcmp(key, "email")) {
	    ps = (char **) ptr;
	    if (pkg->email != NULL) {
		*ps = g_strdup(pkg->email);
	    } else {
		*ps = g_strdup(""); /* ? */
	    }
	} else if (!strcmp(key, "date")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->date);
	} else if (!strcmp(key, "version")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->version);
	} else if (!strcmp(key, "description")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->descrip);
	} else if (!strcmp(key, "help")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->help);
	} else if (!strcmp(key, "gui-help")) {
	    ps = (char **) ptr;
	    handle_optional_string(ps, pkg->gui_help);
	} else if (!strcmp(key, "sample-script")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->sample);
	} else if (!strcmp(key, "help-fname")) {
	    ps = (char **) ptr;
	    handle_optional_string(ps, pkg->help_fname);
	} else if (!strcmp(key, "gui-help-fname")) {
	    ps = (char **) ptr;
	    handle_optional_string(ps, pkg->gui_help_fname);
	} else if (!strcmp(key, "sample-fname")) {
	    ps = (char **) ptr;
	    handle_optional_string(ps, pkg->sample_fname);
	} else if (!strcmp(key, "tags")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->tags);
	} else if (!strcmp(key, "label")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->label);
	} else if (!strcmp(key, "menu-attachment")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->mpath);
	} else if (!strcmp(key, "data-requirement")) {
	    pi = (int *) ptr;
	    *pi = pkg->dreq;
	} else if (!strcmp(key, "model-requirement")) {
	    pi = (int *) ptr;
	    *pi = pkg->modelreq;
	} else if (!strcmp(key, "min-version")) {
	    pi = (int *) ptr;
	    *pi = pkg->minver;
	} else if (!strcmp(key, "lives-in-subdir")) {
	    pi = (int *) ptr;
	    *pi = pkg->uses_subdir;
	} else if (!strcmp(key, "publist")) {
	    plist = (int **) ptr;
	    *plist = function_package_get_list(pkg, PUBLIST, npub);
	} else if (!strcmp(key, "gui-publist")) {
	    plist = (int **) ptr;
	    *plist = function_package_get_list(pkg, GUILIST, npub);
	} else if (!strcmp(key, "privlist")) {
	    plist = (int **) ptr;
	    *plist = function_package_get_list(pkg, PRIVLIST, npriv);
	} else if (!strcmp(key, "gui-main-id")) {
	    pi = (int *) ptr;
	    *pi = pkg_get_special_func_id(pkg, UFUN_GUI_MAIN);
	} else if (!strcmp(key, BUNDLE_PRINT)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_BUNDLE_PRINT);
	} else if (!strcmp(key, BUNDLE_PLOT)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_BUNDLE_PLOT);
	} else if (!strcmp(key, BUNDLE_TEST)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_BUNDLE_TEST);
	} else if (!strcmp(key, BUNDLE_FCAST)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_BUNDLE_FCAST);
	} else if (!strcmp(key, BUNDLE_EXTRA)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_BUNDLE_EXTRA);
	} else if (!strcmp(key, GUI_MAIN)) {	
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_GUI_MAIN);
	} else if (!strcmp(key, GUI_PRECHECK)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_GUI_PRECHECK);
	} else if (!strcmp(key, "gui-attrs")) {
	    unsigned char *s = (unsigned char *) ptr;
	    
	    pkg_get_gui_attrs(pkg, s);
	}
    }

    va_end(ap);

    return err;
}

/* don't tamper with return value! */

const char *function_package_get_string (fnpkg *pkg,
					 const char *id)
{
    if (pkg == NULL || id == NULL) {
	return NULL;
    } else if (!strcmp(id, "help-fname")) {
	return pkg->help_fname;
    } else if (!strcmp(id, "gui-help-fname")) {
	return pkg->gui_help_fname;
    } else if (!strcmp(id, "sample-fname")) {
	return pkg->sample_fname;
    } else if (!strcmp(id, "sample-script")) {
	return pkg->sample;
    } else if (!strcmp(id, "help")) {
	return pkg->help;
    } else if (!strcmp(id, "gui-help")) {
	return pkg->gui_help;
    } else {
	return NULL;
    }
}

char **function_package_get_data_files (fnpkg *pkg, int *n)
{
    char **S = NULL;

    *n = 0;
    if (pkg->datafiles != NULL) {
	S = strings_array_dup(pkg->datafiles, pkg->n_files);
	if (S != NULL) {
	    *n = pkg->n_files;
	}
    }

    return S;
}

int function_package_set_data_files (fnpkg *pkg, char **S, int n)
{
    int err = 0;
    
    if (pkg->datafiles != NULL) {
	strings_array_free(pkg->datafiles, pkg->n_files);
	pkg->datafiles = NULL;
	pkg->n_files = 0;
    }

    if (n > 0) {
	if (S == NULL) {
	    err = E_DATA;
	} else {
	    pkg->datafiles = strings_array_dup(S, n);
	    if (pkg->datafiles == NULL) {
		err = E_ALLOC;
	    } else {
		pkg->n_files = n;
		pkg->uses_subdir = 1;
	    }
	}
    }

    return err;
}

/* quick check to see if there's a gross problem with a package */

static int validate_function_package (fnpkg *pkg)
{
    if (pkg->pub == NULL || pkg->author == NULL ||
	pkg->version == NULL || pkg->date == NULL ||
	pkg->descrip == NULL) {
	return 0;
    }

    return 1;
}

/* for session file use: dump all currently defined functions,
   packaged or not, into a single XML file */

int write_session_functions_file (const char *fname)
{
    FILE *fp;
    int i;

    if (n_ufuns == 0) {
	return 0;
    }

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	return E_FOPEN;
    }

    gretl_xml_header(fp);  
    fputs("<gretl-functions>\n", fp);

    /* first write any loaded function packages */

    for (i=0; i<n_pkgs; i++) {
	if (validate_function_package(pkgs[i])) {
	    real_write_function_package(pkgs[i], fp);
	}
    }

    /* then any unpackaged functions */

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == NULL) {
	    write_function_xml(ufuns[i], fp);
	}
    }

    fputs("</gretl-functions>\n", fp);

    fclose(fp);

    return 0;
}

/* De-allocate a function package: this can be done in either of two
   modes.  If 'full' is non-zero then we destroy all functions that
   are children of the given package, otherwise we leave the
   function-children alone (just 'detaching' them from the parent
   package).
*/

static void real_function_package_free (fnpkg *pkg, int full)
{
    if (pkg != NULL) {
	int i;

	if (full) {
	    for (i=0; i<pkg->n_pub; i++) {
		ufunc_free(pkg->pub[i]);
	    }
	    for (i=0; i<pkg->n_priv; i++) {
		ufunc_free(pkg->priv[i]);
	    }
	} else {
	    for (i=0; i<n_ufuns; i++) {
		if (ufuns[i]->pkg == pkg) {
		    /* remove package info */
		    ufuns[i]->pkg = NULL;
		    set_function_private(ufuns[i], FALSE);
		}
	    }
	}

	if (pkg->datafiles != NULL && pkg->n_files > 0) {
	    strings_array_free(pkg->datafiles, pkg->n_files);
	}

	free(pkg->pub);
	free(pkg->priv);
	free(pkg->fname);
	free(pkg->author);
	free(pkg->email);
	free(pkg->version);
	free(pkg->date);
	free(pkg->descrip);
	free(pkg->help);
	free(pkg->gui_help);
	free(pkg->sample);
	free(pkg->help_fname);
	free(pkg->gui_help_fname);
	free(pkg->sample_fname);
	free(pkg->tags);
	free(pkg->label);
	free(pkg->mpath);
	free(pkg);
    }
}

static void function_package_free (fnpkg *pkg)
{
    real_function_package_free(pkg, 0);
}

static void function_package_free_full (fnpkg *pkg)
{
    real_function_package_free(pkg, 1);
}

/* is the package with filename @fname already in memory? */

static fnpkg *get_loaded_pkg_by_filename (const char *fname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(fname, pkgs[i]->fname)) {
	    return pkgs[i];
	}
    }

    return NULL;
}

/**
 * function_package_unload_by_filename:
 * @fname: package filename.
 *
 * Unloads the specified function package from memory, if it
 * is currently loaded.  The functions 'owned' by the package
 * are not destroyed; they become available for inclusion in
 * other packages.
 */

void function_package_unload_by_filename (const char *fname)
{
    fnpkg *pkg = get_loaded_pkg_by_filename(fname);

    if (pkg != NULL) {
	real_function_package_unload(pkg, 0);
    }
}

/**
 * function_package_unload_full_by_filename:
 * @fname: package filename.
 *
 * Unloads the specified function package from memory, if it
 * is currently loaded.  The functions 'owned' by the package
 * are also unloaded from memory.
 */

void function_package_unload_full_by_filename (const char *fname)
{
    fnpkg *pkg = get_loaded_pkg_by_filename(fname);

    if (pkg != NULL) {
	real_function_package_unload(pkg, 1);
    }
}

/* append a function package to the recorder array of loaded
   packages */

static int function_package_record (fnpkg *pkg)
{
    fnpkg **tmp;
    int err = 0;

    tmp = realloc(pkgs, (n_pkgs + 1) * sizeof *tmp);

    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	pkgs = tmp;
	pkgs[n_pkgs] = pkg;
	n_pkgs++;
    }

    return err;
}

static int broken_package_error (fnpkg *pkg)
{
    gretl_errmsg_sprintf("'%s': package contains "
			 "duplicated function names",
			 pkg->name);
    return E_DATA;
}

#define fn_redef_msg(s) fprintf(stderr, "Redefining function '%s'\n", s)

/* When loading a private function the only real conflict would be
   with a function of the same name owned by the same package.
   Obviously this shouldn't happen but we'll whack it if it does.
*/

static int load_private_function (ufunc *fun)
{
    const char *targ = fun->name;
    int i, err = 0;

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(targ, ufuns[i]->name)) {
	    if (fun->pkg == ufuns[i]->pkg) {
		err = broken_package_error(fun->pkg);
		break;
	    }
	}
    }

    if (!err) {
	err = add_allocated_ufunc(fun);
    }

    return err;
}

/* When loading a public, packaged function we want to avoid conflicts
   with any non-private functions of the same name.  In case we get a
   conflict with a public member of another package we'd best delete
   the entire conflicting package.
*/

static int load_public_function (ufunc *fun)
{
    const char *targ = fun->name;
    int i, done = 0;
    int err = 0;

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(targ, ufuns[i]->name)) {
	    if (ufuns[i]->pkg == fun->pkg) {
		/* name duplication in package */
		err = broken_package_error(fun->pkg);
	    } else if (ufuns[i]->pkg == NULL) {
		/* conflicting unpackaged function */
		ufunc_free(ufuns[i]);
		ufuns[i] = fun;
		done = 1;
	    } else if (!function_is_private(ufuns[i])) {
		/* got a conflicting package */
		fprintf(stderr, "unloading conflicting package '%s'\n",
			ufuns[i]->pkg->name);
		fprintf(stderr, "(filename '%s')\n", ufuns[i]->pkg->fname); 
		function_package_unload_full(ufuns[i]->pkg);
		break;
	    } 
	}
    }

    if (!err && !done) {
	err = add_allocated_ufunc(fun);
    }

    return err;
}

/* A 'real load' is in contrast to just reading some info from a
   package, as is done in various GUI contexts.  We do the real load
   in response to some GUI commands, the "include" command, and also
   when re-opening a gretl session file that contains function
   packages.
*/

static int real_load_package (fnpkg *pkg)
{
    int i, err = 0;

#if PKG_DEBUG
    fprintf(stderr, "real_load_package:\n loading '%s'\n", pkg->fname);
#endif

    if (pkg->priv != NULL) {
	for (i=0; i<pkg->n_priv && !err; i++) {
	    err = load_private_function(pkg->priv[i]);
	}
    }

    if (!err && pkg->pub != NULL) {
	for (i=0; i<pkg->n_pub && !err; i++) {
	    err = load_public_function(pkg->pub[i]);
	}
    } 

    if (!err) {
	/* add to array of loaded packages */
	err = function_package_record(pkg);
    }

    return err;
}

static const char *data_needs_string (FuncDataReq dr)
{
    if (dr == FN_NEEDS_TS) {
	return N_("Time-series data");
    } else if (dr == FN_NEEDS_QM) {
	return N_("Quarterly or monthly data");
    } else if (dr == FN_NEEDS_PANEL) {
	return N_("Panel data");
    } else if (dr == FN_NODATA_OK) {
	return N_("none");
    } else {
	return N_("some sort of dataset");
    }
}

static void print_package_info (const fnpkg *pkg, const char *fname, PRN *prn)
{
    char vstr[8];
    int remote, pdfdoc;

    remote = (strstr(fname, "dltmp.") != NULL);

    if (!remote && g_path_is_absolute(fname)) {
	pprintf(prn, "<@itl=\"File\">: %s\n\n", fname);
    }

    if (pkg->name[0] == '\0' || pkg->author == NULL ||
	pkg->minver <= 0 || pkg->descrip == NULL ||
	pkg->version == NULL || pkg->date == NULL ||
	pkg->help == NULL) {
	pprintf(prn, "\nBroken package! Basic information is missing\n");
	return;
    }

    gretl_version_string(vstr, pkg->minver);
    pdfdoc = has_suffix(pkg->help, ".pdf");

    pprintf(prn, "<@itl=\"Package\">: %s %s (%s)\n", pkg->name, pkg->version,
	    pkg->date);
    pprintf(prn, "<@itl=\"Author\">: %s\n", pkg->author);
    if (pkg->email != NULL && *pkg->email != '\0') {
	pprintf(prn, "<@itl=\"Email\">: %s\n", pkg->email);
    }
    pprintf(prn, "<@itl=\"Required gretl version\">: %s\n", vstr);
    pprintf(prn, "<@itl=\"Data requirement\">: %s\n", _(data_needs_string(pkg->dreq)));
    pprintf(prn, "<@itl=\"Description\">: %s\n", gretl_strstrip(pkg->descrip));

    if (pdfdoc) {
	const char *s = strrchr(pkg->help, ':');

	if (remote) {
	    pprintf(prn, "<@itl=\"Documentation\">: %s\n\n", s + 1);
	} else {
	    gchar *localpdf = g_strdup(fname);
	    gchar *p = strrchr(localpdf, '.');

	    *p = '\0';
	    strncat(p, ".pdf", 4);
	    pprintf(prn, "<@itl=\"Documentation\">: <@adb=\"%s\">\n\n", localpdf);
	    g_free(localpdf);
	}	
    } else {
	pputc(prn, '\n');
    }

    if (pkg->pub != NULL) {
	if (pkg->n_pub == 1) {
	    if (strcmp(pkg->pub[0]->name, pkg->name)) {
		pputs(prn, "<@itl=\"Public interface\">: ");
		pputs(prn, "<mono>\n");
		pprintf(prn, "%s()\n", pkg->pub[0]->name);
		pputs(prn, "</mono>\n\n");
	    }
	} else {
	    int i;
	    
	    pputs(prn, "<@itl=\"Public interfaces\">:\n\n");
	    pputs(prn, "<mono>\n");
	    for (i=0; i<pkg->n_pub; i++) {
		pprintf(prn, "  %s()\n", pkg->pub[i]->name);
	    }
	    pputs(prn, "</mono>\n\n");
	}
    }    

    if (!pdfdoc) {
	pputs(prn, "<@itl=\"Help text\">:\n\n");
	pputs(prn, "<mono>\n");
	pputs(prn, pkg->help);
	pputs(prn, "\n\n");
	pputs(prn, "</mono>\n");
    }

    if (remote && pkg->sample != NULL) {
	pputs(prn, "<@itl=\"Sample script\">:\n\n");
	pputs(prn, "<code>\n");
	pputs(prn, pkg->sample);
	pputs(prn, "</code>\n");
	pputc(prn, '\n');
    }
}

static void print_package_code (const fnpkg *pkg,
				int tabwidth,
				PRN *prn)
{
    int i;

    if (pkg->priv != NULL) {
	pputs(prn, "# private functions\n\n");
	for (i=0; i<pkg->n_priv; i++) {
	    gretl_function_print_code(pkg->priv[i], tabwidth, prn);
	    pputc(prn, '\n');
	}
    }

    if (pkg->pub != NULL) {
	pputs(prn, "# public functions\n\n");
	for (i=0; i<pkg->n_pub; i++) {
	    gretl_function_print_code(pkg->pub[i], tabwidth, prn);
	    pputc(prn, '\n');
	}
    }
}

/* allocate a fnpkg structure and read from XML file into it */

static fnpkg * 
real_read_package (xmlDocPtr doc, xmlNodePtr node, const char *fname, 
		   int *err)
{
    xmlNodePtr cur;
    fnpkg *pkg;
    char *tmp = NULL;
    int id;

#if PKG_DEBUG
    fprintf(stderr, "real_read_package: fname='%s'\n", fname);
#endif

    pkg = function_package_alloc(fname);
    if (pkg == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    gretl_xml_get_prop_as_string(node, "name", &tmp);

    if (tmp == NULL) {
	fprintf(stderr, "real_read_package: package has no name\n");
	*err = E_DATA;
	function_package_free(pkg);
	return NULL;
    }

    strncat(pkg->name, tmp, FN_NAMELEN - 1);
    free(tmp);

    if (gretl_xml_get_prop_as_bool(node, NEEDS_TS)) {
	pkg->dreq = FN_NEEDS_TS;
    } else if (gretl_xml_get_prop_as_bool(node, NEEDS_QM)) {
	pkg->dreq = FN_NEEDS_QM;
    } else if (gretl_xml_get_prop_as_bool(node, NEEDS_PANEL)) {
	pkg->dreq = FN_NEEDS_PANEL;
    } else if (gretl_xml_get_prop_as_bool(node, NO_DATA_OK)) {
	pkg->dreq = FN_NODATA_OK;
    }

    if (gretl_xml_get_prop_as_string(node, "model-requirement", &tmp)) {
	pkg->modelreq = gretl_command_number(tmp);
	free(tmp);
    }    

    if (gretl_xml_get_prop_as_string(node, "minver", &tmp)) {
	pkg->minver = gretl_version_number(tmp);
	free(tmp);
    }

    if (gretl_xml_get_prop_as_int(node, "ID", &id)) {
	pkg->ID = id;
    }

    pkg->uses_subdir = gretl_xml_get_prop_as_bool(node, "lives-in-subdir");

    cur = node->xmlChildrenNode;
    
    while (cur != NULL) {
	if (!xmlStrcmp(cur->name, (XUC) "author")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->author);
	    gretl_xml_get_prop_as_string(cur, "email", &pkg->email);
	} else if (!xmlStrcmp(cur->name, (XUC) "version")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->version);
	} else if (!xmlStrcmp(cur->name, (XUC) "date")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->date);
	} else if (!xmlStrcmp(cur->name, (XUC) "description")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->descrip);
	} else if (!xmlStrcmp(cur->name, (XUC) "help")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->help);
	    gretl_xml_get_prop_as_string(cur, "filename", &pkg->help_fname);
	} else if (!xmlStrcmp(cur->name, (XUC) "gui-help")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->gui_help);
	    gretl_xml_get_prop_as_string(cur, "filename", &pkg->gui_help_fname);
	} else if (!xmlStrcmp(cur->name, (XUC) "sample-script")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->sample);
	    gretl_xml_get_prop_as_string(cur, "filename", &pkg->sample_fname);
	} else if (!xmlStrcmp(cur->name, (XUC) "tags")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->tags);
	} else if (!xmlStrcmp(cur->name, (XUC) "label")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->label);
	} else if (!xmlStrcmp(cur->name, (XUC) "menu-attachment")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->mpath);
	} else if (!xmlStrcmp(cur->name, (XUC) "data-files")) {
	    pkg->datafiles =
		gretl_xml_get_strings_array(cur, doc, &pkg->n_files,
					    0, err);
	}

	cur = cur->next;
    }

    cur = node->xmlChildrenNode;
    while (cur != NULL && !*err) {
        if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
	    *err = read_ufunc_from_xml(cur, doc, pkg);
	}
	cur = cur->next;
    }

#if PKG_DEBUG
    fprintf(stderr, "real_read_package: err = %d\n", *err);
#endif

    return pkg;
}

/* read aggregated file that may contain one or more function packages
   and one or more "loose" functions */

int read_session_functions_file (const char *fname)
{
    fnpkg *pkg = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;
    int err = 0;

#if PKG_DEBUG
    fprintf(stderr, "read_session_functions_file: starting on '%s'\n", fname);
#endif

    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (err) {
	return err;
    }

    /* first get any function packages from this file */
    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "gretl-function-package")) {
	    pkg = real_read_package(doc, cur, fname, &err);
	    if (!err) {
		err = real_load_package(pkg);
	    }
	} 
	cur = cur->next;
    }

    /* then get any unpackaged functions */
    if (!err) {
	cur = node->xmlChildrenNode;
	while (cur != NULL && !err) {
	    if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
		err = read_ufunc_from_xml(cur, doc, NULL);
	    }
	    cur = cur->next;
	}
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

#if PKG_DEBUG
    fprintf(stderr, "read_session_functions_file: returning %d\n", err);
#endif

    return err;
}

/* parse an XML function package file and return an allocated
   package struct with the functions attached */

static fnpkg *read_package_file (const char *fname, int *err)
{
    fnpkg *pkg = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;

#if PKG_DEBUG
    fprintf(stderr, "read_package_file: got '%s'\n", fname);
#endif

    *err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (*err) {
	return NULL;
    }

    cur = node->xmlChildrenNode;
    while (cur != NULL && !*err) {
	if (!xmlStrcmp(cur->name, (XUC) "gretl-function-package")) {
	    pkg = real_read_package(doc, cur, fname, err);
	    break;
	} 
	cur = cur->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    if (!*err && pkg == NULL) {
	*err = E_DATA;
    }

#if PKG_DEBUG
    fprintf(stderr, "read_function_package: err = %d\n", *err);
#endif

    return pkg;
}

/** 
 * function_package_is_loaded:
 * @fname: full path to gfn file.
 *
 * Returns: 1 if the function package with filename @fname is
 * loaded in memory, otherwise 0.
 */

int function_package_is_loaded (const char *fname)
{
    return (get_loaded_pkg_by_filename(fname) != NULL);
}

/** 
 * function_package_is_loaded:
 * @gfnname: basename of gfn file, including suffix.
 *
 * Returns: 1 if the function package with basename @gfnname is
 * loaded in memory, otherwise 0.
 */

int gfn_is_loaded (const char *gfnname)
{
    int ret = 0;
    
    if (strchr(gfnname, SLASH) == NULL && has_suffix(gfnname, ".gfn")) {
	/* plain basename of gfn file */
	gchar *test = g_strdup(gfnname);
	gchar *p = strrchr(test, '.');

	*p = '\0';
	if (get_function_package_by_name(test)) {
	    ret = 1;
	}
	g_free(test);
    }

    return ret;
}

/** 
 * load_function_package_by_filename:
 * @fname: full path to gfn file.
 *
 * Loads the function package located by @fname into
 * memory, if possible. Supports gretl's "include" command
 * for gfn files.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int load_function_package_by_filename (const char *fname, PRN *prn)
{
    fnpkg *pkg = NULL;
    int err = 0;

    if (function_package_is_loaded(fname)) {
	/* already loaded: no-op */
	fprintf(stderr, "load_function_package_by_filename:\n"
		" '%s' is already loaded\n", fname);
    } else {
	pkg = read_package_file(fname, &err);
	if (!err) {
	    err = real_load_package(pkg);
	} 
    }

    if (err) {
	fprintf(stderr, "load function package: failed on %s\n", fname);
    } else if (pkg != NULL && prn != NULL) {
	pprintf(prn, "%s %s, %s (%s)\n", pkg->name, pkg->version, 
		pkg->date, pkg->author);
    }

    return err;
}

/** 
 * get_function_package_by_filename:
 * @fname: gfn filename.
 * @err: location to receive error code.
 *
 * If the package whose filename is @fname is already loaded,
 * returns the package pointer, otherwise attempts to load the
 * package from file. Similar to load_function_package_from_file()
 * but returns the package pointer rather than just a status
 * code.
 *
 * Returns: package-pointer on success, NULL on error.
 */

fnpkg *get_function_package_by_filename (const char *fname, int *err)
{
    fnpkg *pkg = NULL;
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(fname, pkgs[i]->fname)) {
	    pkg = pkgs[i];
	    break;
	}
    }

    if (pkg == NULL) {
	pkg = read_package_file(fname, err);
	if (!*err) {
	    *err = real_load_package(pkg);
	    if (*err) {
		pkg = NULL;
	    }
	}
    } 	

    return pkg;
}

/** 
 * get_function_package_by_name:
 * @pkgname: name of function package.
 *
 * Returns: pointer to function package if a package named
 * @pkgname is already in memory, otherwise NULL.
 */

fnpkg *get_function_package_by_name (const char *pkgname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(pkgname, pkgs[i]->name)) {
	    return pkgs[i];
	}
    }

    return NULL;
}

/** 
 * get_function_package_path_by_name:
 * @pkgname: name of function package.
 *
 * Returns: the full path to the named function package if it 
 * is already in memory, otherwise NULL.
 */

const char *get_function_package_path_by_name (const char *pkgname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(pkgname, pkgs[i]->name)) {
	    return pkgs[i]->fname;
	}
    }

    return NULL;
}

/* Retrieve summary info, sample script, or code listing for a
   function package, identified by its filename.  This is called
   (indirectly) from the GUI (see below for the actual callbacks).
*/

static int real_print_gfn_data (const char *fname, PRN *prn,
				int tabwidth, int task)
{
    fnpkg *pkg;
    int free_pkg = 0;
    int err = 0;

    pkg = get_loaded_pkg_by_filename(fname);

#if PKG_DEBUG
    fprintf(stderr, "real_print_gfn_data: fname='%s', pkg=%p\n", 
	    fname, (void *) pkg);
#endif

    if (pkg == NULL) {
	/* package is not loaded, read it now */
	pkg = read_package_file(fname, &err);
	free_pkg = 1;
    }

    if (!err) {
	if (task == FUNCS_INFO) {
	    print_package_info(pkg, fname, prn);
	} else if (task == FUNCS_SAMPLE) {
	    pputs(prn, pkg->sample);
	} else {
	    print_package_code(pkg, tabwidth, prn);
	}
	if (free_pkg) {
	    function_package_free_full(pkg);
	}
    }

#if PKG_DEBUG
    fprintf(stderr, "real_get_function_file_info: err = %d (free_pkg = %d)\n", 
	    err, free_pkg);
#endif

    return err;
}

/* callback used in the GUI function package browser */

int print_function_package_info (const char *fname, PRN *prn)
{
    return real_print_gfn_data(fname, prn, 0, FUNCS_INFO);
}

/* callback used in the GUI function package browser */

int print_function_package_code (const char *fname, int tabwidth,
				 PRN *prn)
{
    return real_print_gfn_data(fname, prn, tabwidth, FUNCS_CODE);
}

/* callback used in the GUI function package browser */

int print_function_package_sample (const char *fname, int tabwidth,
				   PRN *prn)
{
    return real_print_gfn_data(fname, prn, tabwidth, FUNCS_SAMPLE);
}

/* Read the header from a function package file -- this is used when
   displaying the available packages in the GUI.  We write the
   description into *pdesc and a string representation of version
   number into *pver.
*/

int get_function_file_header (const char *fname, char **pdesc, 
			      char **pver, int *pdfdoc)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr sub;
    int docdone = 0;
    int err = 0;

    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (err) {
	return err;
    }

    if (pdfdoc == NULL) {
	/* not wanted, so count it as "done" already */
	docdone = 1;
    } else {
	*pdfdoc = 0;
    }

    node = node->xmlChildrenNode;
    while (node != NULL) {
	if (!xmlStrcmp(node->name, (XUC) "gretl-function-package")) {
	    sub = node->xmlChildrenNode;
	    while (sub != NULL) {
		if (!xmlStrcmp(sub->name, (XUC) "description")) {
		    gretl_xml_node_get_trimmed_string(sub, doc, pdesc);
		} else if (!xmlStrcmp(sub->name, (XUC) "version")) {
		    gretl_xml_node_get_trimmed_string(sub, doc, pver);
		} else if (pdfdoc != NULL && !xmlStrcmp(sub->name, (XUC) "help")) {
		    char *tmp = NULL;
		    
		    gretl_xml_node_get_trimmed_string(sub, doc, &tmp);
		    if (tmp != NULL) {
			if (!strncmp(tmp, "pdfdoc", 6) || has_suffix(tmp, ".pdf")) {
			    *pdfdoc = 1;
			}
			free(tmp);
		    }
		    docdone = 1;
		}
		if (*pdesc != NULL && *pver != NULL && docdone) {
		    break;
		}
		sub = sub->next;
	    }
	    if (*pdesc != NULL && *pver != NULL) {
		/* already got what we want */
		break;
	    }
	}
	node = node->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    if (*pdesc == NULL) {
	*pdesc = gretl_strdup(_("No description available"));
    }

    if (*pver == NULL) {
	*pver = gretl_strdup("unknown");
    }

    if (*pdesc == NULL || *pver == NULL) {
	err = E_ALLOC;
    }

    return err;
}

int package_has_menu_attachment (const char *fname,
				 char **pkgname,
				 char **attach,
				 char **label)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr sub;
    char *tmp = NULL;
    int found = 0;
    int ftarg = 1;
    int stop = 0;
    int err = 0;

    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (err) {
	return 0;
    }

    if (label != NULL) {
	ftarg++;
    }

    node = node->xmlChildrenNode;
    
    while (!stop && node != NULL) {
	if (!xmlStrcmp(node->name, (XUC) "gretl-function-package")) {
	    if (pkgname != NULL) {
		gretl_xml_get_prop_as_string(node, "name", pkgname);
	    }
	    sub = node->xmlChildrenNode;
	    while (!stop && sub != NULL) {
		if (!xmlStrcmp(sub->name, (XUC) "menu-attachment")) {
		    gretl_xml_node_get_trimmed_string(sub, doc, &tmp);
		    if (tmp != NULL) {
			if (++found == ftarg) {
			    stop = 1;
			}
			if (attach != NULL) {
			    *attach = tmp;
			} else {
			    free(tmp);
			}
		    }
		} else if (label != NULL &&
			   !xmlStrcmp(sub->name, (XUC) "label")) {
		    gretl_xml_node_get_trimmed_string(sub, doc, &tmp);
		    *label = tmp;
		    if (++found == ftarg) {
			stop = 1;
		    }
		} else if (!xmlStrcmp(sub->name, (XUC) "help")) {
		    /* we've overshot */
		    stop = 1;
		}
		sub = sub->next;
	    }
	}
	node = node->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return found;
}

int package_needs_zipping (const char *fname,
			   int *pdfdoc,
			   char ***datafiles,
			   int *n_files)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr sub;
    char *tmp = NULL;
    int stop = 0;
    int retmax = 1;
    int ret = 0;
    int err = 0;

    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (err) {
	return 0;
    }

    if (datafiles != NULL) {
	retmax = 2;
    }

    node = node->xmlChildrenNode;
    
    while (!stop && node != NULL) {
	if (!xmlStrcmp(node->name, (XUC) "gretl-function-package")) {
	    sub = node->xmlChildrenNode;
	    while (!stop && sub != NULL) {
		if (!xmlStrcmp(sub->name, (XUC) "help")) {
		    gretl_xml_node_get_trimmed_string(sub, doc, &tmp);
		    if (tmp != NULL && !strncmp(tmp, "pdfdoc:", 7)) {
			if (pdfdoc != NULL) {
			    *pdfdoc = 1;
			}
			ret++;
		    }
		    free(tmp);
		} else if (!xmlStrcmp(sub->name, (XUC) "data-files")) {
		    if (datafiles != NULL) {
			*datafiles = gretl_xml_get_strings_array(sub, doc, n_files,
								 0, &err);
		    }
		    ret++;
		} else if (!xmlStrcmp(sub->name, (XUC) "gretl-function")) {
		    /* we've overshot */
		    stop = 1;
		}
		if (ret == retmax) {
		    stop = 1;
		} else {
		    sub = sub->next;
		}
	    }
	}
	node = node->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return ret;
}

double function_package_get_version (const char *fname)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node;
    double version = NADBL;
    int err;

    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
 
    if (!err) {
	xmlNodePtr sub;
	int found = 0;
	
	node = node->xmlChildrenNode;
	while (node != NULL && !found) {
	    if (!xmlStrcmp(node->name, (XUC) "gretl-function-package")) {
		sub = node->xmlChildrenNode;
		while (sub != NULL && !found) {
		    if (!xmlStrcmp(sub->name, (XUC) "version")) {
			gretl_xml_node_get_double(sub, doc, &version);
			found = 1;
		    }
		    sub = sub->next;
		}
	    }
	    node = node->next;
	}
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return version;
}

/**
 * gretl_is_public_user_function:
 * @name: name to test.
 *
 * Returns: 1 if @name is the name of a user-defined
 * function that is not a private member of a function
 * package, otherwise 0.
 */

int gretl_is_public_user_function (const char *name)
{
    ufunc *fun = get_user_function_by_name(name);

    return (fun != NULL && !function_is_private(fun));
}

void set_current_function_package (fnpkg *pkg)
{
    current_pkg = pkg;
}

static int skip_private_func (ufunc *ufun)
{
    int skip = 0;

    if (function_is_private(ufun)) {
	/* skip it, unless we're "authorized" to edit it,
	   i.e. coming from the gui package editor */
	skip = 1;
	if (ufun->pkg != NULL && ufun->pkg == current_pkg) {
	    skip = 0;
	}
    }

    return skip;
}

static int check_func_name (const char *fname, ufunc **pfun, PRN *prn)
{
    int i, err = 0;

#if FN_DEBUG
    fprintf(stderr, "check_func_name: '%s'\n", fname);
#endif

    if (!isalpha((unsigned char) *fname)) {
	gretl_errmsg_set(_("Function names must start with a letter"));
	err = 1;
    } else if (gretl_command_number(fname)) {
	gretl_errmsg_sprintf(_("'%s' is the name of a gretl command"),
			     fname);
	err = 1;
    } else if (function_lookup(fname)) {
	gretl_errmsg_sprintf(_("'%s' is the name of a built-in function"),
			     fname);
	err = 1;
    } else {
	for (i=0; i<n_ufuns; i++) {
	    if (!skip_private_func(ufuns[i]) && !strcmp(fname, ufuns[i]->name)) {
#if FN_DEBUG
		fprintf(stderr, "'%s': found an existing function of this name\n", fname);
#endif
		if (pfun != NULL) {
		    clear_ufunc_data(ufuns[i]);
		    *pfun = ufuns[i];
		} else {
		    ufunc_unload(ufuns[i]);
		}
		break;
	    }
	}
    }
    
    return err;
}

static int maybe_delete_function (const char *fname, PRN *prn)
{
    ufunc *fun = get_user_function_by_name(fname);
    int err = 0;

    if (fun == NULL) {
	; /* no-op */
    } else if (function_in_use(fun)) {
	gretl_errmsg_sprintf("%s: function is in use", fname);
	err = 1;
    } else if (fun->pkg != NULL) {
	gretl_errmsg_sprintf("%s: function belongs to package", fname);
	err = 1;
    } else {
	ufunc_unload(fun);
	if (gretl_messages_on()) {
	    pprintf(prn, _("Deleted function '%s'\n"), fname);
	}
    } 

    return err;
}

/* next: apparatus for parsing function definitions */

static int comma_count (const char *s)
{
    int quoted = 0;
    int braced = 0;
    int nc = 0;

    while (*s) {
	if (!quoted) {
	    if (*s == '{') {
		braced++;
	    } else if (*s == '}') {
		braced--;
	    }
	}
	if (*s == '"') {
	    quoted = !quoted;
	} else if (!quoted && !braced && *s == ',') {
	    nc++;
	}
	s++;
    }

    return nc;
}

static int check_parm_min_max (fn_param *p, const char *name, 
			       int *nvals)
{
    int err = 0;

    if (p->type != GRETL_TYPE_DOUBLE && na(p->deflt)) {
	gretl_errmsg_sprintf("%s: parameters of type %s cannot have a "
			     "default value of NA", name,
			     gretl_type_get_name(p->type));
	err = E_DATA;
    }

    if (!err && !na(p->min) && !na(p->max)) {
	if (p->min > p->max) {
	    gretl_errmsg_sprintf("%s: min > max", name);
	    err = E_DATA;
	} else if (!na(p->step) && p->step > p->max - p->min) {
	    gretl_errmsg_sprintf("%s: step > max - min", name);
	    err = E_DATA;
	}
	if (!err) {
	    *nvals = (int) p->max - (int) p->min + 1;
	}
    }

    if (!err && !na(p->deflt) && !default_unset(p)) {
	if (!na(p->min) && p->deflt < p->min) {
	    gretl_errmsg_sprintf("%s: default value out of bounds", name);
	    err = E_DATA;
	} else if (!na(p->max) && p->deflt > p->max) {
	    gretl_errmsg_sprintf("%s: default value out of bounds", name);
	    err = E_DATA;
	}
    } 

    return err;
}

static int colon_count (const char *p)
{
    int n = 0;

    while (*p && *p != ']') {
	if (*p == ':') {
	    n++;
	}
	p++;
    }

    return n;
}

#define VALDEBUG 0

static int read_min_max_deflt (char **ps, fn_param *param, 
			       const char *name, int *nvals)
{
    char *p = *ps;
    int err = 0;

    gretl_push_c_numeric_locale();

    errno = 0;

    if (param->type == GRETL_TYPE_BOOL) {
	if (sscanf(p, "[%lf]", &param->deflt) != 1) {
	    err = E_PARSE;
	}
    } else if (!strncmp(p, "[$xlist]", 8)) {
	param->deflt = INT_USE_XLIST;
    } else {
	int i, len, nf = colon_count(p) + 1;
	char *test, valstr[32];
	char *q = p + 1;
	double x[4] = {NADBL, NADBL, UNSET_VALUE, NADBL}; 

	for (i=0; i<nf && !err; i++) {
	    len = strcspn(q, ":]");
	    if (len > 31) {
		err = E_PARSE;
	    } else if (len > 0) {
		*valstr = '\0';
		strncat(valstr, q, len);
#if VALDEBUG > 1
		fprintf(stderr, "valstr(%d) = '%s'\n", i, valstr);
#endif
		if (!strcmp(valstr, "NA")) {
		    x[i] = NADBL;
		} else {
		    x[i] = strtod(valstr, &test);
		    if (*test != '\0' || errno) {
			err = E_PARSE;
		    }
		}
	    } 
	    q += len + 1;
	}

	if (!err) {
	    if (nf == 1) {
		/* we take a single value with no colons as
		   indicating a default */
		param->deflt = x[0];
	    } else {
		param->min = x[0];
		param->max = x[1];
		param->deflt = x[2];
		param->step = x[3];
	    }

#if VALDEBUG
	    fprintf(stderr, "min %g, max %g, deflt %g, step %g\n", 
		    param->min, param->max, param->deflt, param->step);
#endif

	    err = check_parm_min_max(param, name, nvals);
	}
    }

    gretl_pop_c_numeric_locale();

    if (!err) {
	p = strchr(p, ']');
	if (p == NULL) {
	    err = E_PARSE;
	} else {
	    *ps = p + 1;
	}
    }

    return err;
}

static int read_param_option (char **ps, fn_param *param)
{
    int err = E_PARSE;

#if FNPARSE_DEBUG
    fprintf(stderr, "read_param_option: got '%s'\n", *ps);
#endif

    if (!strncmp(*ps, "[null]", 6)) {
	param->flags |= ARG_OPTIONAL;
	err = 0;
	*ps += 6;
    }

    return err;
}

/* get the descriptive string for a function parameter */

static int read_param_descrip (char **ps, fn_param *param)
{
    char *p = *ps + 1; /* skip opening quote */
    int len = 0;
    int err = E_PARSE;

    while (*p) {
	if (*p == '"') {
	    /* OK, found closing quote */
	    err = 0;
	    break;
	}
	len++;
	p++;
    }

    if (!err && len > 0) {
	p = *ps + 1;
	param->descrip = gretl_strndup(p, len);
	if (param->descrip == NULL) {
	    err = E_ALLOC;
	} else {
	    *ps = p + len + 1;
	}
    }

    return err;
}

/* get the value labels for a function parameter */

static int read_param_labels (char **ps, fn_param *param, 
			      const char *name, int nvals)
{
    char *p = *ps + 1; /* skip opening '{' */
    int len = 0;
    int err = E_PARSE;

    while (*p) {
	if (*p == '}') {
	    /* OK, found closing brace */
	    err = 0;
	    break;
	}
	len++;
	p++;
    }

    if (!err && len > 0) {
	char *tmp;

	p = *ps + 1;
	tmp = gretl_strndup(p, len);

	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    *ps = p + len + 1;
	    param->labels = gretl_string_split_quoted(tmp, &param->nlabels,
						      " ,", &err);
	    free(tmp);
	    if (!err && param->nlabels != nvals) {
		gretl_errmsg_sprintf("%s: found %d values but %d value-labels",
				     name, nvals, param->nlabels);
		err = E_DATA;
	    }
	}
    }

    return err;
}

static void trash_param_info (char *name, fn_param *param)
{
    free(name);
    free(param->descrip);
    param->descrip = NULL;
    if (param->nlabels > 0) {
	strings_array_free(param->labels, param->nlabels);
	param->labels = NULL;
	param->nlabels = 0;
    }
}

static int parse_function_param (char *s, fn_param *param, int i)
{
    GretlType type;
    char tstr[22] = {0};
    char *name;
    int len, nvals = 0;
    int err = 0;

#if FNPARSE_DEBUG
    fprintf(stderr, "parse_function_param: s = '%s'\n", s);
#endif

    while (isspace(*s)) s++;

    /* pick up the "const" flag if present */

    if (!strncmp(s, "const ", 6)) {
	param->flags |= ARG_CONST;
	s += 6;
	while (isspace(*s)) s++;
    }

    /* get parameter type, required */

    len = gretl_namechar_spn(s);
    if (len > 21) {
	err = E_PARSE;
    } else {
	strncat(tstr, s, len);
	s += len;
	while (isspace(*s)) s++;
	if (*s == '*') {
	    strcat(tstr, " *");
	    s++;
	}
	if (*tstr == '\0') {
	    gretl_errmsg_set(_("Expected a type identifier"));
	    err = E_PARSE;
	} else {
	    type = gretl_type_from_string(tstr);
	    if (type == 0) {
		gretl_errmsg_sprintf(_("Unrecognized data type '%s'"), tstr);
		err = E_PARSE;
	    } else if (!ok_function_arg_type(type)) {
		gretl_errmsg_sprintf(_("%s: invalid parameter type"), tstr);
		err = E_INVARG;
	    }		
	} 
    }

    if (err) {
	return err;
    }

    /* now get the required parameter name */
    
    while (isspace(*s)) s++;

    len = gretl_namechar_spn(s);

    if (len == 0) {
	gretl_errmsg_sprintf(_("parameter %d: name is missing"), i + 1);
	err = E_PARSE;
    } else {
	name = gretl_strndup(s, len);
	if (name == NULL) {
	    err = E_ALLOC;
	} else if (gretl_reserved_word(name)) {
	    free(name);
	    err = E_DATA;
	}
    }	

    if (err) {
	return err;
    }

    param->type = type;

    s += len;
    s += strspn(s, " ");

    /* now scan for various optional extras: first we may have
       a [min:max:default] specification, or a [null] spec
       to allow no argument for "pointer"-type arguments
    */

    if (*s == '[') { 
	if (gretl_scalar_type(type)) {
	    err = read_min_max_deflt(&s, param, name, &nvals);
	} else if (arg_may_be_optional(type)) {
	    err = read_param_option(&s, param);
	} else {
	    err = E_PARSE;
	}
	s += strspn(s, " ");
    } 

    /* then we may have a double-quoted descriptive string
       for the parameter */

    if (!err && *s == '"') {
	err = read_param_descrip(&s, param);
	s += strspn(s, " ");
    } 

    /* and finally we may have a set of value-labels enclosed
       in braces, but this is valid only if we have been able
       to determine a definite number of admissible values
       for the parameter
    */

    if (!err && *s == '{') {
	if (nvals == 0) {
	    err = E_PARSE;
	} else {
	    err = read_param_labels(&s, param, name, nvals);
	}
	s += strspn(s, " ");
    }    

    if (!err && *s != '\0') {
	/* got trailing unparseable stuff */
	err = E_PARSE;
    }

    if (!err) {
	param->name = name;
    } else {
	trash_param_info(name, param);
    }

#if FNPARSE_DEBUG
    if (!err) {
	fprintf(stderr, " param[%d] = '%s', ptype = %d\n", 
		i, name, type);
	fprintf(stderr, "  min=%g, max=%g, deflt=%g\n", 
		param->min, param->max, param->deflt);
	if (param->descrip != NULL) {
	    fprintf(stderr, "  comment = '%s'\n", param->descrip);
	}
	if (param->nlabels > 0) {
	    fprintf(stderr, "  value labels: %d\n", param->nlabels);
	}	
    }
#endif

    return err;
}

static void arg_tail_strip (char *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>=0; i--) {
	if (isspace(s[i]) || s[i] == ')') {
	    s[i] = '\0';
	} else {
	    break;
	}
    }
}

static int parse_function_parameters (const char *str,
				      fn_param **pparams,
				      int *pnp,
				      const char *fname,
				      PRN *prn)
{
    fn_param *params = NULL;
    char *p, *s;
    int i, j, np = 0;
    int err = 0;

    s = gretl_strdup(str);
    if (s == NULL) {
	return E_ALLOC;
    }    

    if (!strcmp(s, ")")) {
	/* we got a void function "foo()" */
	free(s);
	return 0;
    }

    /* strip trailing ')' and space */
    arg_tail_strip(s);
    np = comma_count(s) + 1;
    if (np == 1 && !strcmp(s, "void")) {
	free(s);
	return 0;
    }

#if FNPARSE_DEBUG
    fprintf(stderr, "function %s: looking for %d parameters\n",
	    fname, np);
#endif

    params = allocate_params(np);
    if (params == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	int quoted = 0;
	int braced = 0;

	p = s;
	while (*p) {
	    if (!quoted) {
		if (*p == '{') {
		    braced++;
		} else if (*p == '}') {
		    braced--;
		}
	    }
	    if (*p == '"') {
		quoted = !quoted;
	    } else if (!quoted && !braced && *p == ',') {
		*p = '\0';
	    }
	    p++;
	}
	p = s;
	if (braced != 0) {
	    err = E_PARSE;
	}
	for (i=0; i<np && !err; i++) {
	    err = parse_function_param(p, &params[i], i);
	    p += strlen(p) + 1;
	}
    }

    free(s);

    for (i=0; i<np && !err; i++) {
	for (j=i+1; j<np && !err; j++) {
	    if (!strcmp(params[i].name, params[j].name)) {
		gretl_errmsg_sprintf(_("%s: duplicated parameter name '%s'"),
				     fname, params[i].name);
		err = E_DATA;
	    }
	}
    }

    if (err) {
	free_params_array(params, np);
    } else {
	*pparams = params;
	*pnp = np;
    }
    
    return err;
}

static int get_two_words (const char *s, char *w1, char *w2,
			  int *err)
{
    int nf = 0;

    if (strncmp(s, "function ", 9)) {
	*err = E_PARSE;
    } else {
	int n;

	*w1 = *w2 = '\0';

	s += 9;
	s += strspn(s, " ");
	n = strcspn(s, " (");

	if (n == 0 || n >= FN_NAMELEN) {
	    *err = E_PARSE;
	} else {
	    strncat(w1, s, n);
	    nf++;
	    s += n;
	    s += strspn(s, " ");
	    n = strcspn(s, " (");
	}

	if (n > 0 && n < FN_NAMELEN) {
	    strncat(w2, s, n);
	    nf++;
	}
    }

    return nf;
}

/**
 * gretl_start_compiling_function:
 * @line: command line.
 * @prn: printing struct for feedback.
 *
 * Responds to a command of the form "function ...".  In most
 * cases, embarks on compilation of a function, but this
 * also handles the construction "function foo delete".
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_start_compiling_function (const char *line, PRN *prn)
{
    ufunc *fun = NULL;
    fn_param *params = NULL;
    int nf, n_params = 0;
    int rettype = 0;
    char s1[FN_NAMELEN];
    char s2[FN_NAMELEN];
    char *p = NULL, *fname = NULL;
    int err = 0;

    if (gretl_function_depth() > 0) {
	return E_FNEST;
    }

    nf = get_two_words(line, s1, s2, &err);

    if (err) {
	return err;
    } else if (nf < 2) {
	gretl_errmsg_set("A function definition must have a return type and name");
	return E_PARSE;
    }

    if (!strcmp(s2, "clear") || !strcmp(s2, "delete")) {
	return maybe_delete_function(s1, prn);
    }

    /* If we didn't get a special such as "function foo delete",
       then @s1 should contain the return type and @s2 the
       name.
    */

    fname = s2;
    rettype = return_type_from_string(s1, &err);

    if (!err) {
	err = check_func_name(fname, &fun, prn);
    }

    if (!err) {
	/* now for the args bit */
	p = strchr(line, '(');
	if (p == NULL || strchr(p, ')') == NULL) {
	    err = E_PARSE;
	} else {
	    p++; /* skip the left paren */
	}
    }

    if (!err) {
	err = parse_function_parameters(p, &params, &n_params, 
					fname, prn);
    }

    if (err) {
	pprintf(prn, "> %s\n", line);
    }

    if (!err && fun == NULL) {
	fun = add_ufunc(fname);
	if (fun == NULL) {
	    free_params_array(params, n_params);
	    err = E_ALLOC;
	}
    }

#if GLOBAL_TRACE
    fprintf(stderr, "started compiling function %s (err = %d)\n",
	    fname, err);
#endif

    if (!err) {
	strcpy(fun->name, fname);
	fun->params = params;
	fun->n_params = n_params;
	fun->rettype = rettype;
	current_fdef = fun;
	set_compiling_on();
    } else {
	current_fdef = NULL;
    }
    
    return err;
}

#define NEEDS_IF(c) (c == ELSE || c == ELIF || c == ENDIF)

static void python_check (const char *line)
{
    char s1[8], s2[16];

    if (sscanf(line, "%7s %15s", s1, s2) == 2) {
	if (!strcmp(s1, "foreign") && strstr(s2, "ytho")) {
	    compiling_python = 1;
	} else if (!strcmp(s1, "end") && !strcmp(s2, "foreign")) {
	    compiling_python = 0;
	}
    }
}

/**
 * gretl_function_append_line:
 * @line: line of code to append.
 *
 * Continuation of definition of user-function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_function_append_line (const char *line)
{
    static CMD *cmd;
    static int ifdepth;
    char tmpline[MAXLINE];
    ufunc *fun = current_fdef;
    int blank = 0;
    int abort = 0;
    int ignore = 0;
    int err = 0;

    if (fun == NULL) {
	fprintf(stderr, "gretl_function_append_line: fun is NULL\n");
	return 1;
    }

#if FNPARSE_DEBUG
    fprintf(stderr, "gretl_function_append_line: '%s' (idx = %d)\n", 
	    line, fun->line_idx);
#endif

    if (cmd == NULL) {
	cmd = gretl_cmd_new();
	if (cmd == NULL) {
	    return E_ALLOC;
	}
    }

    blank = string_is_blank(line);

    /* below: append line to function, carrying out some
       basic structural checks as we go 
    */

    if (!blank) {
	/* note: avoid losing comment lines */
	strcpy(tmpline, line);
	err = get_command_index(tmpline, FUNC, cmd);
	if (!err) {
	    if (cmd->ci == QUIT) {
		abort = 1;
	    } else if (cmd->flags & CMD_ENDFUN) {
		if (fun->n_lines == 0) {
		    gretl_errmsg_sprintf("%s: empty function", fun->name);
		    err = 1;
		}
		set_compiling_off();
	    } else if (cmd->ci == FUNC) {
		err = E_FNEST;
	    } else if (cmd->ci == IF) {
		ifdepth++;
	    } else if (NEEDS_IF(cmd->ci) && ifdepth == 0) {
		gretl_errmsg_sprintf("%s: unbalanced if/else/endif", fun->name);
		err = E_PARSE;
	    } else if (cmd->ci == ENDIF) {
		ifdepth--;
	    } else if (cmd->ci < 0) {
		ignore = 1;
	    }
	}
    }

    if (abort || err) {
	set_compiling_off();
    }    

    if (compiling) {
	/* actually add the line */
	int i = fun->n_lines;
	
	err = push_function_line(fun, line);
	if (!err) {
	    if (ignore) {
		fun->lines[i].ignore = 1;
	    } else if (!blank) {
		python_check(line);
	    }
	}
    } else {
	/* finished compilation */
	if (!abort && !err && ifdepth != 0) {
	    gretl_errmsg_sprintf("%s: unbalanced if/else/endif", fun->name);
	    err = E_PARSE;
	}
#if GLOBAL_TRACE	
	fprintf(stderr, "finished compiling function %s\n", fun->name);
#endif	
	/* reset static vars */
	gretl_cmd_destroy(cmd);
	cmd = NULL;
	ifdepth = 0;
    }

    if (abort || err) {
	ufunc_unload(fun);
    }	

    return err;
}

/* next block: handling function arguments */

/* Given a list of variables supplied as an argument to a function,
   copy the list under the name assigned by the function and
   make the variables referenced in that list accessible within
   the function. We also handle the case where the given arg
   was not a real list, but we can make a list out of it for the
   purposes of the function.
*/

static int localize_list (fncall *call, fn_arg *arg,
			  fn_param *fp, DATASET *dset)
{
    int *list = NULL;
    int err = 0;

    if (arg == NULL || arg->type == GRETL_TYPE_NONE) {
	/* empty arg -> gives an empty list */
	int tmp[] = {0};

	list = copy_list_as_arg(fp->name, tmp, &err);
    } else if (arg->type == GRETL_TYPE_LIST) {
	/* actual list arg -> copy to function level */
	list = arg->val.list;
	err = copy_as_arg(fp->name, GRETL_TYPE_LIST, list);
    } else if (arg->type == GRETL_TYPE_USERIES) {
	/* series arg -> becomes a singleton list */
	int tmp[] = {1, arg->val.idnum};

	list = copy_list_as_arg(fp->name, tmp, &err);
    } else {
	/* "can't happen" */
	err = E_DATA;
    }

    if (!err && list == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	int i, vi, level = fn_executing + 1;

	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi > 0) {
		if (!in_gretl_list(call->listvars, vi)) {
		    gretl_list_append_term(&call->listvars, vi);
		}
		series_set_stack_level(dset, vi, level);
	    }
	}
    }

#if UDEBUG
    fprintf(stderr, "localize_list (%s): returning %d\n", 
	    fp->name, err);
#endif

    return err;
}

static void maybe_set_arg_const (fn_arg *arg, fn_param *fp)
{
    if (fp->flags & ARG_CONST) {
	/* param is marked CONST directly */
	arg->name = fp->name;
	arg->flags |= ARG_CONST;
    } else if (object_is_const(arg->upname)) {
	/* param is CONST by inheritance */
	arg->name = fp->name;
	arg->flags |= ARG_CONST;
    }
}

/* handle the case where the GUI passed an anonymous
   bundle as a "bundle *" argument */

static int localize_bundle_as_shell (fn_arg *arg, 
				     fn_param *fp)
{
    int err = 0;

    err = arg_add_as_shell(fp->name, GRETL_TYPE_BUNDLE, arg->val.b);

    if (!err) {
	arg->name = fp->name;
    }

    return err;
}

static int localize_const_matrix (fn_arg *arg, fn_param *fp)
{
    user_var *u = get_user_var_by_data(arg->val.m);
    int err = 0;

    if (u == NULL) {
	/* the const argument is an anonymous matrix */
	err = arg_add_as_shell(fp->name, GRETL_TYPE_MATRIX, arg->val.m);
    } else {
	/* a named matrix: in view of its "const-ness" we
	   don't need to copy the data
	*/
	user_var_adjust_level(u, 1);
	user_var_set_name(u, fp->name);
    }

    if (!err) {
	arg->name = fp->name;
	arg->flags |= ARG_CONST;
    }

    return err;
}

static int localize_series_ref (fncall *call, fn_arg *arg, 
				fn_param *fp, DATASET *dset)
{
    int v = arg->val.idnum;

    series_increment_stack_level(dset, v);
    strcpy(dset->varname[v], fp->name);

    if (!in_gretl_list(call->ptrvars, v)) {
	gretl_list_append_term(&call->ptrvars, v);
    }

    maybe_set_arg_const(arg, fp);

    return 0;
}

static int argval_get_int (fn_param *param, double x, int *err)
{
    int ret = gretl_int_from_double(x, err);
    
    if (*err) {
	gretl_errmsg_sprintf(_("%s: expected an integer but found %g"),
			     param->name, x);
    }

    return ret;
}

static int real_add_scalar_arg (fn_param *param, double x)
{
    int err = 0;

    if (na(x)) {
	/* always allow NA for scalar args (?) */
	err = copy_as_arg(param->name, GRETL_TYPE_DOUBLE, &x);
    } else {
	if (param->type == GRETL_TYPE_BOOL) {
	    if (x != 0.0) {
		x = 1.0;
	    }
	} else if (param->type == GRETL_TYPE_INT ||
		   param->type == GRETL_TYPE_OBS) {
	    x = argval_get_int(param, x, &err);
	}

	if (!err) {
	    if ((!na(param->min) && x < param->min) ||
		(!na(param->max) && x > param->max)) {
		gretl_errmsg_sprintf(_("%s: argument value %g is out of bounds"),
				     param->name, x);
		err = E_DATA;
	    } else {
		err = copy_as_arg(param->name, GRETL_TYPE_DOUBLE, &x);
	    }
	}
    }

    return err;
}

/* Scalar function arguments only: if the arg is not supplied, use the
   default that is found in the function specification, if any.
*/

static int add_scalar_arg_default (fn_param *param)
{
    if (default_unset(param)) {
	/* should be impossible here, but... */
	return E_DATA;
    } else {
	return real_add_scalar_arg(param, param->deflt);
    }
}

/* Handle the case of a scalar parameter for which a 1 x 1 matrix 
   was given as argument.
*/

static int do_matrix_scalar_cast (fn_arg *arg, fn_param *param)
{
    gretl_matrix *m = arg->val.m;

    return real_add_scalar_arg(param, m->val[0]);
}

static void fncall_finalize_listvars (fncall *call)
{
    int i, v;

    for (i=call->listvars[0]; i>0; i--) {
	v = call->listvars[i];
	if (in_gretl_list(call->ptrvars, v)) {
	    gretl_list_delete_at_pos(call->listvars, i);
	}
    }

    if (call->listvars[0] == 0) {
	free(call->listvars);
	call->listvars = NULL;
    }
}

static int duplicated_pointer_arg_check (fn_arg *args, int argc)
{
    fn_arg *ai, *aj;
    int i, j, err = 0;

    /* note: a caller cannot be allowed to supply a given 
       variable in pointer form for more than one argument 
       slot in a function call (although ordinary arguments
       may be repeated)
    */

    for (i=0; i<argc && !err; i++) {
	ai = &args[i];
	if (gretl_ref_type(ai->type)) {
	    for (j=i+1; j<argc && !err; j++) {
		aj = &args[j];
		if (gretl_ref_type(aj->type) && 
		    !strcmp(ai->upname, aj->upname)) {
		    gretl_errmsg_set(_("Duplicated pointer argument: not allowed"));
		    err = E_DATA;
		}
	    }
	}
    }

    return err;
}

static int allocate_function_args (fncall *call, DATASET *dset)
{
    ufunc *fun = call->fun;
    fn_arg *arg;
    fn_param *fp;
    int i, err;

    err = duplicated_pointer_arg_check(call->args, call->argc);

    for (i=0; i<call->argc && !err; i++) {
	arg = &call->args[i];
	fp = &fun->params[i];

#if UDEBUG
	fprintf(stderr, "arg[%d], param type %s, arg type %s\n",
		i, gretl_type_get_name(fp->type),
		gretl_type_get_name(arg->type));
#endif

	if (gretl_scalar_type(fp->type)) {
	    if (arg->type == GRETL_TYPE_NONE) {
		err = add_scalar_arg_default(fp);
	    } else if (arg->type == GRETL_TYPE_MATRIX) {
		err = do_matrix_scalar_cast(arg, fp);
	    } else {
		err = real_add_scalar_arg(fp, arg->val.x);
	    }
	} else if (fp->type == GRETL_TYPE_SERIES) {
	    if (arg->type == GRETL_TYPE_USERIES) {
		err = dataset_copy_series_as(dset, arg->val.idnum, fp->name);
	    } else {
		err = dataset_add_series_as(dset, arg->val.px, fp->name);
	    }	    
	} else if (fp->type == GRETL_TYPE_MATRIX) {
	    if (arg->type != GRETL_TYPE_NONE) {
		if (fp->flags & ARG_CONST) {
		    err = localize_const_matrix(arg, fp);
		} else {
		    err = copy_as_arg(fp->name, fp->type, arg->val.m);
		}
	    }
	} else if (fp->type == GRETL_TYPE_BUNDLE) {
	    err = copy_as_arg(fp->name, fp->type, arg->val.b);
	} else if (gretl_array_type(fp->type)) {
	    err = copy_as_arg(fp->name, fp->type, arg->val.a);
	} else if (fp->type == GRETL_TYPE_LIST) {
	    err = localize_list(call, arg, fp, dset);
	} else if (fp->type == GRETL_TYPE_STRING) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = copy_as_arg(fp->name, fp->type, arg->val.str); 
	    }
	} else if (fp->type == GRETL_TYPE_SERIES_REF) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = localize_series_ref(call, arg, fp, dset);
	    }
	} else if (gretl_ref_type(fp->type)) {
	    if (arg->type != GRETL_TYPE_NONE) {
		if (fp->type == GRETL_TYPE_BUNDLE_REF && arg->upname[0] == '\0') {
		    err = localize_bundle_as_shell(arg, fp);
		} else {
		    err = user_var_localize(arg->upname, fp->name, fp->type);
		}
		if (!err) {
		    maybe_set_arg_const(arg, fp);
		}
	    }
	}

	if (!err && arg->type == GRETL_TYPE_USERIES) {
	    if (fp->type == GRETL_TYPE_LIST) {
		/* FIXME ? */
		arg->upname[0] = '\0';
	    }
	}	
    }

    /* now for any parameters without matching arguments */

    for (i=call->argc; i<fun->n_params && !err; i++) {
	fp = &fun->params[i];
	if (gretl_scalar_type(fp->type)) {
	    err = add_scalar_arg_default(fp);
	} else if (fp->type == GRETL_TYPE_LIST) {
	    err = localize_list(call, NULL, fp, dset);
	}
    }

    if (call->listvars != NULL) {
	if (err) {
	    free(call->listvars);
	    call->listvars = NULL;
	} else {
	    fncall_finalize_listvars(call);
	}
    }

    if (!err) {
	set_listargs_from_call(call, dset);
    }

#if UDEBUG
    fprintf(stderr, "allocate_function args: returning %d\n", err);
#endif

    return err;
}

/**
 * check_function_needs:
 * @dset: pointer to dataset info.
 * @dreq: function data requirements flag.
 * @minver: function minimum program version requirement.
 *
 * Checks whether the requirements in @dreq and @minver
 * are jointly satisfied by the current dataset and gretl
 * version.
 *
 * Returns: 0 if all is OK; E_DATA if the current dataset
 * (or lack thereof) does not satisfy @dreq; 1 otherwise.
 */

int check_function_needs (const DATASET *dset, FuncDataReq dreq,
			  int minver)
{
    static int thisver = 0;

    if (thisver == 0) {
	thisver = gretl_version_number(GRETL_VERSION);
    }

    if (minver > thisver) {
	char vstr[8];

	gretl_version_string(vstr, minver);
	gretl_errmsg_sprintf("This function needs gretl version %s", vstr);
	return 1;
    }

    if ((dset == NULL || dset->v == 0) && dreq != FN_NODATA_OK) {
	gretl_errmsg_set("This function needs a dataset in place");
	return E_DATA;
    }

    if (dreq == FN_NEEDS_TS && !dataset_is_time_series(dset)) {
	gretl_errmsg_set("This function needs time-series data");
	return E_DATA;
    }

    if (dreq == FN_NEEDS_PANEL && !dataset_is_panel(dset)) {
	gretl_errmsg_set("This function needs panel data");
	return E_DATA;
    }

    if (dreq == FN_NEEDS_QM && 
	(!dataset_is_time_series(dset) || 
	 (dset->pd != 4 && dset->pd != 12))) {
	gretl_errmsg_set("This function needs quarterly or monthly data");
	return E_DATA;
    } 

    return 0;
}

/**
 * package_version_ok:
 * @minver: function package minimum program version requirement.
 * @reqstr: location to write required version string (should be
 * at least 8 bytes), or NULL.
 *
 * Returns: 1 if the running version of gretl satisfies the
 * minimum version requirement in @minver; otherwise returns 0,
 * in which case the string representation of @minver is written
 * @reqstr if this is non-NULL.
 */

int package_version_ok (int minver, char *reqstr)
{
    static int thisver = 0;
    int ret = 0;

    if (thisver == 0) {
	thisver = gretl_version_number(GRETL_VERSION);
    }

    ret = thisver >= minver;

    if (!ret && reqstr != NULL) {
	gretl_version_string(reqstr, minver);
    }

    return ret;
}

static int maybe_check_function_needs (const DATASET *dset,
				       const ufunc *fun)
{
    if (fun->pkg == NULL) {
	return 0;
    } else {
	return check_function_needs(dset, fun->pkg->dreq, 
				    fun->pkg->minver);
    }
}

/* next block: handling function return values */

static int handle_scalar_return (const char *vname, void *ptr)
{
    int err = 0;

    if (gretl_is_scalar(vname)) {
	*(double *) ptr = gretl_scalar_get_value(vname, NULL);
    } else {
	*(double *) ptr = NADBL;
	err = E_UNKVAR; /* FIXME */
    }

    return err;
}

static int handle_series_return (const char *vname, void *ptr,
				 DATASET *dset, int copy)
{
    int v = series_index(dset, vname);
    double *x = NULL;
    int err = 0;

    if (!copy && v == 0) {
	copy = 1;
    }
    if (v >= 0 && v < dset->v) {
	if (copy) {
	    x = copyvec(dset->Z[v], dset->n);
	    if (x == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    x = dset->Z[v];
	    dset->Z[v] = NULL;
	}
    } else {
	err = E_UNKVAR;
    }

    *(double **) ptr = x;

    return err;
}

static int handle_matrix_return (const char *name, void *ptr,
				 int copy)
{
    gretl_matrix *ret = NULL;
    int err = 0;

    if (copy) {
	gretl_matrix *m = get_matrix_by_name(name);

	if (m != NULL) {
	    ret = gretl_matrix_copy(m);
	    if (ret == NULL) {
		err = E_ALLOC;
	    }
	} 
    } else {
	ret = steal_matrix_by_name(name);
    }

    if (ret == NULL && !err) {
	err = E_UNKVAR;
    }

    *(gretl_matrix **) ptr = ret;

    return err;
}

static int handle_bundle_return (fncall *call, void *ptr, int copy)
{
    const char *name = call->retname;
    gretl_bundle *ret = NULL;
    int err = 0;

    if (copy) {
	gretl_bundle *b = get_bundle_by_name(name);

	if (b != NULL) {
	    ret = gretl_bundle_copy(b, &err);
	} 
    } else {
	ret = gretl_bundle_pull_from_stack(name, &err);
    }

    if (ret != NULL && call->fun->pkg != NULL &&
	gretl_function_depth() == 1) {
	gretl_bundle_set_creator(ret, call->fun->pkg->name);
    }

    *(gretl_bundle **) ptr = ret;

    return err;
}

static int handle_array_return (fncall *call, void *ptr, int copy)
{
    const char *name = call->retname;
    gretl_array *ret = NULL;
    int err = 0;

    if (copy) {
	gretl_array *a = get_array_by_name(name);

	if (a != NULL) {
	    ret = gretl_array_copy(a, &err);
	} 
    } else {
	ret = gretl_array_pull_from_stack(name, &err);
    }

    *(gretl_array **) ptr = ret;

    return err;
}

static void replace_caller_series (int targ, int src, DATASET *dset)
{
    int t;

    /* replace data values */
    for (t=dset->t1; t<=dset->t2; t++) {
	dset->Z[targ][t] = dset->Z[src][t];
    }

    /* replace variable info? */
    series_set_label(dset, targ, series_get_label(dset, src));
    series_set_display_name(dset, targ, series_get_display_name(dset, src));
}

/* Deal with a list that exists at the level of a user-defined
   function whose execution is now terminating.  Note that this list
   may be the direct return value of the function, or it may have been
   given to the function as an argument, or it may have been constructed
   on the fly. This is flagged by @arg (which will be non-NULL only if
   the list was supplied as a function argument).
*/

static int unlocalize_list (const char *lname, fn_arg *arg,
			    void *ret, DATASET *dset)
{
    int *list = get_list_by_name(lname);
    int d = gretl_function_depth();
    int upd = d - 1;
    int i, vi;

#if UDEBUG
    fprintf(stderr, "\n*** unlocalize_list: '%s', function depth = %d\n", lname, d);
    printlist(list, lname);
    fprintf(stderr, " dset = %p, dset->v = %d\n", (void *) dset, dset->v);
    fprintf(stderr, " list is direct return value? %s\n", arg == NULL ? "yes" : "no");
#endif	    

    if (list == NULL) {
	return E_DATA;
    }

    /* 
       If the list we're looking at was given as a function argument
       we simply shunt all its members to the prior stack level.  But
       if the list is the direct return value from the function, we
       need to overwrite any variables at caller level that have been
       redefined within the function.  If any series have been
       redefined in that way we also need to adjust the list itself,
       replacing the ID numbers of local variables with the
       caller-level IDs -- all supposing the direct return value is
       being assigned by the caller.

       On the other hand, if the list was a temporary construction
       (a list parameter was required, but a single series was given as
       argument), we destroy the list.
    */

    if (arg == NULL && ret == NULL) {
	; /* no-op: we'll trash the list later */
    } else if (arg == NULL) {
	/* handle list as direct return value */
	int j, lev, overwrite;
	const char *vname;

	for (i=1; i<=list[0]; i++) {
	    overwrite = 0;
	    vi = list[i];
	    vname = dset->varname[vi];
	    series_unset_flag(dset, vi, VAR_LISTARG);
	    if (vi > 0 && vi < dset->v && series_get_stack_level(dset, vi) == d) {
		for (j=1; j<dset->v; j++) {
		    lev = series_get_stack_level(dset, j);
		    if (lev == upd && !strcmp(dset->varname[j], vname)) {
			overwrite = 1;
			break;
		    }
		    if (lev == d && j < vi && series_is_listarg(dset, j) &&
			!strcmp(dset->varname[j], vname)) {
			overwrite = 1;
			break;
		    }			
		}
		if (overwrite) {
		    replace_caller_series(j, vi, dset);
		    /* replace ID number in list */
		    list[i] = j;
		} else {
		    series_set_stack_level(dset, vi, upd);
		}
	    }
#if UDEBUG
	    fprintf(stderr, " list-member var %d, '%s': ", vi, vname); 
	    if (overwrite) {
		fprintf(stderr, "found match in caller, overwrote var %d\n", j);
	    } else {
		fprintf(stderr, "no match in caller\n");
	    }
#endif
	}
    } else {
	/* list was given as argument to function */
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi == LISTSEP) {
		continue;
	    }
	    if (series_is_listarg(dset, vi)) {
		series_unset_flag(dset, vi, VAR_LISTARG);
		series_set_stack_level(dset, vi, upd);
	    }
	}
	if (arg->type != GRETL_TYPE_LIST) {
	    /* the list was constructed on the fly */
	    user_var_delete_by_name(lname, NULL);
	}
    }

    return 0;
}

static int handle_string_return (const char *sname, void *ptr)
{
    const char *s = get_string_by_name(sname);
    char *ret = NULL;
    int err = 0;

    if (s == NULL) {
	err = E_DATA;
    } else {
	ret = gretl_strdup(s);
	if (ret == NULL) {
	    err = E_ALLOC;
	}
    }

    *(char **) ptr = ret;

    return err;
}

static void 
maybe_set_return_description (fncall *call, int rtype, 
			      DATASET *dset, 
			      char **descrip)
{
    if (rtype == GRETL_TYPE_SERIES) {
	int v = series_index(dset, call->retname);

	if (v < dset->v) {
	    *descrip = gretl_strdup(series_get_label(dset, v));
	}
    }
}

static int is_pointer_arg (fncall *call, int rtype)
{
    ufunc *u = call->fun;

    if (call->retname != NULL) {
	int i;

	for (i=0; i<call->argc; i++) {
	    if (u->params[i].type == gretl_type_get_ref_type(rtype) &&
		!strcmp(u->params[i].name, call->retname)) {
		return 1;
	    }
	}
    }

    return 0;
}

#define null_return(t) (t == GRETL_TYPE_VOID || t == GRETL_TYPE_NONE)

#define needs_dataset(t) (t == GRETL_TYPE_SERIES || \
			  t == GRETL_TYPE_LIST ||   \
			  t == GRETL_TYPE_SERIES_REF || \
			  t == GRETL_TYPE_LISTS || \
			  t == GRETL_TYPE_LISTS_REF)

static int 
function_assign_returns (fncall *call, int rtype, 
			 DATASET *dset, void *ret,
			 char **descrip, PRN *prn,
			 int *perr)
{
    ufunc *u = call->fun;
    int i, err = 0;

#if UDEBUG
    fprintf(stderr, "function_assign_returns: rtype = %d, call->retname = %s\n", 
	    rtype, call->retname);
#endif

    if (*perr == 0 && !null_return(rtype) && call->retname == NULL) {
	/* missing return value */
	gretl_errmsg_sprintf("Function %s did not provide the specified return value\n"
			     "(expected %s)", u->name, gretl_type_get_name(rtype));
	*perr = err = E_UNKVAR;
    } else if (*perr == 0 && needs_dataset(rtype) && dset == NULL) {
	/* "can't happen" */
	*perr = err = E_DATA;
    }

    if (*perr == 0) {
	/* first we work on the value directly returned by the
	   function (but only if there's no error) 
	*/
	int copy = is_pointer_arg(call, rtype);

	if (rtype == GRETL_TYPE_DOUBLE) {
	    err = handle_scalar_return(call->retname, ret);
	} else if (rtype == GRETL_TYPE_SERIES) {
	    err = handle_series_return(call->retname, ret, dset, copy);
	} else if (rtype == GRETL_TYPE_MATRIX) {
	    err = handle_matrix_return(call->retname, ret, copy);
	} else if (rtype == GRETL_TYPE_LIST) {
	    /* note: in this case the job is finished in
	       stop_fncall(); here we just adjust the info on the
	       listed variables so they don't get deleted
	    */
	    err = unlocalize_list(call->retname, NULL, ret, dset);
	} else if (rtype == GRETL_TYPE_BUNDLE) {
	    err = handle_bundle_return(call, ret, copy);
	} else if (gretl_array_type(rtype)) {
	    err = handle_array_return(call, ret, copy);
	} else if (rtype == GRETL_TYPE_STRING) {
	    err = handle_string_return(call->retname, ret);
	} 

	if (err == E_UNKVAR) {
	    pprintf(prn, "Function %s did not provide the specified return value\n",
		    u->name);
	}

	*perr = err;

	if (!err && dset != NULL && descrip != NULL) {
	    maybe_set_return_description(call, rtype, dset, descrip);
	}
    }

    /* "Indirect return" values and other pointerized args: 
       these should be handled even if the function bombed.
    */

    for (i=0; i<call->argc; i++) {
	fn_arg *arg = &call->args[i];
	fn_param *fp = &u->params[i];
	
	if (needs_dataset(fp->type) && dset == NULL) {
	    err = E_DATA;
	} else if (gretl_ref_type(fp->type)) {
	    if (arg->type == GRETL_TYPE_BUNDLE_REF &&
		arg->upname[0] == '\0') {
		; /* pointer to anonymous bundle: no-op */
	    } else if (arg->type == GRETL_TYPE_SERIES_REF) {
		int v = arg->val.idnum;

		series_decrement_stack_level(dset, v);
		strcpy(dset->varname[v], arg->upname);
	    } else if (gretl_ref_type(arg->type)) {
		user_var_unlocalize(fp->name, arg->upname, fp->type);
	    } 
	} else if (fp->type == GRETL_TYPE_MATRIX && 
		   (fp->flags & ARG_CONST) && 
		   arg->upname[0] != '\0') {
	    /* non-pointerized const matrix argument */
	    user_var *u = get_user_var_by_data(arg->val.m);

	    if (u != NULL) {
		user_var_adjust_level(u, -1);
		user_var_set_name(u, arg->upname);
	    }
	} else if (fp->type == GRETL_TYPE_LIST) {
	    unlocalize_list(fp->name, arg, NULL, dset);
	} 
    }

#if UDEBUG
    fprintf(stderr, "function_assign_returns: returning %d\n", err);
#endif

    return err;
}

/* make a record of the sample information at the time a function is
   called */

static void record_obs_info (obsinfo *o, DATASET *dset)
{
    o->changed = 0;

    if (dset != NULL) {
	o->structure = dset->structure;
	o->pd = dset->pd;
	o->t1 = dset->t1;
	o->t2 = dset->t2;
	strcpy(o->stobs, dset->stobs);
    }
}

/* on function exit, restore the sample information that was in force
   on entry */

static int restore_obs_info (obsinfo *oi, DATASET *dset)
{
    gretlopt opt = OPT_NONE;

    if (oi->structure == CROSS_SECTION) {
	opt = OPT_X;
    } else if (oi->structure == TIME_SERIES) {
	opt = OPT_T;
    } else if (oi->structure == STACKED_TIME_SERIES) {
	opt = OPT_S;
    } else if (oi->structure == SPECIAL_TIME_SERIES) {
	opt = OPT_N;
    } 

    return simple_set_obs(dset, oi->pd, oi->stobs, opt);
}

/* do the basic housekeeping that is required when a function exits:
   destroy local variables, restore previous sample info, etc.
*/

static int stop_fncall (fncall *call, int rtype, void *ret,
			DATASET *dset, PRN *prn, int orig_v,
			int redir_level)
{
    int i, d = gretl_function_depth();
    int delv, anyerr = 0;
    int err = 0;

#if FN_DEBUG
    fprintf(stderr, "stop_fncall: terminating call to "
	    "function '%s' at depth %d, dset->v = %d\n", 
	    call->fun->name, d, (dset != NULL)? dset->v : 0);
#endif

    /* below: delete series local to the function, taking care not to
       delete any that have been "promoted" to caller level via their 
       inclusion in a returned list; also catch any "listarg" series
       that are no longer attached to a list
    */

    if (dset != NULL) {
	for (i=orig_v, delv=0; i<dset->v; i++) {
	    if (series_get_stack_level(dset, i) == d) {
		delv++;
	    }
	}
	if (delv > 0) {
	    if (delv == dset->v - orig_v) {
		/* deleting all added series */
		anyerr = dataset_drop_last_variables(dset, delv);
		if (anyerr && !err) {
		    err = anyerr;
		}
	    } else {
		for (i=dset->v-1; i>=orig_v; i--) {
		    if (series_get_stack_level(dset, i) == d) {
			anyerr = dataset_drop_variable(i, dset);
			if (anyerr && !err) {
			    err = anyerr;
			}
		    } 
		}
	    }    
	}
	if (call->listvars != NULL) {
	    int vi;
	    
	    for (i=1; i<=call->listvars[0]; i++) {
		vi = call->listvars[i];
		series_unset_flag(dset, vi, VAR_LISTARG);
		series_set_stack_level(dset, vi, d - 1);
	    }
	}
    }

    /* direct list return: write the possibly revised list to the
       return pointer.  Note that we can't do this earlier, because
       the ID numbers of the variables in the return list may be
       changed due to the deletion of function-local variables.
    */

    if (!err && rtype == GRETL_TYPE_LIST && ret != NULL) {
	int *lret = gretl_list_copy(get_list_by_name(call->retname));

	if (lret != NULL) {
	    *(int **) ret = lret;
	} else {
	    err = E_ALLOC;
	}
    }

    /* now we're ready to trash function-local vars */
    anyerr = destroy_user_vars_at_level(d);

    if (anyerr && !err) {
	err = anyerr;
    }    

    /* if the function defined a Kalman filter, clean that up */
    delete_kalman(NULL);

    /* if any anonymous equations system was defined: clean up */
    delete_anonymous_equation_system(d);

    pop_program_state();

    if (dset != NULL && call->obs.changed) {
	restore_obs_info(&call->obs, dset);
    }

    set_executing_off(call, dset, prn);

    if (print_redirection_level(prn) > redir_level) {
	gretl_errmsg_set("Incorrect use of 'outfile' in function");
	err = 1;
    }

    return err;
}

#if LOOPSAVE

/* reset any saved "genr" structs in context of saved loops */

static void reset_saved_loops (ufunc *u)
{
    int i;

    for (i=0; i<u->n_lines; i++) {
	if (u->lines[i].loop != NULL) {
# if GLOBAL_TRACE
	    fprintf(stderr, "at exit from %s, reset_saved_loop %p\n",
		    u->name, (void *) u->lines[i].loop);
# endif
	    loop_reset_genrs(u->lines[i].loop);
	}
    }
}

#endif 

static void set_pkgdir (fnpkg *pkg)
{
    const char *p = strrchr(pkg->fname, SLASH);

    if (p != NULL) {
	char *pkgdir = gretl_strndup(pkg->fname, p - pkg->fname);
	
	gretl_insert_builtin_string("pkgdir", pkgdir);
	free(pkgdir);
    }
}

static int start_fncall (fncall *call, DATASET *dset, PRN *prn)
{
    GList *tmp = callstack;
    fncall *prevcall;

    while (tmp != NULL) {
	prevcall = tmp->data;
	if (prevcall->fun == call->fun) {
	    call->recursing = 1;
	    break;
	}
	tmp = tmp->next;
    }    
    
    fn_executing++;
    push_program_state();

    callstack = g_list_append(callstack, call);

#if EXEC_DEBUG
    fprintf(stderr, "start_fncall: added call to %s, depth now %d\n", 
	    call->fun->name, g_list_length(callstack));
#endif

    switch_uservar_hash(fn_executing);

    record_obs_info(&call->obs, dset);

    if (gretl_debugging_on() || call->fun->debug) {
	set_gretl_echo(1);
	set_gretl_messages(1);
	pputs(prn, "*** ");
	bufspace(gretl_function_depth(), prn);
	pprintf(prn, "executing function %s\n", call->fun->name);
    } else {
	set_gretl_echo(0);
	set_gretl_messages(0);
    }

    if (fn_executing == 1 && call->fun->pkg != NULL) {
	set_pkgdir(call->fun->pkg);
    } 

    return 0;
}

static void func_exec_callback (ExecState *s, void *ptr,
				GretlObjType type)
{
    if (s->cmd->ci == FLUSH || is_plotting_command(s->cmd)) {
	/* we permit "reach-back" into the GUI for these */
	EXEC_CALLBACK gc = get_gui_callback();

	if (gc != NULL) {
	    gc(s, NULL, 0);
	}
    }
}

static double arg_get_double_val (fn_arg *arg)
{
    if (gretl_scalar_type(arg->type)) {
	return arg->val.x;
    } else if (arg->type == GRETL_TYPE_MATRIX) {
	return arg->val.m->val[0];
    } else {
	return NADBL;
    }
}

static int check_function_args (fncall *call, PRN *prn)
{
    ufunc *u = call->fun;    
    fn_arg *arg;
    fn_param *fp;
    double x;
    int i, err = 0;

    for (i=0; i<call->argc && !err; i++) {
	arg = &call->args[i];
	fp = &u->params[i];

	if ((fp->flags & ARG_OPTIONAL) && arg->type == GRETL_TYPE_NONE) {
	    ; /* this is OK */
	} else if (gretl_scalar_type(fp->type) && arg->type == GRETL_TYPE_DOUBLE) {
	    ; /* OK */
	} else if (fp->type == GRETL_TYPE_SERIES && arg->type == GRETL_TYPE_USERIES) {
	    ; /* OK */
	} else if (gretl_scalar_type(fp->type) && 
		   arg->type == GRETL_TYPE_MATRIX &&
		   gretl_matrix_is_scalar(arg->val.m)) {
	    ; /* OK */
	} else if (fp->type == GRETL_TYPE_LIST && arg->type == GRETL_TYPE_USERIES) {
	    ; /* OK */
	} else if (fp->type == GRETL_TYPE_LIST && arg->type == GRETL_TYPE_NONE) {
	    ; /* OK ("null" was passed as argument) */
	} else if (fp->type != arg->type) {
	    pprintf(prn, _("%s: argument %d is of the wrong type (is %s, should be %s)\n"), 
		    u->name, i + 1, gretl_type_get_name(arg->type), 
		    gretl_type_get_name(fp->type));
	    err = E_TYPES;
	}

	if (!err && fp->type == GRETL_TYPE_DOUBLE) {
	    x = arg_get_double_val(arg);
	    if ((!na(fp->min) && x < fp->min) ||
		(!na(fp->max) && x > fp->max)) {
		pprintf(prn, _("%s, argument %d: value %g is out of bounds\n"), 
			u->name, i + 1, x);
		err = E_INVARG;
	    }
	}
    }

    for (i=call->argc; i<u->n_params && !err; i++) {
	/* do we have defaults for any empty args? */
	fp = &u->params[i];
	if (!(fp->flags & ARG_OPTIONAL) && no_scalar_default(fp)) {
	    pprintf(prn, _("%s: not enough arguments\n"), u->name);
	    err = E_ARGS;
	}
    }

#if UDEBUG
    if (err) {
	fprintf(stderr, "CHECK_FUNCTION_ARGS: err = %d\n", err);
    }
#endif    

    return err;
}

/**
 * user_function_set_debug:
 * @name: the name of the function.
 * @debug: boolean, if 1 then start debugging function, if
 * 0 then stop debugging.
 *
 * Enables or disables debugging for a user-defined function.
 *
 * Returns: 0 on success, non-zero if no function is specified
 * or if the function is not found.
 */

int user_function_set_debug (const char *name, int debug)
{
    ufunc *fun;

    if (name == NULL || *name == '\0') {
	return E_ARGS;
    } 

    fun = get_user_function_by_name(name);

    if (fun == NULL) {
	return E_UNKVAR;
    } else {
	fun->debug = debug;
	return 0;
    }
}

#define debug_cont(c) (c->ci == FUNDEBUG && (c->opt & OPT_C))
#define debug_next(c) (c->ci == FUNDEBUG && (c->opt & OPT_N))

#define set_debug_cont(c) (c->ci = FUNDEBUG, c->opt = OPT_C)
#define set_debug_next(c) (c->ci = FUNDEBUG, c->opt = OPT_N)

/* loop for stepping through function commands; returns
   1 if we're debugging, otherwise 0 
*/

static int debug_command_loop (ExecState *state, 
			       DATASET *dset,
			       DEBUG_READLINE get_line,
			       DEBUG_OUTPUT put_func,
			       int *errp)
{
    int brk = 0, err = 0;

    state->flags |= DEBUG_EXEC;

    while (!brk) {
#if DDEBUG
	fprintf(stderr, "--- debug_command_loop calling get_line\n"); 
#endif
	*errp = (*get_line)(state);
	if (*errp) {
	    /* the debugger failed: get out right away */
	    state->flags &= ~DEBUG_EXEC;
	    return 0;
	}

	err = parse_command_line(state->line, state->cmd, 
				 dset, NULL);
	if (err) {
	    if (!strcmp(state->line, "c")) {
		/* short for 'continue' */
		set_debug_cont(state->cmd);
		err = 0;
	    } else if (!strcmp(state->line, "n")) {
		/* short for 'next' */
		set_debug_next(state->cmd);
		err = 0;
	    }
	}

	if (err || debug_cont(state->cmd) || debug_next(state->cmd)) {
	    brk = 1;
	} else {
	    /* execute interpolated command */
	    err = gretl_cmd_exec(state, dset);
	    if (put_func != NULL) {
#if DDEBUG		
		fprintf(stderr, "--- debug_command_loop calling put_func\n"); 
#endif
		(*put_func)(state);
	    }		
	}
    }

    if (!debug_next(state->cmd)) {
	state->flags &= ~DEBUG_EXEC;
    }

    return 1 + debug_next(state->cmd);
}

static const char *get_funcerr_message (ExecState *state)
{
    const char *p = state->cmd->param;

    if (p != NULL && strlen(p) == gretl_namechar_spn(p)) {
	const char *s = get_string_by_name(p);

	if (s != NULL) {
	    return s;
	}
    } else if (p == NULL) {
	return "none";
    }

    return p;
}

static void set_func_error_message (int err, ufunc *u, 
				    ExecState *state,
				    const char *line,
				    int lineno)
{
    if (err == E_FUNCERR) {
	/* let the function writer set the message */
	gretl_errmsg_sprintf(_("Error message from %s():\n %s"), u->name,
			     get_funcerr_message(state));
    } else if (err == E_STOP) {
	; /* no-op */
    } else {
	/* we'll handle this here */
	const char *msg = gretl_errmsg_get();

	if (*msg == '\0') {
	    msg = errmsg_get_with_default(err);
	    gretl_errmsg_set(msg);
	}

	if (lineno < 0) {
	    /* indicates that we don't have a real line number
	       available (the error occurred inside a loop)
	    */
	    if (*line != '\0' && err != E_FNEST) {
		gretl_errmsg_sprintf("*** error within loop in function %s\n> %s", 
				     u->name, line);
	    } else {
		gretl_errmsg_sprintf("*** error within loop in function %s\n", 
				     u->name, lineno);
	    }
	} else {
	    if (*line != '\0' && err != E_FNEST) {
		gretl_errmsg_sprintf("*** error in function %s, line %d\n> %s", 
				     u->name, lineno, line);
	    } else {
		gretl_errmsg_sprintf("*** error in function %s, line %d\n", 
				     u->name, lineno);
	    }
	}

	if (gretl_function_depth() > 1) {
	    GList *tmp = g_list_last(callstack);

	    if (tmp != NULL) {
		tmp = g_list_previous(tmp);
		if (tmp != NULL) {
		    fncall *call = tmp->data;

		    gretl_errmsg_sprintf(" called by function %s", call->fun->name);
		}
	    }
	}
    }
}

static int handle_return_statement (fncall *call,
				    ExecState *state,
				    DATASET *dset,
				    int lineno)
{
    const char *s = state->line + 6; /* skip "return" */
    ufunc *fun = call->fun;
    int err = 0;

#if EXEC_DEBUG
    fprintf(stderr, "%s: return: s = '%s'\n", fun->name, s);
#endif

    s += strspn(s, " ");

    if (fun->rettype == GRETL_TYPE_VOID) {
	if (*s == '\0') {
	    ; /* plain "return" from void function: OK */
	} else if (lineno < 0) {
	    gretl_errmsg_sprintf("%s: non-null return value '%s' is not valid", 
				 fun->name, s);
	    err = E_TYPES;
	} else {
	    gretl_errmsg_sprintf("%s, line %d: non-null return value '%s' is not valid", 
				 fun->name, lineno, s);
	    err = E_TYPES;
	}
	return err; /* handled */
    }

    if (*s == '\0') {
	if (lineno < 0) {
	    gretl_errmsg_sprintf("%s: return value is missing", fun->name);
	} else {
	    gretl_errmsg_sprintf("%s, line %d: return value is missing",
				 fun->name, lineno);
	}	    
	err = E_TYPES;
    } else {
	int len = gretl_namechar_spn(s);

	if (len == strlen(s)) {
	    /* returning a named variable */
	    call->retname = gretl_strndup(s, len);
	} else {
	    char formula[MAXLINE];
	    
	    sprintf(formula, "$retval=%s", s);
	    err = generate(formula, dset, fun->rettype, OPT_P, NULL);
	    if (err) {
		set_func_error_message(err, fun, state, s, lineno);
	    } else {
		call->retname = gretl_strdup("$retval");
	    }
	}
    }

    if (!err && call->retname == NULL) {
	err = E_ALLOC;
    }

    return err;
}

static int do_debugging (ExecState *s)
{
    return s->cmd->ci > 0 && !s->cmd->context &&
	!gretl_compiling_loop();
}

#define CDEBUG 1

/* Construct bundle of arguments to send to C function;
   load plugin; send arguments; and try to retrieve
   return value, if any.
*/

static int handle_plugin_call (fncall *call,
			       DATASET *dset,
			       void *retval,
			       PRN *prn)
{
    ufunc *u = call->fun;
    int (*cfunc) (gretl_bundle *, PRN *);
    void *handle;
    gretl_bundle *argb;
    int i, err = 0;

#if CDEBUG
    fprintf(stderr, "handle_plugin_call: function '%s'\n", u->name);
#endif

    if (u->pkg == NULL) {
	return E_DATA;
    }

    cfunc = get_packaged_C_function(u->pkg->name, u->name, &handle);

    if (cfunc == NULL) {
	fputs(I_("Couldn't load plugin function\n"), stderr);
	return E_FOPEN;
    }    

    /* bundle to serve as wrapper for arguments */
    argb = gretl_bundle_new();

    if (argb == NULL) {
	close_plugin(handle);
	return E_ALLOC;
    }

    for (i=0; i<call->argc && !err; i++) {
	fn_param *fp = &u->params[i];
	fn_arg *arg = &call->args[i];
	const char *key = fp->name;

	if (!type_can_be_bundled(fp->type)) {
	    fprintf(stderr, "type %d: cannot be bundled\n", fp->type);
	    err = E_TYPES;
	    break;
	}

	if (gretl_scalar_type(fp->type)) {
	    double x;

	    if (arg->type == GRETL_TYPE_NONE) {
		x = fp->deflt;
	    } else if (arg->type == GRETL_TYPE_MATRIX) {
		x = arg->val.m->val[0];
	    } else {
		x = arg->val.x;
	    }
	    if (fp->type == GRETL_TYPE_INT || fp->type == GRETL_TYPE_OBS) {
		x = (int) x;
	    } else if (fp->type == GRETL_TYPE_BOOL) {
		x = (x != 0.0);
	    }
	    err = gretl_bundle_set_data(argb, key, &x, GRETL_TYPE_DOUBLE, 0);
	} else if (fp->type == GRETL_TYPE_MATRIX) {
	    gretl_matrix *m = arg->val.m;
	    
	    err = gretl_bundle_set_data(argb, key, m, fp->type, 0);
	} else if (fp->type == GRETL_TYPE_SERIES) {
	    int size = sample_size(dset);
	    double *px;

	    if (arg->type == GRETL_TYPE_USERIES) {
		px = dset->Z[arg->val.idnum];
	    } else {
		px = arg->val.px;
	    }
	    err = gretl_bundle_set_data(argb, key, px + dset->t1, 
					GRETL_TYPE_SERIES, size);
	} else if (fp->type == GRETL_TYPE_BUNDLE) {
	    gretl_bundle *b = arg->val.b;

	    err = gretl_bundle_set_data(argb, key, b, 
					GRETL_TYPE_BUNDLE, 0);
	} else {
	    /* FIXME strings and maybe other types */
	    err = E_TYPES;
	}
#if CDEBUG
	fprintf(stderr, "arg[%d] (\"%s\") type = %d, err = %d\n", 
		i, key, fp->type, err);
#endif
    }

    if (!err) {
	/* call the plugin function */
	err = (*cfunc) (argb, prn);
    }

    close_plugin(handle);

    if (!err && u->rettype != GRETL_TYPE_VOID && retval != NULL) {
	GretlType type;
	int size;
	void *ptr;

	ptr = gretl_bundle_steal_data(argb, "retval", &type, &size, &err);

	if (!err) {
#if CDEBUG
	    fprintf(stderr, "%s: stole return value of type %d\n", 
		    u->name, type);
#endif
	    if (type != u->rettype) {
		fprintf(stderr, "handle_plugin_call: type doesn't match u->rettype\n");
		err = E_TYPES;
	    }
	}
	if (!err) {
	    if (type == GRETL_TYPE_DOUBLE) {
		double *px = ptr;
		
		*(double *) retval = *px;
		free(ptr);
	    } else if (type == GRETL_TYPE_MATRIX) {
		*(gretl_matrix **) retval = ptr;
	    } else if (type == GRETL_TYPE_BUNDLE) {
		*(gretl_bundle **) retval = ptr;
	    } else if (gretl_array_type(type)) {
		*(gretl_array **) retval = ptr;
	    }
	}
    }

    /* free the arguments wrapper */
    gretl_bundle_destroy(argb);

    return err;
}

int current_function_size (void)
{
    ufunc *u = currently_called_function();

    return (u != NULL)? u->n_lines : 0;
}

#if LOOPSAVE

int attach_loop_to_function (void *ptr)
{
    fncall *call = current_function_call();
    ufunc *u = NULL;
    int err = 0;

    if (call != NULL && !call->recursing) {
	u = call->fun;
    }

    if (u == NULL) {
	return E_DATA;
    } else {
	if (u->lines[u->line_idx].loop != NULL) {
	    /* there should not already be a loop attached */
	    err = E_DATA;
	} else {
	    /* OK, attach it */
	    u->lines[u->line_idx].loop = ptr;
#if LSDEBUG || GLOBAL_TRACE
	    fprintf(stderr, "attaching loop at %p to function %s, line %d\n",
		    ptr, u->name, u->line_idx);
#endif
	}
    }

    return err;
}

int detach_loop_from_function (void *ptr)
{
    fncall *call = current_function_call();
    ufunc *u = NULL;
    int err = 0;

    if (call != NULL) {
	u = call->fun;
    }

    if (u == NULL) {
	return E_DATA;
    } else {
	int i;

	for (i=0; i<u->n_lines; i++) {
	    if (u->lines[i].loop == ptr) {
#if LSDEBUG || GLOBAL_TRACE
		fprintf(stderr, "detaching loop at %p from function %s\n",
			ptr, u->name);
#endif
		u->lines[i].loop = NULL;
		break;
	    }
	}
    }

    return err;
}

#endif /* LOOPSAVE */

static gchar *return_line;

/* to be called when a "return" statement occurs within a
   loop that's being called by a function
*/

int set_function_should_return (const char *line)
{
    if (gretl_function_depth() > 0) {
	g_free(return_line);
	return_line = g_strdup(line);
	return 0;
    } else {
	gretl_errmsg_set("return: can only be used in a function");
	return E_PARSE;
    }
}

static int get_return_line (ExecState *state)
{
    if (return_line == NULL) {
	return 0;
    } else {
	strcpy(state->line, return_line);
	g_free(return_line);
	return_line = NULL;
	return 1;
    }
}

int gretl_function_exec (ufunc *u, int rtype, DATASET *dset,
			 void *ret, char **descrip, PRN *prn)
{
    DEBUG_READLINE get_line = NULL;
    DEBUG_OUTPUT put_func = NULL;
    ExecState state;
    fncall *call = NULL;
    MODEL *model = NULL;
    char line[MAXLINE];
    CMD cmd;
    int orig_v = 0;
    int orig_n = 0;
    int orig_t1 = 0;
    int orig_t2 = 0;
    int indent0 = 0, started = 0;
    int redir_level = 0;
    int retline = -1;
    int debugging = u->debug;
#if LOOPSAVE
    int loopstart = 0;
#endif
    int i, err = 0;

#if EXEC_DEBUG || GLOBAL_TRACE
    fprintf(stderr, "gretl_function_exec: starting %s\n", u->name);
#endif

    err = maybe_check_function_needs(dset, u);
    if (err) {
	function_clear_args(u);
	return err;
    }

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: argc = %d\n", u->argc);
    fprintf(stderr, "u->n_params = %d\n", u->n_params);
#endif

    if (dset != NULL) {
	orig_v = dset->v;
	orig_n = dset->n;
	orig_t1 = dset->t1;
	orig_t2 = dset->t2;
    }

    if (!function_is_plugin(u)) {
	/* precaution */
	function_state_init(&cmd, &state, &indent0);
	call = fncall_new(u);
	if (call == NULL) {
	    fprintf(stderr, "fncall_new() returned NULL\n");
	    return E_ALLOC;
	}
    }

    err = check_function_args(call, prn);

    if (function_is_plugin(u)) {
	if (!err) {
	    err = handle_plugin_call(call, dset, ret, prn);
	}
	fncall_free(call);
	return err;
    }

    if (!err) {
	err = allocate_function_args(call, dset);
    }

    if (err) {
	/* get out before allocating further storage */
	fncall_free(call);
	return err;
    }  

    model = allocate_working_model();
    if (model == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_cmd_init(&cmd);
    }

    if (!err) {
	*line = '\0';
	gretl_exec_state_init(&state, FUNCTION_EXEC, line, &cmd, 
			      model, prn);
	if (dset != NULL) {
	    if (dset->submask != NULL) {
		state.submask = copy_datainfo_submask(dset, &err);
	    }
	    state.padded = dset->padmask != NULL;
	}
	state.callback = func_exec_callback;
    }

    if (!err) {
	err = start_fncall(call, dset, prn);
	if (!err) {
	    started = 1;
	    redir_level = print_redirection_level(prn);
	}
    }

#if EXEC_DEBUG
    fprintf(stderr, "start_fncall: err = %d\n", err);
#endif

    if (debugging) {
	get_line = get_debug_read_func();
	put_func = get_debug_output_func();
	if (get_line != NULL) {
	    debugging = 2;
	}
    }

    /* get function lines in sequence and check, parse, execute */

    for (i=0; i<u->n_lines && !err; i++) {
	int lineno = u->lines[i].idx;

	strcpy(line, u->lines[i].s);

	if (debugging) {
	    pprintf(prn, "%s> %s\n", u->name, line);
	} else if (gretl_echo_on()) {
	    pprintf(prn, "? %s\n", line);
	}

	if (u->lines[i].ignore) {
	    continue;
	}

#if LOOPSAVE
	if (u->lines[i].loop != NULL && !call->recursing) {
# if LSDEBUG	    
	    fprintf(stderr, "%s: got loop %p on line %d (%s)\n", u->name,
		    (void *) u->lines[i].loop, i, line);
# endif
	    if (!gretl_if_state_false()) {
		/* not blocked, so execute the loop code */
		err = gretl_loop_exec(&state, dset, u->lines[i].loop);
		if (err) {
		    set_func_error_message(err, u, &state, state.line, -1);
		    break;
		} else if (get_return_line(&state)) {
		    /* a "return" statement was encountered in the loop */
		    err = handle_return_statement(call, &state, dset, -1);
		    break;
		}
	    }
	    /* skip to the matching 'endloop' */
	    i = u->lines[i].next_idx;
	    continue;
	} else {
	    err = maybe_exec_line(&state, dset, &loopstart);
	    if (loopstart) {
		u->line_idx = i;
		loopstart = 0;
	    }
	}
#else
	err = maybe_exec_line(&state, dset, NULL);
#endif /* LOOPSAVE or not */

	if (!err && !gretl_compiling_loop() && state.cmd->ci == FUNCRET) {
	    err = handle_return_statement(call, &state, dset, lineno);
	    if (i < u->n_lines - 1) {
		retline = i;
	    }
	    break;
	}

	if (debugging > 1 && do_debugging(&state)) {
	    if (put_func != NULL) {
		pprintf(prn, "-- debugging %s, line %d --\n", u->name, lineno);
		(*put_func)(&state);
	    } else {
		pprintf(prn, "-- debugging %s, line %d --\n", u->name, lineno);
	    }
	    debugging = debug_command_loop(&state, dset,
					   get_line, put_func, 
					   &err);
	}

	if (err) {
	    set_func_error_message(err, u, &state, line, lineno);
	    break;
	}

	if (state.cmd->ci == SETOBS) {
	    /* set flag for reverting on exit */
	    call->obs.changed = 1;
	}

	if (gretl_execute_loop()) {
#if LOOPSAVE
	    /* mark the ending point of an (outer) loop */
	    u->lines[u->line_idx].next_idx = i;
#endif
	    err = gretl_loop_exec(&state, dset, NULL);
	    if (err) {
		/* note that @lineno will point at the end of a loop here */
		set_func_error_message(err, u, &state, state.line, -1);
	    }
	    if (get_return_line(&state)) {
		/* a "return" statement was encountered in the loop */
		err = handle_return_statement(call, &state, dset, -1);
		break;
	    }
	}
    }

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: %s: finished main exec, "
	    "err = %d, dset->v = %d\n", u->name, err, 
	    (dset != NULL)? dset->v : 0);
#endif

    if (dset != NULL) {
	/* restore the sample that was in place on entry */
	if (dataset_is_complex_subsampled(dset)) {
	    if (state.submask == NULL) {
		/* we were not sub-sampled on entry */
		restore_full_sample(dset, NULL);
	    } else if (submask_cmp(state.submask, dset->submask)) {
		/* we were sub-sampled differently on entry */
		gretlopt opt = state.padded ? OPT_B : OPT_NONE;

		restore_full_sample(dset, NULL);
		restrict_sample_from_mask(state.submask, dset, opt);
	    } 
	} else if (dset->n > orig_n) {
	    /* some observations were added inside the function; note
	       that this is not allowed if the dataset is subsampled
	       on entry
	    */
	    dataset_drop_observations(dset, dset->n - orig_n);
	}
	dset->t1 = orig_t1;
	dset->t2 = orig_t2;
    }

    if (err) {
	if (gretl_function_depth() == 1) {
	    gretl_if_state_clear();
	} else {
	    gretl_if_state_reset(indent0);
	}
    } else if (retline >= 0) {
	/* we returned prior to the end of the function */
	gretl_if_state_reset(indent0);
    } else {
	err = gretl_if_state_check(indent0);
    }

    function_assign_returns(call, rtype, dset, ret, 
			    descrip, prn, &err);

#if LOOPSAVE
    if (!err && !call->recursing) {
	reset_saved_loops(call->fun);
    }
#endif    

    gretl_exec_state_clear(&state);

    if (started) {
	int stoperr = stop_fncall(call, rtype, ret, dset, prn,
				  orig_v, redir_level);

	if (stoperr && !err) {
	    err = stoperr;
	}
    } else {
	fncall_free(call);
    }

#if EXEC_DEBUG || GLOBAL_TRACE
    fprintf(stderr, "gretl_function_exec (%s) finished: err = %d\n",
	    u->name, err);
#endif

    return err;    
}

/* look up name of supplied argument based on name of variable
   inside function */

char *gretl_func_get_arg_name (const char *argvar, int *err)
{
    fncall *call = current_function_call();
    char *ret = NULL;

    *err = E_DATA;

    if (call != NULL) {
	ufunc *u = call->fun;
	int i, n = call->argc;

	for (i=0; i<n; i++) {
	    if (!strcmp(argvar, u->params[i].name)) {
		*err = 0;
		ret = gretl_strdup(call->args[i].upname); 
		if (ret == NULL) {
		    *err = E_ALLOC;
		}
		break;
	    }
	}
    }

    return ret;
}

/**
 * object_is_const:
 * @name: name of object (e.g. matrix).
 *
 * Checks whether the named object currently has 'const' status,
 * by virtue of its being made available as a const argument
 * to a user-defined function.
 *
 * Returns: non-zero if the object is const, 0 if it is not.
 */

int object_is_const (const char *name)
{
    fncall *call = current_function_call();

    if (call != NULL) {
	int i, n = call->argc;

	for (i=0; i<n; i++) {
	    const char *aname = call->args[i].name;

	    if (aname != NULL && !strcmp(name, aname)) {
		return call->args[i].flags & ARG_CONST;
	    }
	}
    }

    return 0;
}

/**
 * object_is_function_arg:
 * @name: name of object (e.g. matrix).
 *
 * Checks whether the named object has been made available
 * as a function argument.
 *
 * Returns: 1 if the object has arg status, 0 otherwise.
 */

int object_is_function_arg (const char *name)
{
    fncall *call = current_function_call();

    if (call != NULL) {
	int i, n = call->argc;

	for (i=0; i<n; i++) {
	    const char *aname = call->args[i].name;

	    if (aname != NULL && !strcmp(name, aname)) {
		return 1;
	    }
	}
    }

    return 0;
}

/**
 * sample_range_get_extrema:
 * @dset: dataset info.
 * @t1: location to receive earliest possible starting observation.
 * @t2: location to receive latest possible ending observation.
 *
 * Fills out @t1 and @t2, making allowance for the possibility
 * that we're currently executing a function, on entry to 
 * which the sample range was restricted: within the function,
 * we are not allowed to overstep the bounds set on entry. 
 */

void sample_range_get_extrema (const DATASET *dset, int *t1, int *t2)
{
    fncall *call = current_function_call();

    if (call != NULL) {
	*t1 = call->obs.t1;
	*t2 = call->obs.t2;
    } else {
	*t1 = 0;
	*t2 = dset->n - 1;
    }
}

/**
 * gretl_functions_cleanup:
 * 
 * For internal use: frees all resources associated with
 * user-defined functions and function packages.
 */

void gretl_functions_cleanup (void)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	ufunc_free(ufuns[i]);
    }

    free(ufuns);
    ufuns = NULL;
    n_ufuns = 0;

    for (i=0; i<n_pkgs; i++) {
	function_package_free(pkgs[i]);
    }

    free(pkgs);
    pkgs = NULL;
    n_pkgs = 0;
}

/* generate help output for a packaged function, either for
   display on the console, or with markup for display in a
   GtkTextView window in the GUI (opt & OPT_M)
*/

static void real_user_function_help (ufunc *fun, gretlopt opt, PRN *prn)
{
    fnpkg *pkg = fun->pkg;
    int markup = (opt & OPT_M);
    int i;

    if (markup) {
	pprintf(prn, "<@itl=\"Function\">: %s\n", fun->name);
    } else {
	pprintf(prn, "Function: %s\n", fun->name);
    }

    if (pkg != NULL) {
	if (markup) {
	    pprintf(prn, "<@itl=\"Package\">: %s %s (%s)\n", pkg->name, pkg->version,
		    pkg->date);
	    pprintf(prn, "<@itl=\"Author\">: %s\n", pkg->author? pkg->author : "unknown");
	    if (pkg->email != NULL && *pkg->email != '\0') {
		pprintf(prn, "<@itl=\"Email\">: %s\n", pkg->email);
	    }
	} else {
	    pprintf(prn, "Package: %s %s (%s)\n", pkg->name, pkg->version,
		    pkg->date);
	    pprintf(prn, "Author: %s\n", pkg->author? pkg->author : "unknown");
	    if (pkg->email != NULL && *pkg->email != '\0') {
		pprintf(prn, "Email:  %s\n", pkg->email);
	    }
	}
	pputc(prn, '\n');
    }

    if ((opt & OPT_G) && pkg != NULL && pkg->gui_help != NULL) {
	/* GUI-specific help is preferred, and is available */
	if (markup) {
	    pputs(prn, "<@itl=\"Help text\">:\n\n");
	    pputs(prn, "<mono>\n");
	} else {
	    pputs(prn, "Help text:\n");
	}   	
	pputs(prn, pkg->gui_help);
	if (markup) {
	    pputs(prn, "\n</mono>");
	}
	pputs(prn, "\n\n");
	return;
    }

    if (markup) {
	pputs(prn, "<@itl=\"Parameters\">: ");
    } else {
	pputs(prn, "Parameters: ");
    }    

    if (fun->n_params > 0) {
	pputc(prn, '\n');
	for (i=0; i<fun->n_params; i++) {
	    pprintf(prn, " %s (%s", 
		    fun->params[i].name, gretl_type_get_name(fun->params[i].type));
	    if (fun->params[i].descrip != NULL) {
		pprintf(prn, ": %s)\n", fun->params[i].descrip);
	    } else {
		pputs(prn, ")\n");
	    }
	}
	pputc(prn, '\n');
    } else {
	pputs(prn, "none\n\n");
    }

    if (markup) {
	pputs(prn, "<@itl=\"Return value\">: ");
    } else {
	pputs(prn, "Return value: ");
    }      

    if (fun->rettype != GRETL_TYPE_NONE && fun->rettype != GRETL_TYPE_VOID) {
	pprintf(prn, "%s\n\n", gretl_type_get_name(fun->rettype));
    } else {
	pputs(prn, "none\n\n");
    }

    if (pkg != NULL && pkg->help != NULL) {
	if (markup) {
	    pputs(prn, "<@itl=\"Help text\">:\n\n");
	    pputs(prn, "<mono>\n");
	} else {
	    pputs(prn, "Help text:\n");
	}   	
	if (has_suffix(pkg->help, ".pdf")) {
	    const char *s = strrchr(pkg->help, ':');

	    if (s != NULL) {
		if (markup) {
		    pprintf(prn, "See <@pdf=\"%s\">", s + 1);
		} else {
		    pprintf(prn, "See %s", s + 1);
		}
	    } else {
		pputs(prn, pkg->help);
	    }
	} else {	
	    pputs(prn, pkg->help);
	}
	if (markup) {
	    pputs(prn, "\n</mono>");
	}
	pputs(prn, "\n\n");
    }

    if (pkg != NULL && pkg->sample != NULL) {
	if (markup) {
	    pputs(prn, "<@itl=\"Sample script\">:\n\n");
	    pputs(prn, "<code>\n");
	} else {
	    pputs(prn, "Sample script:\n\n");
	}
	pputs(prn, pkg->sample);
	if (markup) {
	    pputs(prn, "\n</code>\n");
	} else {
	    pprintf(prn, "\n\n");
	}
    }	
}

/**
 * user_function_help:
 * @fnname: name of function.
 * @opt: may include OPT_M for adding markup, OPT_G
 * for preferring GUI-specific help, if available.
 * @prn: printing struct.
 * 
 * Looks for a function named @fnname and prints
 * as much help information as possible.  
 *
 * Returns: 0 on success, non-zero if the function is not
 * found.
 */

int user_function_help (const char *fnname, gretlopt opt, PRN *prn)
{
    ufunc *fun = get_user_function_by_name(fnname);
    int err = 0;

    if (fun == NULL) {
	pprintf(prn, _("\"%s\" is not defined.\n"), fnname);
	err = 1;
    } else {
	real_user_function_help(fun, opt, prn);
    }

    return err;
}

/**
 * function_package_has_PDF_doc:
 * @pkg: function-package pointer.
 * @pdfname: location to receive basename of PDF file,
 * if applicable (or NULL if this is not wanted).
 * 
 * Checks whether @pkg is documented in the form of a PDF file.
 *
 * Returns: 1 if so (in which case the name of that file is
 * returned via @pdfname), 0 otherwise.
 */

int function_package_has_PDF_doc (fnpkg *pkg, char **pdfname)
{
    int ret = 0;

    if (pkg->help != NULL && !strncmp(pkg->help, "pdfdoc", 6)) {
	ret = 1;
	if (pdfname != NULL) {
	    *pdfname = switch_ext_new(pkg->fname, "pdf");
	    if (*pdfname == NULL) {
		ret = 0;
	    }
	}
    }

    return ret;
}

/**
 * function_package_has_gui_help:
 * @pkg: function-package pointer.
 * 
 * Checks whether @pkg has GUI-specific help text.
 *
 * Returns: 1 if so, 0 otherwise.
 */

int function_package_has_gui_help (fnpkg *pkg)
{
    return pkg->gui_help != NULL && !string_is_blank(pkg->gui_help);
}

/**
 * delete_function_package:
 * @gfnname: full path to gfn file.
 * 
 * Deletes @gfnname, taking care of deleting its enclosing
 * specific subdirectory and all of its contents if applicable.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int delete_function_package (const char *gfnname)
{
    char *p = strrchr(gfnname, SLASH);
    gchar *pkgname = NULL;
    gchar *pkgdir = NULL;
    gchar *pkgsub = NULL;
    int err = 0;

    if (p != NULL) {
	pkgname = g_strdup(p + 1);
	p = strrchr(pkgname, '.');
	if (p != NULL) {
	    *p = '\0';
	}
	pkgdir = g_strdup(gfnname);
	p = strrchr(pkgdir, SLASH);
	*p = '\0';
	p = strrchr(pkgdir, SLASH);
	if (p != NULL) {
	    pkgsub = g_strdup(p + 1);
	}
    }

    if (pkgname != NULL && pkgdir != NULL && pkgsub != NULL &&
	!strcmp(pkgname, pkgsub)) {
	err = gretl_deltree(pkgdir);
    } else {
	/* should not happen! */
	err = gretl_remove(gfnname);
    }

    g_free(pkgname);
    g_free(pkgdir);
    g_free(pkgsub);

    return err;
}

/**
 * uninstall_function_package:
 * @package: name of package (with or without .gfn extension).
 * @opt: option flag(s).
 * @prn: gretl printer.
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int uninstall_function_package (const char *package, gretlopt opt,
				PRN *prn)
{
    gchar *gfnname = NULL;
    gchar *pkgname = NULL;
    char *p, fname[MAXLEN];
    int err;

    if (has_suffix(package, ".gfn")) {
	gfnname = g_strdup(package);
	pkgname = g_strdup(package);
	p = strrchr(pkgname, '.');
	*p = '\0';
    } else {
	gfnname = g_strdup_printf("%s.gfn", package);
	pkgname = g_strdup(package);
    }

    *fname = '\0';
    err = get_full_filename(gfnname, fname, OPT_I);

    if (!err && !gretl_file_exists(fname)) {
	gretl_errmsg_sprintf("Couldn't find %s", gfnname);
	err = E_FOPEN;
    }

    if (!err) {
	function_package_unload_full_by_filename(fname);
	if (opt & OPT_P) {
	    err = delete_function_package(fname);
	}	
    }

    if (!err && gretl_messages_on()) {
	if (opt & OPT_P) {
	    pprintf(prn, "Purged %s\n", pkgname);
	} else {
	    pprintf(prn, "Removed %s\n", pkgname);
	}
    }

    g_free(gfnname);
    g_free(pkgname);

    return err;
}
