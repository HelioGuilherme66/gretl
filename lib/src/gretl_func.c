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

#include "libgretl.h"
#include "gretl_func.h"

#define CALLSTACK_DEPTH 8

typedef struct ufunc_ ufunc;
typedef struct fncall_ fncall;

struct ufunc_ {
    char name[32];
    int n_lines;
    char **lines;
};

struct fncall_ {
    ufunc *fun;
    int argc;
    char **argv;
};

static int n_ufuns;
static ufunc **ufuns;

static fncall **callstack;

static void free_fncall (fncall *call);

/* record of state, and communication of state with outside world */

static int compiling;
static int executing;

int gretl_compiling_function (void)
{
    return compiling;
}

static void set_compiling_on (void)
{
    compiling = 1;
}

static void set_compiling_off (void)
{
    compiling = 0;
}

int gretl_executing_function (void)
{
    return executing;
}

static void set_executing_on (void)
{
    executing++;
}

static void set_executing_off (void)
{
    executing--;
}

/* function call stack mechanism */

static int callstack_init (void)
{
    int i, err = 0;

    if (callstack != NULL) {
	return 0;
    }

    callstack = malloc(CALLSTACK_DEPTH * sizeof *callstack);
    if (callstack == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<CALLSTACK_DEPTH; i++) {
	    callstack[i] = NULL;
	}
    }

    return err;
}

static void callstack_destroy (void)
{
    int i;

    if (callstack != NULL) {
	for (i=0; i<CALLSTACK_DEPTH; i++) {
	    if (callstack[i] != NULL) {
		free_fncall(callstack[i]);
	    }
	}
    }

    free(callstack);
    callstack = NULL;
}

static int callstack_get_n_calls (void)
{
    int i, n = 0;

    for (i=0; i<CALLSTACK_DEPTH; i++) {
	if (callstack[i] != NULL) n++;
	else break;
    }

    return n;
}

static int push_fncall (fncall *call)
{
    int i, nc;

    if (callstack == NULL && callstack_init()) {
	return E_ALLOC;
    }

    nc = callstack_get_n_calls();
    if (nc == CALLSTACK_DEPTH) {
	strcpy(gretl_errmsg, "Function call stack depth exceeded");
	return 1;
    }

    for (i=nc; i>0; i--) {
	callstack[i] = callstack[i-1];
    }	

    callstack[0] = call;

    return 0;
}

static fncall *pop_fncall (void)
{
    int i, nc;
    fncall *call;

    if (callstack == NULL) return NULL;

    call = callstack[0];
    nc = callstack_get_n_calls();

    for (i=0; i<nc; i++) {
	if (i == nc - 1) {
	    callstack[i] = NULL;
	} else {
	    callstack[i] = callstack[i+1];
	}
    }

    return call;
}

static fncall *current_call (void)
{
    if (callstack == NULL) {
	return NULL;
    } else {
	return callstack[0];
    }
}

/* constructors */

static ufunc *ufunc_new (void)
{
    ufunc *func = malloc(sizeof *func);

    if (func == NULL) return NULL;

    func->name[0] = '\0';
    func_nlines = 0;
    func->lines = NULL;

    return func;
}

static fncall *fncall_new (ufunc *fun, int argc, char **argv)
{
    fncall *call = malloc(sizeof *call);

    if (call == NULL) {
	if (argc > 0) {
	    int i;

	    for (i=0; i<argc; i++) {
		free(argv[i]);
	    }
	    free(argv);
	}
	return NULL;
    }

    call->fun = fun;
    call->argc = argc;
    call->argv = argv;

    return call;
}

static ufunc *add_ufunc (void)
{
    int nf = n_ufuns;
    ufunc **myfuns;

    myfuns = realloc(ufuns, (nf + 1) * sizeof *myfuns);
    if (myfuns == NULL) {
	return NULL;
    }
    ufuns = myfuns;

    ufuns[nf] = ufunc_new();
    if (ufuns[nf] == NULL) {
	return NULL;
    }

    n_ufuns++;

    return ufuns[nf];
}

/* destructors */

static void free_ufunc (ufunc *fun)
{
    int i;

    for (i=0; i<fun->n_lines; i++) {
	free(fun->lines[i]);
    }
    free(fun->lines);

    free(fun);
}

static void free_fncall (fncall *call)
{
    int i;

    for (i=0; i<call->argc; i++) {
	free(call->argv[i]);
    }
    free(call->argv);

    free(call);
}

static ufunc *get_ufunc_by_name (const char *fname)
{
    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(fname, (ufuns[i])->name)) {
	    return ufuns[i];
	}
    }

    return NULL;
}

int gretl_is_user_function (const char *s)
{
    char fname[32];

    sscanf(s, "%31s", fname);
    if (get_ufunc_by_name(fname) != NULL) {
	return 1;
    } else {
	return 0;
    }
}

static int delete_ufunc_by_name (const char *fname)
{
    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(fname, (ufuns[i])->name)) {
	    free_ufunc(ufuns[i]);
	    ufuns[i] = NULL;
	    return 0;
	}
    }

    return 1;
}

static int check_func_name (const char *fname)
{
    int i;

    if (!isalpha((unsigned char) *fname)) {
	strcpy(gretl_errmsg, "function names must start with a letter");
	return 1;
    }

    if (gretl_command_number(fname)) {
	sprintf(gretl_errmsg, "'%s' is the name of a gretl command",
		fname);
	return 1;
    }

    /* or should we overwrite? */

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(fname, (ufuns[i])->name)) {
	    sprintf(gretl_errmsg, "'%s': function is already defined",
		    fname);
	    return 1;
	}
    }

    return 0;
}

static int comma_count (const char *s)
{
    int nc = 0;

    while (*s) {
	if (*s == ',') nc++;
	s++;
    }

    return nc;
}

static char **parse_args (const char *s, int *argc, int *err)
{
    char **argv = NULL;
    int i, na;

    *argc = 0;
    *err = 0;

    /* advance to first arg, if any */

    while (*s) {
	if (*s == ' ') break;
	s++;
    }

    if (*s == '\0') return NULL;

    while (*s) {
	if (*s != ' ') break;
	s++;
    } 

    if (*s == '\0') return NULL;

    /* count comma-separated arguments */
    na = comma_count(++s) + 1;

    argv = malloc(na * sizeof *argv);
    if (argv == NULL) {
	*err = 1;
	return NULL;
    }

    for (i=0; i<na; i++) {
	char *arg;
	int len;

	if (i < na - 1) {
	    len = strcspn(s, ',');
	} else {
	    len = strlen(s);
	}

	arg = gretl_strndup(s, len);
	if (arg == NULL) {
	    na = i;
	    *err = 1;
	    break;
	}
	argv[i] = arg;
    }

    if (*err) {
	for (i=0; i<na; i++) {
	    free(argv[i]);
	}
	free(argv);
	argv = NULL;
    } else {
	*argc = na;
    }

    return argv;
}

int gretl_start_compiling_function (const char *line)
{
    char fname[32];
    ufunc *fun = NULL;
    int err = 0;

    if (!sscanf(line, "function %31s", fname)) {
	return E_PARSE;
    }

    if (check_func_name(fname)) {
	return 1;
    }

    fun = add_ufunc();
    if (fun == NULL) {
	return E_ALLOC;
    }

    strcpy(fun->name, fname);

    set_compiling_on();
    
    return 0;
}

static ufunc *get_latest_ufunc (void)
{
    if (n_ufuns > 0) {
	return ufuns[n_ufuns - 1];
    } else {
	return NULL;
    }
}

int gretl_function_append_line (const char *line)
{
    ufunc *fun = get_latest_ufunc();
    char **lines;
    int nl;

    if (fun == NULL) return 1;

    if (!strncmp(line, "end ", 4)) {
	set_compiling_off();
	return 0;
    }

    nl = fun->n_lines;
    lines = realloc(fun->lines, (nl + 1) * sizeof *lines);
    if (lines == NULL) {
	return E_ALLOC;
    }

    fun->lines = lines;

    fun->lines[nl] = gretl_strdup(line);
    if (fun->lines[nl] == NULL) {
	return E_ALLOC;
    }

    fun->n_lines += 1;

    return 0;
}

static void 
add_args_to_func (ufunc *fun, int argc, char **argv)
{
    if (fun->argc > 0) {
	int i;

	for (i=0; i<fun->argc; i++) {
	    free(fun->argv[i]);
	}
	free(fun->argv);
    }

    fun->argc = argc;
    fun->argv = argv;
}

int gretl_function_start_exec (const char *line)
{
    char **argv;
    char fname[32];
    ufunc *fun;
    fncall *call;
    int argc;
    int err = 0;

    sscanf(line, "%31s", fname);
    fun = get_ufunc_by_name(fname);

    if (fun == NULL) {
	return 1;
    }

    argv = parse_args(line + 1, &argc, &err);

    if (err) {
	return E_ALLOC;
    }

    call = fncall_new(fun, argc, argv);
    if (call == NULL) {
	return E_ALLOC;
    } 

    set_executing_on();
    
    push_fncall(call);

    return 0;
}

static int dollar_term_length (const char *s)
{
    int len = 1;

    s++;
    while (*s) {
	if (isdigit((unsigned char) *s)) len++;
	else break;
	s++;
    }

    return len;
}

static int 
substitute_dollar_terms (char *s, int argc, const char **argv)
{
    char *p;
    int pos, err = 0;

    while ((p = strstr(s, "$")) != NULL) {
	int dlen;

	/* got a positional parameter? */
	if (!sscanf(p, "$%d", pos)) {
	    continue;
	}
	dlen = dollar_term_length(p);
	if (pos >= argc) {
	    /* blank the field */
	    memmove(p, p + dlen, strlen(p + dlen) + 1);
	} else {
	    /* make the substitution */
	    char *q = malloc(strlen(p));

	    if (q == NULL) {
		err = 1;
		break;
	    }
	    strcpy(q, p + dlen);
	    strcpy(p, argv[pos - 1]);
	    strcpy(p + strlen(argv[pos - 1]), q);
	    free(q);	
	}
    }

    return err;
}

int gretl_function_substitute_args (char *s)
{
    fncall *call = current_call();
    int err;

    if (call == NULL || call->fun == NULL) {
	strcpy(gretl_errmsg, "Couldn't find function");
	return 1;
    }

    err = substitute_dollar_terms(s, call->fun->argc,
				  call->fun->argv);

    return err;
}
