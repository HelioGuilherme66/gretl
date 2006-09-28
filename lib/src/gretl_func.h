/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2004 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef GRETL_FUNC_H
#define GRETL_FUNC_H

typedef enum {
    ARG_NONE = 0,
    ARG_SCALAR,
    ARG_SERIES,
    ARG_LIST,
    ARG_MATRIX,
    ARG_BOOL,
    ARG_INT
} FuncArgTypes;

typedef enum {
    FN_NEEDS_TS = 1, /* function requires time-series data */
    FN_NEEDS_QM,     /* function requires quarterly or monthly data */
    FN_NEEDS_PANEL   /* function requires panel data */
} FuncDataReq;

#define NEEDS_TS    "needs-time-series-data"
#define NEEDS_QM    "needs-qm-data"
#define NEEDS_PANEL "needs-panel-data"

typedef struct ufunc_ ufunc;
typedef struct fnpkg_ fnpkg;
typedef struct fnargs_ fnargs;

struct fnargs_ {
    char *types;
    int nx, nX, nM, nl;
    double *x;
    double **X;
    gretl_matrix **M;
    char **lists;
};

int n_user_functions (void);

const ufunc *get_user_function_by_name (const char *name);

const ufunc *get_user_function_by_index (int idx);

int fn_n_params (const ufunc *fun);

int fn_param_type (const ufunc *fun, int i);

const char *fn_param_name (const ufunc *fun, int i);

double fn_param_default (const ufunc *fun, int i);

double fn_param_minval (const ufunc *fun, int i);

double fn_param_maxval (const ufunc *fun, int i);

int user_func_get_return_types (const ufunc *fun,
				int *n_returns,
				char **return_types);

int user_func_first_return_type (const ufunc *fun);

const char *user_function_name_by_index (int i);

int user_function_index_by_name (const char *name);

int gretl_compiling_function (void);

int gretl_executing_function (void);

int gretl_start_compiling_function (const char *line, PRN *prn);

int gretl_function_append_line (const char *line);

int gretl_is_user_function (const char *line);

int gretl_is_public_user_function (const char *name);

int gretl_get_user_function (const char *line, char **fnname);

int is_user_matrix_function (const char *word);

int gretl_function_start_exec (const char *line, const char *fname,
			       double ***pZ, DATAINFO *pdinfo);

char *gretl_function_get_line (char *line, int len,
			       double ***pZ, DATAINFO **ppdinfo,
			       int *err);

int gretl_function_stack_depth (void);

void gretl_function_stop_on_error (double ***pZ, DATAINFO **ppdinfo, PRN *prn);

int gretl_function_flagged_error (const char *s, PRN *prn);

int gretl_function_set_info (int i, const char *help);

int gretl_function_get_info (int i, const char *key, char const **value);

int gretl_function_print_code (int i, PRN *prn);

void gretl_function_set_private (int i, int priv);

int write_function_package (fnpkg *pkg,
			    const char *fname,
			    int pub, 
			    const int *privlist, 
			    const char *author,
			    const char *version,
			    const char *date,
			    const char *descrip,
			    FuncDataReq dreq,
			    float minver);

int function_package_get_info (const char *fname,
			       fnpkg **ppkg,
			       int *pub, 
			       int **privlist,
			       char **author,
			       char **version,
			       char **date,
			       char **descrip,
			       FuncDataReq *dreq,
			       float *minver);

int check_function_needs (const DATAINFO *pdinfo, FuncDataReq dreq,
			  float minver);

int write_user_function_file (const char *fname);

int user_function_file_is_loaded (const char *fname);

int load_user_function_file (const char *fname);

int get_function_file_info (const char *fname, PRN *prn, char **pname);

int get_function_file_code (const char *fname, PRN *prn, char **pname);

char *get_function_file_header (const char *fname, int *err);

int update_function_from_script (const char *fname, int idx);

int user_function_help (const char *fnname, PRN *prn);

void gretl_functions_cleanup (void);

#endif /* GRETL_FUNC_H */
