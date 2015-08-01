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

#include "libgretl.h"
#include "matrix_extra.h"
#include "version.h"

/* run the vif regression for regressor k */

static double get_vif (MODEL *mod, const int *xlist, 
		       int *vlist, int k,
		       DATASET *dset,
		       int *err)
{
    double vk = NADBL;
    int i, j;

    vlist[1] = xlist[k]; /* dep. var. is regressor k */
    /* position 2 in vlist holds 0 = const */
    j = 3;
    for (i=1; i<=xlist[0]; i++) {
	if (i != k) {
	    vlist[j++] = xlist[i];
	}
    }

    *mod = lsq(vlist, dset, OLS, OPT_A); 
    *err = mod->errcode;

    if (!*err && !xna(mod->rsq) && mod->rsq != 1.0) {
	vk = 1.0 / (1.0 - mod->rsq);
    }

    clear_model(mod);

    return vk;
}

/* run regressions of each x_i on the other x_j's */

static double *model_vif_vector (MODEL *pmod, const int *xlist,
				 DATASET *dset, int *err)
{
    MODEL tmpmod;
    double *vif = NULL;
    int *vlist = NULL;
    int nvif = xlist[0];
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int i;

    if (nvif <= 1) {
	gretl_errmsg_set(_("The statistic you requested is not meaningful "
			   "for this model"));
	return NULL;
    }

    vif = malloc(nvif * sizeof *vif);
    if (vif == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* vlist is the list for the vif regressions:
       allow space for the constant */
    vlist = gretl_list_new(nvif + 1);
    if (vlist == NULL) {
	*err = E_ALLOC;
	free(vif);
	return NULL;
    }

    /* impose original model sample */
    dset->t1 = pmod->t1;
    dset->t2 = pmod->t2;

    for (i=1; i<=xlist[0] && !*err; i++) {
	vif[i-1] = get_vif(&tmpmod, xlist, vlist, i, dset, err);
    }

    /* reinstate sample */
    dset->t1 = save_t1;
    dset->t2 = save_t2;

    free(vlist);

    if (*err) {
	free(vif);
	vif = NULL;
    }

    return vif;
}

gretl_matrix *bkw_matrix (const gretl_matrix *VCV, int *err)
{
    gretl_matrix *Vi = NULL;
    gretl_matrix *S = NULL;
    gretl_matrix *Q = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *lambda = NULL;
    gretl_matrix *BKW = NULL;
    double x, y;
    int k = VCV->rows;
    int i, j;

    /* invert the covariance matrix */
    Vi = gretl_matrix_copy(VCV);
    if (Vi == NULL) {
	*err = E_ALLOC;
	return NULL;
    }
    
    *err = gretl_invert_symmetric_matrix(Vi);
    if (*err) {
	goto bailout;
    }

    /* allocate workspace */
    S = gretl_identity_matrix_new(k);
    Q = gretl_matrix_alloc(k, k);
    BKW = gretl_matrix_alloc(k, k+2);
    
    if (S == NULL || Q == NULL || BKW == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<k; i++) {
	x = gretl_matrix_get(Vi, i, i);
	gretl_matrix_set(S, i, i, 1/sqrt(x));
    }

    *err = gretl_matrix_qform(S, GRETL_MOD_TRANSPOSE,
			      Vi, Q, GRETL_MOD_NONE);

    if (!*err) {
	*err = gretl_matrix_SVD(Q, NULL, &lambda, &V);
    }
    
    if (*err) {
	goto bailout;
    }

    /* S = (1/lambda) ** ones(k, 1) */
    for (j=0; j<k; j++) {
	x = lambda->val[j];
	for (i=0; i<k; i++) {
	    gretl_matrix_set(S, i, j, 1/x);
	}
    }

    for (i=0; i<k; i++) {
	for (j=0; j<k; j++) {
	    x = gretl_matrix_get(V, j, i);
	    y = gretl_matrix_get(S, i, j);
	    gretl_matrix_set(Q, i, j, x * x * y);
	}
    }    

    for (i=0; i<k; i++) {
	/* compute row sums */
	y = 0.0;
	for (j=0; j<k; j++) {
	    y += gretl_matrix_get(Q, i, j);
	}
	for (j=0; j<k; j++) {
	    x = gretl_matrix_get(Q, i, j);
	    gretl_matrix_set(V, j, i, x/y);
	}	
    }

    y = lambda->val[0];

    /* assemble the matrix to return */
    for (i=0; i<k; i++) {
	x = lambda->val[i];
	gretl_matrix_set(BKW, i, 0, x);
	gretl_matrix_set(BKW, i, 1, sqrt(y / x));
	for (j=0; j<k; j++) {
	    x = gretl_matrix_get(V, i, j);
	    gretl_matrix_set(BKW, i, j+2, x);
	}
    }

 bailout:

    gretl_matrix_free(Vi);
    gretl_matrix_free(S);
    gretl_matrix_free(Q);
    gretl_matrix_free(V);
    gretl_matrix_free(lambda);

    if (*err) {
	gretl_matrix_free(BKW);
	BKW = NULL;
    }

    return BKW;
}

static void BKW_print (gretl_matrix *B, PRN *prn)
{
    const char *strs[] = {
	N_("Belsley-Kuh-Welsch collinearity diagnostics"),
	N_("variance proportions"),
	N_("eigenvalues of X'X, largest to smallest"),
	N_("condition index"),
	N_("note: variance proportions columns sum to 1.0")
    };

    pprintf(prn, "\n%s:\n\n", _(strs[0]));
    bufspace(25, prn);
    pprintf(prn, "--- %s ---\n", _(strs[1]));
    gretl_matrix_print_with_format(B, "%10.3f", 0, 0, prn);
    pprintf(prn, "\n  lambda = %s\n", _(strs[2]));
    pprintf(prn, "  cond   = %s\n", _(strs[3]));
    pprintf(prn, "  %s\n\n", _(strs[4]));
}

static void maybe_truncate_param_name (char *s)
{
    int n = strlen(s);

    if (n > 9) {
	char tmp[VNAMELEN];

	tmp[0] = '\0';
	strncat(tmp, s, 8);
	strcat(tmp, "~");
	strcpy(s, tmp);
    }
}

int print_vifs (MODEL *pmod, DATASET *dset, PRN *prn)
{
    gretl_matrix *V = NULL;
    gretl_matrix *B = NULL;
    double *vif = NULL;
    int *xlist;
    double vj;
    int vi, i, n;
    int maxlen = 0;
    int err = 0;

    /* fetch list of regressors */
    xlist = gretl_model_get_x_list(pmod);
    if (xlist == NULL) {
	return E_DATA;
    }

    /* drop the constant if present in xlist */
    for (i=1; i<=xlist[0]; i++) {
	if (xlist[i] == 0) {
	    gretl_list_delete_at_pos(xlist, i);
	    break;
	}
    }

    vif = model_vif_vector(pmod, xlist, dset, &err);
    if (err) {
	return err;
    }

    pprintf(prn, "\n%s\n", _("Variance Inflation Factors"));
    pprintf(prn, "%s\n", _("Minimum possible value = 1.0"));
    pprintf(prn, "%s\n", _("Values > 10.0 may indicate a collinearity problem"));
    pputc(prn, '\n');

    for (i=1; i<=xlist[0]; i++) {
	vi = xlist[i];
	vj = vif[i-1];
	if (!na(vj)) {
	    n = strlen(dset->varname[vi]);
	    if (n > maxlen) {
		maxlen = n;
	    }
	}
    }

    maxlen = maxlen < 12 ? 12 : maxlen;

    for (i=1; i<=xlist[0]; i++) {
	vi = xlist[i];
	vj = vif[i-1];
	if (!na(vj)) {
	    pprintf(prn, "%*s %8.3f\n", maxlen, dset->varname[vi], vj);
	}
    }
    pputc(prn, '\n');

    pputs(prn, _("VIF(j) = 1/(1 - R(j)^2), where R(j) is the "
		 "multiple correlation coefficient\nbetween "
		 "variable j and the other independent variables"));
    pputc(prn, '\n');

    /* now for some more sophisticated diagnostics */

    V = gretl_vcv_matrix_from_model(pmod, NULL, &err);

    if (!err) {
	B = bkw_matrix(V, &err);
    }

    if (!err) {
	int k = pmod->ncoeff + 2;
	char **S = strings_array_new_with_length(k, VNAMELEN);

	if (S != NULL) {
	    int i;

	    strcpy(S[0], "lambda");
	    strcpy(S[1], "cond");
	    for (i=0; i<pmod->ncoeff; i++) {
		gretl_model_get_param_name(pmod, dset, i, S[i+2]);
		maybe_truncate_param_name(S[i+2]);
	    }
	    gretl_matrix_set_colnames(B, S);
	    BKW_print(B, prn);
	}
    }

    gretl_matrix_free(B);
    gretl_matrix_free(V);

    free(vif);
    free(xlist);

    return 0;
}
