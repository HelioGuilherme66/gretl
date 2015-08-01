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

#define ADF_DEBUG 0
#define KPSS_DEBUG 0

#include "libgretl.h" 
#include "transforms.h"

/**
 * SECTION:adf_kpss
 * @short_description: unit root and cointegration tests
 * @title: ADF, KPSS, Engle-Granger 
 * @include: libgretl.h
 *
 * Implementations of the (Augmented) Dickey-Fuller test
 * and the Kwiatkowski, Phillips, Schmidt and Shin test for
 * the presence of a unit root in a time series, along with
 * the Engle-Granger test for cointegration of two or more
 * time series.
 *
 * The Johansen cointegration test is also provided in
 * libgretl; see johansen_test() and johansen_test_simple().
 */

/* codes for deterministic regressors */

typedef enum {
    UR_NO_CONST = 1,
    UR_CONST,
    UR_TREND,
    UR_QUAD_TREND,
    UR_MAX
} DetCode;

/* flags for "special stuff" going on */

typedef enum {
    ADF_EG_TEST   = 1 << 0, /* doing Engle-Granger test */
    ADF_EG_RESIDS = 1 << 1, /* final stage of the above */
    ADF_PANEL     = 1 << 2, /* working on panel data */
    ADF_OLS_FIRST = 1 << 3  /* Perron-Qu, 2007 */
} AdfFlags;

/* automatic lag selection methods */

enum {
    k_AIC = 1,
    k_BIC,
    k_TSTAT
};

typedef struct adf_info_ adf_info;
typedef struct kpss_info_ kpss_info;

struct adf_info_ {
    int v;           /* ID number of series to test (in/out) */
    int order;       /* lag order for ADF (in/out) */
    int kmax;        /* max. order (for testing down) */
    int altv;        /* ID of modified series (detrended) */
    int niv;         /* number of (co-)integrated vars (Engle-Granger) */
    AdfFlags flags;  /* bitflags: see above */
    DetCode det;     /* code for deterministics */
    int nseas;       /* number of seasonal dummies */
    int T;           /* number of obs used in test */
    int df;          /* degrees of freedom, test regression */
    double b0;       /* coefficient on lagged level */
    double tau;      /* test statistic */
    double pval;     /* p-value of test stat */
    int *list;       /* regression list */
    int *slist;      /* list of seasonal dummies, if applicable */
    const char *vname; /* name of series tested */
    gretl_matrix *g; /* GLS coefficients (if applicable) */
};

struct kpss_info_ {
    int T;
    double test;
    double pval;
};

/* replace @y with demeaned or detrended y via GLS */

static int GLS_demean_detrend (double *y, int offset,
			       int T, DetCode det,
			       adf_info *ainfo)
{
    gretl_matrix *yd = NULL;
    gretl_matrix *Xd = NULL;
    gretl_matrix *b = NULL;
    double c;
    int t, xcols;
    int err = 0;

    xcols = (det == UR_TREND)? 2 : 1;

    if (T - xcols <= 0) {
	return E_DF;
    }

    yd = gretl_column_vector_alloc(T);
    Xd = gretl_matrix_alloc(T, xcols);
    b = gretl_column_vector_alloc(xcols);

    if (yd == NULL || Xd == NULL || b == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    c = (det == UR_CONST)? (1.0 - 7.0/T) : (1.0 - 13.5/T);

    /* invalidate pre-offset observations */
    for (t=0; t<offset; t++) {
	y[t] = NADBL;
    }
    y += offset;

    gretl_vector_set(yd, 0, y[0] /* (1 - c) * y[0] ?? */);
    for (t=1; t<T; t++) {
	gretl_vector_set(yd, t, y[t] - c * y[t-1]);
    }

    gretl_matrix_set(Xd, 0, 0, 1);
    if (xcols == 2) {
	gretl_matrix_set(Xd, 0, 1, 1);
    }

    for (t=1; t<T; t++) {
	gretl_matrix_set(Xd, t, 0, 1 - c);
	if (xcols == 2) {
	    gretl_matrix_set(Xd, t, 1, t+1 - t*c);
	}
    }

    err = gretl_matrix_ols(yd, Xd, b, NULL, NULL, NULL);

    if (!err) {
	for (t=0; t<T; t++) {
	    y[t] -= b->val[0];
	    if (xcols == 2) {
		y[t] -= b->val[1] * (t+1);
	    }
	}
    }

 bailout:

    gretl_matrix_free(yd);
    gretl_matrix_free(Xd);

    if (err) {
	gretl_matrix_free(b);
    } else {
	ainfo->g = b;
    }
    
    return err;
}

/* replace @y with demeaned or detrended y via OLS */

static int OLS_demean_detrend (double *y, int offset,
			       int T, DetCode det)
{
    gretl_matrix *yd = NULL;
    gretl_matrix *Xd = NULL;
    gretl_matrix *b = NULL;
    int t, xcols = 1;
    int err = 0;

    if (det == UR_QUAD_TREND) {
	xcols = 3;
    } else if (det == UR_TREND) {
	xcols = 2;
    }

    if (T - xcols <= 0) {
	return E_DF;
    }

    yd = gretl_column_vector_alloc(T);
    Xd = gretl_matrix_alloc(T, xcols);
    b = gretl_column_vector_alloc(xcols);

    if (yd == NULL || Xd == NULL || b == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* invalidate pre-offset observations */
    for (t=0; t<offset; t++) {
	y[t] = NADBL;
    }
    y += offset;    

    for (t=0; t<T; t++) {
	gretl_vector_set(yd, t, y[t]);
	gretl_matrix_set(Xd, t, 0, 1.0);
	if (xcols > 1) {
	    gretl_matrix_set(Xd, t, 1, t+1);
	}
	if (xcols > 2) {
	    gretl_matrix_set(Xd, t, 2, (t+1) * (t+1));
	}
    }

    err = gretl_matrix_ols(yd, Xd, b, NULL, NULL, NULL);

    if (!err) {
	for (t=0; t<T; t++) {
	    y[t] -= b->val[0];
	    if (xcols > 1) {
		y[t] -= b->val[1] * (t+1);
	    }
	    if (xcols > 2) {
		y[t] -= b->val[2] * (t+1) * (t+1);
	    }
	}
    }

 bailout:

    gretl_matrix_free(yd);
    gretl_matrix_free(Xd);
    gretl_matrix_free(b);
    
    return err;
}

static int real_adf_form_list (adf_info *ainfo,
			       DATASET *dset)
{
    int v, save_t1 = dset->t1;
    int k, j, err = 0;

    /* using the original var, or transformed? */
    v = ainfo->altv > 0 ? ainfo->altv : ainfo->v;

    /* temporararily reset sample */
    dset->t1 = 0;

    /* generate the first difference of series @v:
       this will be the LHS variable in the test
    */
    ainfo->list[1] = diffgenr(v, DIFF, dset);
    if (ainfo->list[1] < 0) {
	dset->t1 = save_t1;
	return E_DATA;
    }

    /* generate lag 1 of series @v: the basic RHS series */
    ainfo->list[2] = laggenr(v, 1, dset); 
    if (ainfo->list[2] < 0) {
	dset->t1 = save_t1;
	return E_DATA;
    }

    /* generate lagged differences for augmented test */
    j = 3;
    for (k=1; k<=ainfo->order && !err; k++) {
	int vk = laggenr(ainfo->list[1], k, dset);

	if (vk < 0) {
	    fprintf(stderr, "Error generating lag variable\n");
	    err = E_DATA;
	} else {
	    ainfo->list[j++] = vk;
	} 
    }

    if (!err && ainfo->nseas > 0) {
	/* should we center these? */
	ainfo->slist = seasonals_list(dset, dset->pd, 0, &err);
	if (err) {
	    ainfo->nseas = 0;
	} else {
	    ainfo->nseas = ainfo->slist[0];
	}
    }

    /* restore incoming sample */
    dset->t1 = save_t1;

    return err;
}

/* The "offset" will determine where the data start for the
   initial detrending regression. I suppose we don't want
   to include excessive "pre-sample" data if the user has
   explicitly moved the sample start off zero, but we need
   the data to start early enough to provide @order + 1
   lags, to be in sync with the ADF regressions.
*/

static int adf_y_offset (adf_info *ainfo, int v, DATASET *dset)
{
    int min_offset = dset->t1 - (ainfo->order + 1);
    int t, offset = 0;

    /* copy original data into series v */
    for (t=0; t<=dset->t2; t++) {
	dset->Z[v][t] = dset->Z[ainfo->v][t];
	if (na(dset->Z[v][t])) {
	    offset = t+1;
	}
    }

    if (offset < min_offset) {
	offset = min_offset;
    }

#if ADF_DEBUG
    fprintf(stderr, "adf_y_offset: dset->t1=%d, order=%d, offset=%d\n",
	    dset->t1, ainfo->order, offset);
#endif
    
    return offset;
}

/* Generate the various differences and lags required for
   the ADF test, detrending first if required.
*/

static int adf_prepare_vars (adf_info *ainfo, DATASET *dset,
			     gretlopt opt)
{
    int err = 0;

    if (ainfo->v == 0) {
	return E_DATA;
    }

    if (ainfo->list == NULL) {
	/* the max number of terms (quadratic trend case) */
	int len = 5 + ainfo->nseas + ainfo->order;
	
	ainfo->list = gretl_list_new(len);
	if (ainfo->list == NULL) {
	    return E_ALLOC;
	}
    }

    if (ainfo->flags & ADF_OLS_FIRST) {
	/* OLS adjustment is wanted (first pass) */
	DetCode det = UR_CONST;
	int v = dset->v;

	if (opt & OPT_R) {
	    det = UR_QUAD_TREND;
	} else if (opt & OPT_T) {
	    det = UR_TREND;
	}

	err = dataset_add_series(dset, 1);
	
	if (!err) {
	    int offset = adf_y_offset(ainfo, v, dset);
	    int T = dset->t2 - offset + 1;

	    err = OLS_demean_detrend(dset->Z[v], offset, T, det);
	}
	
	if (!err) {
	    /* replace with OLS-detrended version */
	    strcpy(dset->varname[v], "ydetols");
	    ainfo->altv = v;
	}	
    } else if (opt & OPT_G) {
	/* GLS adjustment is wanted */
	DetCode det = (opt & OPT_T)? UR_TREND : UR_CONST;
	int v;

	if (ainfo->altv > 0) {
	    /* if we already did detrending, re-use the
	       storage we added earlier */
	    v = ainfo->altv;
	} else {
	    v = dset->v;
	    err = dataset_add_series(dset, 1);
	}
	
	if (!err) {
	    int offset = adf_y_offset(ainfo, v, dset);
	    int T = dset->t2 - offset + 1;

	    err = GLS_demean_detrend(dset->Z[v], offset, T,
				     det, ainfo);
	}
	
	if (!err) {
	    /* replace with GLS-detrended version */
	    strcpy(dset->varname[v], "ydetrend");
	    ainfo->altv = v;
	}
    }

    if (!err) {
	err = real_adf_form_list(ainfo, dset);
    }

#if ADF_DEBUG
    printlist(ainfo->list, "adf initial list");
#endif

    return err;
}

#if 0 /* based on a much larger replication of ERS:
	 activate this after putting an account of
	 the simulation in place */

static void get_df_gls_ct_cval (int T, double *c)
{
    static double b[4][3] = {
	/* b0       b(1/T)    b(1/T^2) */
	{ -2.56073, -17.5434,  68.4750 }, /* 10% */
	{ -2.84864, -17.6702,  36.9221 }, /* 5% */
	{ -3.10420, -18.2513,  9.99274 }, /* 2.5% */
	{ -3.40846, -19.4237, -23.9869 }  /* 1% */
    };
    double T2 = T * T;
    int i;

    for (i=0; i<4; i++) {
	c[i] = b[i][0] + b[i][1] / T + b[i][2] / T2;
    }
}

#else

static void get_df_gls_ct_cval (int T, double *c)
{
    /* Elliott, Rothenberg and Stock (1996), Table 1 */
    static double df_gls_ct_cvals[4][4] = {
	/* 10%     5%    2.5%    1% */
	{ -2.89, -3.19, -3.46, -3.77 }, /* T = 50  */
	{ -2.74, -3.03, -3.29, -3.58 }, /* T = 100 */
	{ -2.64, -2.93, -3.18, -3.46 }, /* T = 200 */
	{ -2.57, -2.89, -3.15, -3.48 }  /* \infty  */
    };
    int j, i = 3;

    if (T <= 50){
	i = 0;
    } else if (T <= 100){
	i = 1;
    } else if (T <= 200){
	i = 2;
    }

    for (j=0; j<4; j++) {
	c[j] = df_gls_ct_cvals[i][j];
    }
}

#endif

/* display an F-test for the joint significance of the lagged
   \delta y terms in ADF test */

static void show_lags_test (MODEL *pmod, int order, PRN *prn)
{
    int *llist = gretl_list_new(order);
    double F;
    int i;

    if (llist != NULL) {
	for (i=0; i<order; i++) {
	    /* lagged differences */
	    llist[i+1] = pmod->list[pmod->ifc + 3 + i];
	}

	F = wald_omit_F(llist, pmod);

	if (!na(F)) {
	    pprintf(prn, "  %s: F(%d, %d) = %.3f [%.4f]\n",
		    _("lagged differences"), order, pmod->dfd, F, 
		    snedecor_cdf_comp(order, pmod->dfd, F));
	}

	free(llist);
    }
}

static const char *test_down_string (int i, gretlopt opt)
{
    if (opt & OPT_U) {
	/* perron-qu */
	if (i == k_BIC) {
	    return _("modified BIC, Perron-Qu");
	} else {
	    return _("modified AIC, Perron-Qu");
	}
    } else {
	int gls = (opt & OPT_G);
	
	if (i == k_BIC) {
	    return gls ? _("modified BIC") : _("BIC");
	} else if (i == k_TSTAT) {
	    return _("t-statistic");
	} else {
	    /* the default */
	    return gls ? _("modified AIC") : _("AIC");
	}
    }
}

static void DF_header (const char *s, int p, int pmax,
		       int test_down, gretlopt opt, 
		       PRN *prn)
{
    pputc(prn, '\n');

    if (p <= 0 && pmax == 0) {
	if (opt & OPT_G) {
	    pprintf(prn, _("Dickey-Fuller (GLS) test for %s\n"), s);
	} else {
	    pprintf(prn, _("Dickey-Fuller test for %s\n"), s);
	}
    } else {
	if (opt & OPT_G) {
	    pprintf(prn, _("Augmented Dickey-Fuller (GLS) test for %s\n"), s);
	} else {
	    pprintf(prn, _("Augmented Dickey-Fuller test for %s\n"), s);
	}
	if (p == 1) {
	    pprintf(prn, _("including one lag of (1-L)%s"), s);
	} else {
	    pprintf(prn, _("including %d lags of (1-L)%s"), p, s);
	}
	if (pmax >= p) {
	    const char *critstr = test_down_string(test_down, opt);

	    pputc(prn, '\n');
	    pprintf(prn, _("(max was %d, criterion %s)"), 
		    pmax, critstr);
	}
	pputc(prn, '\n');
    }
}

static const char *DF_model_string (int i)
{
    const char *models[] = {
	"(1-L)y = (a-1)*y(-1) + e",
	"(1-L)y = b0 + (a-1)*y(-1) + e",
	"(1-L)y = b0 + b1*t + (a-1)*y(-1) + e",
	"(1-L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + e"
    };

    if (i >= 0 && i < 4) {
	return models[i];
    } else {
	return "";
    }
}

static const char *ADF_model_string (int i)
{
    const char *models[] = {
	"(1-L)y = (a-1)*y(-1) + ... + e",
	"(1-L)y = b0 + (a-1)*y(-1) + ... + e",
	"(1-L)y = b0 + b1*t + (a-1)*y(-1) + ... + e",
	"(1-L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + ... + e"
    };

    if (i >= 0 && i < 4) {
	return models[i];
    } else {
	return "";
    }
}

static const char *DF_test_string (int i)
{
    const char *tests[] = {
	N_("test without constant"),
	N_("test with constant"),
	N_("with constant and trend"),
	N_("with constant and quadratic trend")
    };

    if (i >= 0 && i < 4) {
	return tests[i];
    } else {
	return "";
    }
}

static void print_adf_results (adf_info *ainfo, MODEL *dfmod,
			       int *blurb_done, gretlopt opt,
			       int test_down, PRN *prn)
{
    const char *urcstrs[] = {
	"nc", "c", "ct", "ctt"
    };
    char pvstr[48];
    char taustr[16];
    int i;

    if (prn == NULL) return;

    /* convert deterministics code to 0-base */
    i = ainfo->det - 1;

    if (na(ainfo->pval)) {
	sprintf(pvstr, "%s %s", _("p-value"), _("unknown"));
    } else {
	int asy = (ainfo->order > 0 || (opt & OPT_G));

	sprintf(pvstr, "%s %.4g", 
		(asy)? _("asymptotic p-value") : _("p-value"), 
		ainfo->pval);
    } 

    if (*blurb_done == 0) {
	DF_header(ainfo->vname, ainfo->order, ainfo->kmax,
		  test_down, opt, prn);
	pprintf(prn, _("sample size %d\n"), ainfo->T);
	if (ainfo->flags & ADF_PANEL) {
	    pputc(prn, '\n');
	} else {
	    pputs(prn, _("unit-root null hypothesis: a = 1"));
	    pputs(prn, "\n\n");
	}
	*blurb_done = 1;
    }

    if (ainfo->flags & ADF_EG_RESIDS) {
	/* last step of Engle-Granger test */
	pprintf(prn, "  %s: %s\n", _("model"), 
		(ainfo->order > 0)? ADF_model_string(0) :
		DF_model_string(0));
    } else {
	pprintf(prn, "  %s ", _(DF_test_string(i)));
	if (ainfo->nseas > 0 && i > 0) {
	    pputs(prn, _("plus seasonal dummies"));
	}
	pputc(prn, '\n');
	pprintf(prn, "  %s: %s\n", _("model"), 
		(ainfo->order > 0)? ADF_model_string(i) :
		DF_model_string(i));
    }

    if (opt & OPT_G) {
	strcpy(taustr, "tau");
    } else {
	sprintf(taustr, "tau_%s(%d)", urcstrs[i], ainfo->niv);
    }

    pprintf(prn, "  %s: %g\n"
	    "  %s: %s = %g\n",
	    _("estimated value of (a - 1)"), ainfo->b0,
	    _("test statistic"), taustr, ainfo->tau);

    if ((opt & OPT_G) && i+1 == UR_TREND) {
	double c[4];

	get_df_gls_ct_cval(ainfo->T, c);
	pprintf(prn, "\n  %*s    ", TRANSLATED_WIDTH(_("Critical values")), " ");
	pprintf(prn, "%g%%     %g%%     %g%%     %g%%\n", 10.0, 5.0, 2.5, 1.0);
	pprintf(prn, "  %s: %.2f   %.2f   %.2f   %.2f\n", 
		_("Critical values"), c[0], c[1], c[2], c[3]);
    } else {
	pprintf(prn, "  %s\n", pvstr);
    }

    if (!na(dfmod->rho)) {
	pprintf(prn, "  %s: %.3f\n", _("1st-order autocorrelation coeff. for e"), 
		dfmod->rho);
    }

    if (ainfo->order > 1) {
	show_lags_test(dfmod, ainfo->order, prn);
    } 
}

/* test the lag order down using the t-statistic criterion */

static int t_adjust_order (adf_info *ainfo, DATASET *dset,
			   int *err, PRN *prn)
{
    gretlopt kmod_opt = (OPT_A | OPT_Z);
    MODEL kmod;
    int kmax = ainfo->kmax;
    double tstat, pval;
    int k, pos;

    for (k=kmax; k>0; k--) {
	kmod = lsq(ainfo->list, dset, OLS, kmod_opt);
	if (!kmod.errcode && kmod.dfd == 0) {
	    kmod.errcode = E_DF;
	}
	if (kmod.errcode) {
	    fprintf(stderr, "t_adjust_order: k = %d, err = %d\n", k,
		    kmod.errcode);
	    *err = kmod.errcode;
	    clear_model(&kmod);
	    k = -1;
	    break;
	}
#if ADF_DEBUG
	printmodel(&kmod, dset, OPT_NONE, prn);
#endif
	pos = k + kmod.ifc;
	tstat = kmod.coeff[pos] / kmod.sderr[pos];
	clear_model(&kmod);
	pval = normal_pvalue_2(tstat);

	if (pval > 0.10) {
#if ADF_DEBUG
	    pprintf(prn, "\nt_adjust_order: lagged difference not "
		    "significant at order %d (t = %g)\n\n", k, tstat);
#endif
	    gretl_list_delete_at_pos(ainfo->list, k + 2);
	} else {
#if ADF_DEBUG
	    pprintf(prn, "\nt_adjust_order: lagged difference is "
		    "significant at order %d (t = %g)\n\n", k, tstat);
#endif
	    break;
	}
    }

    return k;
}

/* compute modified information criterion */

static double get_MIC (MODEL *pmod, int k, double sum_ylag2,
		       int kmethod, const DATASET *dset)
{
    double g, CT, ttk, s2k = 0;
    int t, T = pmod->nobs;

    g = pmod->coeff[pmod->ifc];

    for (t=pmod->t1; t<=pmod->t2; t++) {
	s2k += pmod->uhat[t] * pmod->uhat[t];
    }
    
    s2k /= pmod->nobs;
    ttk = g * g * sum_ylag2 / s2k;
    CT = kmethod == k_BIC ? log(T) : 2.0;

    return log(s2k) + CT * (ttk + k)/T;
}

/* component calculation for modified IC methods */

static double get_sum_y2 (adf_info *ainfo, MODEL *pmod,
			  const DATASET *dset)
{
    const double *ylag = dset->Z[ainfo->list[2]];
    double sumy2 = 0;
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	sumy2 += ylag[t] * ylag[t];
    }

    return sumy2;
}

/* Using modified information criterion, as per Ng and Perron,
   "Lag Length Selection and the Construction of Unit Root Tests 
   with Good Size and Power", Econometrica 69/6, Nov 2001, pp. 
   1519-1554, for the GLS case -- otherwise plain IC (as of
   2015-03-31).
*/

static int ic_adjust_order (adf_info *ainfo, int kmethod,
			    DATASET *dset, gretlopt opt,
			    int test_num, int *err,
			    PRN *prn)
{
    MODEL kmod;
    gretlopt kmod_opt = (OPT_A | OPT_Z);
    double IC, ICmin = 0;
    double sum_ylag2 = 0;
    int kmax = ainfo->kmax;
    int k, kopt = kmax;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int use_MIC = 0;
    int *tmplist;

    tmplist = gretl_list_copy(ainfo->list);
    if (tmplist == NULL) {
	*err = E_ALLOC;
	return -1;
    }

    if (opt & OPT_G) {
	/* ADF-GLS */
	use_MIC = 1;
    }

    for (k=kmax; k>=0; k--) {
	kmod = lsq(tmplist, dset, OLS, kmod_opt);
	if (!kmod.errcode && kmod.dfd == 0) {
	    kmod.errcode = E_DF;
	}
	if (kmod.errcode) {
	    fprintf(stderr, "ic_adjust_order: k = %d, err = %d\n", k,
		    kmod.errcode);
	    *err = kmod.errcode;
	    clear_model(&kmod);
	    kopt = -1;
	    break;
	}
	if (use_MIC) {
	    if (k == kmax) {
		/* this need only be done once */
		sum_ylag2 = get_sum_y2(ainfo, &kmod, dset);
	    }
	    IC = get_MIC(&kmod, k, sum_ylag2, kmethod, dset);
	} else if (kmethod == k_BIC) {
	    IC = kmod.criterion[C_BIC];
	} else {
	    IC = kmod.criterion[C_AIC];
	}
	if (k == kmax) {
	    /* ensure a uniform sample */
	    dset->t1 = kmod.t1;
	    dset->t2 = kmod.t2;
	    ICmin = IC;
	} else if (IC < ICmin) {
	    ICmin = IC;
	    kopt = k;
	}
#if ADF_DEBUG
	printmodel(&kmod, dset, OPT_NONE, prn);
#endif
	if (opt & OPT_V) {
	    const char *tag;

	    if (use_MIC) {
		tag = (kmethod == k_BIC) ? "MBIC" : "MAIC";
	    } else {
		tag = (kmethod == k_BIC) ? "BIC" : "AIC";
	    }

	    if (k == kmax && test_num == 1) {
		pputc(prn, '\n');
	    }
	    pprintf(prn, "  k = %2d: %s = %#g\n", k, tag, IC);
	}	    
	clear_model(&kmod);
	gretl_list_delete_at_pos(tmplist, k + 2);
    }

    if ((opt & OPT_V) && test_num > 1) {
	pputc(prn, '\n');
    }

    free(tmplist);

    if (kopt >= 0) {
	/* now trim the "real" list to @kopt lags */
	for (k=kmax; k>kopt; k--) {
	    gretl_list_delete_at_pos(ainfo->list, k + 2);
	}
    }

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return kopt;
}

/* targ must be big enough to accept all of src! */

static void copy_list_values (int *targ, const int *src)
{
    int i;

    for (i=0; i<=src[0]; i++) {
	targ[i] = src[i];
    }
}

/**
 * get_urc_pvalue:
 * @tau: test statistic.
 * @n: sample size (or 0 for asymptotic result).
 * @niv: number of potentially cointegrated variables
 * (1 for simple unit-root test).
 * @itv: code: 1, 2, 3, 4 for nc, c, ct, ctt models
 * respectively.
 * @opt: give OPT_G if GLS adjustment was applied in
 * the test from which @tau was obtained.
 *
 * Retrieves the p-value for @tau from the Dickey–Fuller 
 * unit-root test or the Engle–Granger cointegration 
 * test, as per James MacKinnon (1996).
 *
 * Returns: p-value, or %NADBL on failure.
 */

double get_urc_pvalue (double tau, int n, int niv, int itv,
		       gretlopt opt)
{
    char datapath[FILENAME_MAX];
    double (*mackinnon_pvalue)(double, int, int, int, char *);
    double pval = NADBL;
    static int nodata;
    
    if (nodata) {
	return pval;
    }

    mackinnon_pvalue = get_plugin_function("mackinnon_pvalue");
    if (mackinnon_pvalue == NULL) {
	nodata = 1;
        return pval;
    }

    strcpy(datapath, gretl_lib_path());
#ifdef WIN32
    append_dir(datapath, "plugins");
#endif

    if ((opt & OPT_G) && itv == UR_CONST) {
	itv = UR_NO_CONST;
    }

    pval = (*mackinnon_pvalue)(tau, n, niv, itv, datapath);

#if ADF_DEBUG
    fprintf(stderr, "getting pval: tau=%g, n=%d, niv=%d, itv=%d: pval=%g\n",
	    tau, n, niv, itv, pval);
#endif

    if (*datapath == '\0') {
	nodata = 1;
    } 

    return pval;
}

#define test_opt_not_set(o) (!(o & OPT_N) && !(o & OPT_C) && \
                             !(o & OPT_T) && !(o & OPT_R))

static int test_wanted (DetCode det, gretlopt opt)
{
    int ret = 0;

    switch (det) {
    case UR_NO_CONST:
	ret = (opt & OPT_N);
	break;
    case UR_CONST:
	ret = (opt & OPT_C);
	break;
    case UR_TREND:
	ret = (opt & OPT_T);
	break;
    case UR_QUAD_TREND:
	ret = (opt & OPT_R);
	break;
    default:
	break;
    }

    return ret;
}

static DetCode engle_granger_itv (gretlopt opt)
{
    DetCode itv = UR_CONST;

    if (opt & OPT_N) {
	itv = UR_NO_CONST;
    } else if (opt & OPT_T) {
	itv = UR_TREND;
    } else if (opt & OPT_R) {
	itv = UR_QUAD_TREND;
    }

    return itv;
}

static int gettrend (DATASET *dset, int square)
{
    int idx, t, v = dset->v;
    double x;

    idx = series_index(dset, (square)? "timesq" : "time");

    if (idx < v) {
	return idx;
    }

    if (dataset_add_series(dset, 1)) {
	return 0; /* error: valid value cannot == 0 */
    }

    for (t=0; t<dset->n; t++) {
	x = (double) t + 1; 
	dset->Z[v][t] = (square)? x * x : x;
    }

    if (square) {
	strcpy(dset->varname[v], "timesq");
	series_set_label(dset, v, _("squared time trend variable"));
    } else {
	strcpy(dset->varname[v], "time");
	series_set_label(dset, v, _("time trend variable"));
    }
	    
    return idx;
}

static void print_df_model (adf_info *ainfo, MODEL *pmod,
			    int dfnum, DATASET *dset,
			    PRN *prn)
{
    pmod->aux = (ainfo->order > 0)? AUX_ADF : AUX_DF;
    
    if (!na(ainfo->pval)) {
	gretl_model_set_int(pmod, "dfnum", dfnum);
	gretl_model_set_double(pmod, "dfpval", ainfo->pval);
    }
    
    if (ainfo->flags & ADF_EG_RESIDS) {
	gretl_model_set_int(pmod, "eg-resids", 1);
    }

    if (ainfo->g != NULL) {
	gretl_model_set_int(pmod, "dfgls", 1);
    }
    
    printmodel(pmod, dset, OPT_NONE, prn);

    if (ainfo->g != NULL) {
	pputs(prn, "  ");
	if (ainfo->g->rows == 2) {
	    pprintf(prn, _("GLS detrending: b0 = %g, b1 = %g\n"),
		    ainfo->g->val[0], ainfo->g->val[1]);
	} else {
	    pprintf(prn, _("GLS estimate of b0: %g\n"),
		    ainfo->g->val[0]);
	}
	pputc(prn, '\n');
    }    
}

static int set_deterministic_terms (adf_info *ainfo,
				    DATASET *dset)
{
    int i, j;

    /* Note that list[1] and list[2], plus the @order 
       lagged differences, are in common for all 
       specifications 
    */

    ainfo->list[0] = 1 + ainfo->order + ainfo->det;

    if (ainfo->det >= UR_TREND) {
	i = 3 + ainfo->order;
	ainfo->list[i] = gettrend(dset, 0);
	if (ainfo->list[i] == 0) {
	    return E_ALLOC;
	}
    }

    if (ainfo->det == UR_QUAD_TREND) {
	i = 4 + ainfo->order;
	ainfo->list[i] = gettrend(dset, 1);
	if (ainfo->list[i] == 0) {
	    return E_ALLOC;
	}
    }

    if (ainfo->det != UR_NO_CONST) {
	i = ainfo->list[0];
	ainfo->list[0] += ainfo->nseas;
	/* stick constant on end of list */
	ainfo->list[ainfo->list[0]] = 0;
	/* preceded by seasonal dummies if wanted */
	for (j=0; j<ainfo->nseas; j++) {
	    ainfo->list[i++] = ainfo->slist[j+1];
	}	    
    } 

    return 0;
}

/* When we're done with testing down in the DF-GLS context
   we may need to regenerate the detrended data: this
   applies if (a) we used Perron-Qu OLS detrending in the
   test-down phase or (b) the sample start is currently
   set after the start of the dataset.
*/

static int reset_detrended_data (adf_info *ainfo,
				 DATASET *dset,
				 gretlopt opt)
{
    if (ainfo->flags & ADF_OLS_FIRST) {
	/* testing down using Perron-Qu */
	return 1;
    } else if ((opt & OPT_G) && dset->t1 > 0) {
	/* GLS + testing down + t1 > 0 */
	return 1;
    } else {
	return 0;
    }
}

static int handle_test_down_option (adf_info *ainfo,
				    gretlopt opt,
				    int *err)
{
    int kmethod = 0;
    const char *s;

    if (ainfo->flags & (ADF_EG_TEST | ADF_EG_RESIDS)) {
	s = get_optval_string(COINT, OPT_E);
    } else {
	s = get_optval_string(ADF, OPT_E);
    }

    if (s == NULL || *s == '\0') {
	/* the default */
	kmethod = k_AIC;
    } else if (!strcmp(s, "MAIC") || !strcmp(s, "AIC")) {
	kmethod = k_AIC;
    } else if (!strcmp(s, "MBIC") || !strcmp(s, "BIC")) {
	kmethod = k_BIC;
    } else if (!strcmp(s, "tstat")) {
	kmethod = k_TSTAT;
    } else {
	gretl_errmsg_set(_("Invalid option"));
	*err = E_DATA;
    }

    if (!*err) {
	/* take the given order to be the max */
	ainfo->kmax = ainfo->order;
	if (opt & OPT_U) {
	    /* --perron-qu */
	    if (kmethod != k_AIC && kmethod != k_BIC) {
		*err = E_BADOPT;
	    } else {
		ainfo->flags |= ADF_OLS_FIRST;
	    }
	}
    }

    return kmethod;
}

static int check_adf_options (gretlopt opt)
{
    int err = 0;

    if (opt & OPT_G) {
	/* we can only have one of the basic deterministics options */
	err = incompatible_options(opt, OPT_N | OPT_C | OPT_T | OPT_R);
	/* options incompatible with --gls: no-const, seasonals,
	   and quadratic trend
	*/
	if (opt & (OPT_N | OPT_D | OPT_R)) {
	    err = E_BADOPT;
	}
    } else if (opt & OPT_U) {
	/* option dependent on --gls: Perron-Qu modified AIC/BIC */
	err = E_BADOPT;
    }

    return err;
}

static int real_adf_test (adf_info *ainfo, DATASET *dset,
			  gretlopt opt, PRN *prn)
{
    MODEL dfmod;
    gretlopt eg_opt = OPT_NONE;
    gretlopt df_mod_opt = (OPT_A | OPT_Z);
    int *biglist = NULL;
    int orig_nvars = dset->v;
    int blurb_done = 0;
    int test_down = 0;
    int test_num = 0;
    int i, err;

    /* (most of) this may have been done already
       but it won't hurt to check here */
    err = check_adf_options(opt);
    if (err) {
	return err;
    }

    /* safety-first initializations */
    ainfo->nseas = ainfo->kmax = ainfo->altv = 0;
    ainfo->vname = dset->varname[ainfo->v];
    ainfo->list = ainfo->slist = NULL;
    ainfo->det = 0;

#if ADF_DEBUG
    fprintf(stderr, "real_adf_test: got order = %d\n", ainfo->order);
#endif

    if (gretl_isconst(dset->t1, dset->t2, dset->Z[ainfo->v])) {
	gretl_errmsg_sprintf(_("%s is a constant"), ainfo->vname);
	return E_DATA;
    }    

    if (opt & OPT_E) {
	/* --test-down[=...] */
	test_down = handle_test_down_option(ainfo, opt, &err);
	if (err) {
	    return err;
	}
    }

    if (ainfo->order < 0) {
	test_down = k_AIC;
	ainfo->order = ainfo->kmax = -ainfo->order;
    }

#if ADF_DEBUG
    fprintf(stderr, "real_adf_test: order = %d, test_down = %d\n",
	    ainfo->order, test_down);
#endif

    if (ainfo->flags & ADF_EG_RESIDS) {
	/* Final step of Engle-Granger test: the (A)DF test
	   regression will contain no deterministic terms, 
	   but the selection of the p-value is based on the
	   deterministic terms in the cointegrating
	   regression, represented by @eg_opt.
	*/
	int verbose = (opt & OPT_V);
	int silent = (opt & OPT_I);

	eg_opt = opt;
	opt = OPT_N;
	if (silent) {
	    opt |= OPT_I;
	} else if (verbose) {
	    opt |= OPT_V;
	}
    }

    if (opt & OPT_F) {
	/* difference the target series before testing */
	int t1 = dset->t1;

	dset->t1 = 0;
	ainfo->v = diffgenr(ainfo->v, DIFF, dset);
	dset->t1 = t1;
	if (ainfo->v < 0) {
	    return E_DATA;
	}
	ainfo->vname = dset->varname[ainfo->v];
    }

    if ((opt & OPT_D) && dset->pd > 1) {
	/* arrange to add seasonal dummies */
	ainfo->nseas = dset->pd - 1;
    }

    if (test_opt_not_set(opt)) {
	/* default model(s) */
	if (opt & OPT_G) {
	    opt |= OPT_C;
	} else {
	    opt |= (OPT_C | OPT_T);
	}
    }

    err = adf_prepare_vars(ainfo, dset, opt);
    if (err) {
	return err;
    }

    if (test_down) {
	ainfo->list[0] = ainfo->order + 5;
	biglist = gretl_list_copy(ainfo->list);
	if (biglist == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    gretl_model_init(&dfmod, dset);

    /* Now loop across the wanted deterministics cases:
       in many instances we'll actually be doing only one
       case.
    */

    for (i=UR_NO_CONST; i<UR_MAX; i++) {
	int b0pos = (i > UR_NO_CONST);

	ainfo->det = i;

	if (!test_wanted(ainfo->det, opt)) {
	    continue;
	}

	if (test_down) {
	    /* re-establish max order before testing down */
	    ainfo->order = ainfo->kmax;
	    copy_list_values(ainfo->list, biglist);
	}

	if (opt & OPT_G) {
	    /* DF-GLS: skip deterministics */
	    ainfo->list[0] = ainfo->order + 2;
	    b0pos = 0;
	} else {
	    err = set_deterministic_terms(ainfo, dset);
	    if (err) {
		goto bailout;
	    }
	}

	test_num++;

	if (test_down) {
	    /* determine the optimal lag order */
	    if (test_down == k_TSTAT) {
		ainfo->order = t_adjust_order(ainfo, dset, &err, prn);
	    } else {
		ainfo->order = ic_adjust_order(ainfo, test_down,
					       dset, opt, test_num,
					       &err, prn);
	    }	    
	    if (err) {
		clear_model(&dfmod);
		goto bailout;
	    } else if (reset_detrended_data(ainfo, dset, opt)) {
		/* swap out the detrended data */
		ainfo->flags &= ~ADF_OLS_FIRST;
		err = adf_prepare_vars(ainfo, dset, opt);
	    }	    
	}

#if ADF_DEBUG
	printlist(ainfo->list, "final ADF regression list");
#endif

	/* run the actual test regression */
	dfmod = lsq(ainfo->list, dset, OLS, df_mod_opt);

	if (!dfmod.errcode && dfmod.dfd == 0) {
	    /* we can't tolerate an exact fit here */
	    dfmod.errcode = E_DF;
	}
	if (dfmod.errcode) {
	    fprintf(stderr, "adf_test: dfmod.errcode = %d\n", 
		    dfmod.errcode);
	    err = dfmod.errcode;
	    clear_model(&dfmod);
	    goto bailout;
	}

	/* transcribe info from test regression */
	ainfo->b0 = dfmod.coeff[b0pos];
	ainfo->tau = ainfo->b0 / dfmod.sderr[b0pos];
	ainfo->T = dfmod.nobs;
	ainfo->df = dfmod.dfd;

	if (ainfo->flags & ADF_EG_RESIDS) {
	    ainfo->det = engle_granger_itv(eg_opt);
	} 

	if (getenv("DFGLS_NO_PVALUE")) {
	    /* to speed up monte carlo stuff */
	    ainfo->pval = NADBL;
	} else if ((opt & OPT_G) && ainfo->det == UR_TREND) {
	    /* DF-GLS with trend: MacKinnon p-values won't work */
	    ainfo->pval = NADBL;
	} else {
	    /* Use asymp. p-value in augmented case; also the
	       finite-sample MacKinnon p-values are not correct
	       in the GLS case.
	    */
	    int asymp = (ainfo->order > 0 || (opt & OPT_G));

	    ainfo->pval = get_urc_pvalue(ainfo->tau, asymp ? 0 : dfmod.nobs, 
					 ainfo->niv, ainfo->det, opt);
	}

	if (!(opt & (OPT_Q | OPT_I)) && !(ainfo->flags & ADF_PANEL)) {
	    print_adf_results(ainfo, &dfmod, &blurb_done,
			      opt, test_down, prn);
	}

	if ((opt & OPT_V) && !(ainfo->flags & ADF_PANEL)) {
	    /* verbose */
	    print_df_model(ainfo, &dfmod, b0pos, dset, prn);
	} else if (!(opt & OPT_Q) && !(ainfo->flags & (ADF_EG_RESIDS | ADF_PANEL))) {
	    pputc(prn, '\n');
	}

	clear_model(&dfmod);
    }

    if (!err) {
	if (!(ainfo->flags & (ADF_EG_TEST | ADF_PANEL)) ||
	    (ainfo->flags & ADF_EG_RESIDS)) {
	    record_test_result(ainfo->tau, ainfo->pval, "Dickey-Fuller");
	}
    }

 bailout:

    free(ainfo->list);
    ainfo->list = NULL;

    free(ainfo->slist);
    ainfo->slist = NULL;

    gretl_matrix_free(ainfo->g);
    ainfo->g = NULL;
    
    free(biglist);

    dataset_drop_last_variables(dset, dset->v - orig_nvars);

    return err;
}

/* print critical value for ADF-IPS or KPSS */

static void print_critical_values (double *a, double *cv, 
				   int ci, PRN *prn)
{
    const char *label = N_("Critical values");
    int figs = (ci == ADF)? 2 : 3;

    pprintf(prn, "%*s    ", TRANSLATED_WIDTH(_(label)), " ");
    pprintf(prn, "%g%%      %g%%      %g%%\n", 
	    100*a[0], 100*a[1], 100*a[2]);
    pprintf(prn, "%s: %.*f   %.*f   %.*f\n", 
	    _(label), figs, cv[0], figs, cv[1], figs, cv[2]);
}

static int panel_adjust_ADF_opt (gretlopt *opt)
{
    int err = 0;

    /* Has the user selected an option governing the 
       deterministic terms to be included? If so, don't 
       mess with it.
    */
    if (*opt & (OPT_N | OPT_C | OPT_R | OPT_T)) {
	; /* no-op */
    } else {
	/* panel default: test with constant */
	*opt |= OPT_C;
    }

    return err;
}

static int DF_index (gretlopt opt)
{
    if (opt & OPT_N) {
	return 0;
    } else if (opt & OPT_C) {
	return 1;
    } else if (opt & OPT_T) {
	return 2;
    } else {
	return 3;
    }
}

/* See Im, Pesaran and Shin, "Testing for unit roots in
   heterogeneous panels", Journal of Econometrics 115 (2003),
   53-74.
*/

static int do_IPS_test (double tbar, int n, const int *Ti, 
			int order, const int *Oi,
			gretlopt opt, PRN *prn)
{
    int (*get_IPS_critvals) (int, int, int, double *);
    int (*IPS_tbar_moments) (int, double *, double *);
    int (*IPS_tbar_rho_moments) (int, int, int, double *, double *);
    int Tmin = Ti[1], Tmax = Ti[1];
    int i, T, err = 0;

    for (i=2; i<=n; i++) {
	if (Ti[i] > Tmax) {
	    Tmax = Ti[i];
	}
	if (Ti[i] < Tmin) {
	    Tmin = Ti[i];
	}
    }

    if (Oi != NULL || order > 0) {
	/* non-zero lag order: use IPS's W_{tbar} statistic */
	double E, V, Wtbar, Etbar = 0, Vtbar = 0;
	int order_i;

	IPS_tbar_rho_moments = get_plugin_function("IPS_tbar_rho_moments");

	if (IPS_tbar_rho_moments != NULL) {
	    for (i=0; i<n && !err; i++) {
		T = Ti[i+1];
		order_i = (Oi != NULL)? Oi[i+1] : order;
		err = IPS_tbar_rho_moments(order_i, T, (opt & OPT_T), &E, &V);
		Etbar += E;
		Vtbar += V;
	    }

	    if (!err) {
		Etbar /= n;
		Vtbar /= n;
		Wtbar = sqrt(n) * (tbar - Etbar) / sqrt(Vtbar);
		pprintf(prn, "N = %d, Tmin = %d, Tmax = %d\n", n, Tmin, Tmax);
		pprintf(prn, "Im-Pesaran-Shin W_tbar = %g [%.4f]\n", Wtbar, 
			normal_pvalue_1(-Wtbar));
	    }
	}
    } else if (Tmax > Tmin) {
	/* sample sizes differ: use IPS's Z_{tbar} */
	double E, V, Ztbar, Etbar = 0, Vtbar = 0;

	IPS_tbar_moments = get_plugin_function("IPS_tbar_moments");

	if (IPS_tbar_moments != NULL) {
	    for (i=0; i<n && !err; i++) {
		T = Ti[i+1];
		err = IPS_tbar_moments(T, &E, &V);
		Etbar += E;
		Vtbar += V;
	    }

	    if (!err) {
		Etbar /= n;
		Vtbar /= n;
		Ztbar = sqrt(n) * (tbar - Etbar) / sqrt(Vtbar);
		pprintf(prn, "N = %d, Tmin = %d, Tmax = %d\n", n, Tmin, Tmax);
		pprintf(prn, "Im-Pesaran-Shin Z_tbar = %g [%.4f]\n", Ztbar, 
			normal_pvalue_1(-Ztbar));
	    }
	}
    } else {
	/* simple case: use tbar with exact critical values */
	pprintf(prn, "N,T = (%d,%d)\n", n, Tmax);
	pprintf(prn, "Im-Pesaran-Shin t-bar = %g\n", tbar);

	get_IPS_critvals = get_plugin_function("get_IPS_critvals");

	if (get_IPS_critvals != NULL) {
	    double a[] = { 0.1, 0.05, 0.01 };
	    double cv[3];
		
	    err = (*get_IPS_critvals) (n, Tmax, (opt & OPT_T), cv);
	    if (!err) {
		print_critical_values(a, cv, ADF, prn);
	    }
	}
    }

    return err;
}

/* See In Choi, "Unit root tests for panel data", Journal of
   International Money and Finance 20 (2001), 249-272.
*/

static void do_choi_test (double ppv, double zpv, double lpv, 
			  int n, PRN *prn)
{
    double P = -2 * ppv;
    double Z = zpv / sqrt((double) n);
    int tdf = 5 * n + 4;
    double k = (3.0*tdf)/(M_PI*M_PI*n*(5*n+2));
    double L = sqrt(k) * lpv;

    pprintf(prn, "%s\n", _("Choi meta-tests:"));
    pprintf(prn, "   %s(%d) = %g [%.4f]\n", _("Inverse chi-square"),
	    2*n, P, chisq_cdf_comp(2*n, P));
    pprintf(prn, "   %s = %g [%.4f]\n", _("Inverse normal test"),
	    Z, normal_pvalue_1(-Z));
    pprintf(prn, "   %s: t(%d) = %g [%.4f]\n", _("Logit test"),
	    tdf, L, student_pvalue_1(tdf, -L));
}

static void panel_unit_DF_print (adf_info *ainfo, int i, PRN *prn)
{
    pprintf(prn, "%s %d, T = %d, %s = %d\n", _("Unit"), i, 
	    ainfo->T, _("lag order"), ainfo->order);
    pprintf(prn, "   %s: %g\n"
	    "   %s = %g", 
	    _("estimated value of (a - 1)"), ainfo->b0,
	    _("test statistic"), ainfo->tau);
    if (na(ainfo->pval)) {
	pputs(prn, "\n\n");
    } else {
	pprintf(prn, " [%.4f]\n\n", ainfo->pval);
    }
}

static int panel_DF_test (int v, int order, DATASET *dset, 
			  gretlopt opt, PRN *prn)
{
    int u0 = dset->t1 / dset->pd;
    int uN = dset->t2 / dset->pd;
    int quiet = (opt & OPT_Q);
    int verbose = (opt & OPT_V);
    double ppv = 0.0, zpv = 0.0, lpv = 0.0;
    double pval, tbar = 0.0;
    int *Ti = NULL, *Oi = NULL;
    int i, n, err;

    err = panel_adjust_ADF_opt(&opt);
    if (err) {
	return err;
    }

    if (opt & OPT_G) {
	/* GLS option: can't do IPS t-bar test */
	tbar = NADBL;
    } else {
	Ti = gretl_list_new(uN - u0 + 1);
	if (Ti == NULL) {
	    return E_ALLOC;
	}
	if ((opt & OPT_E) || order < 0) {
	    /* testing down: lag order may vary by unit */
	    Oi = gretl_list_new(uN - u0 + 1);
	    if (Oi == NULL) {
		free(Ti);
		return E_ALLOC;
	    }
	}	    
    }

    if (!quiet) {
	int j = DF_index(opt);

	DF_header(dset->varname[v], (opt & OPT_E)? 0 : order, 
		  0, 0, opt, prn);
	pprintf(prn, "   %s ", _(DF_test_string(j)));
	pputc(prn, '\n');
	pprintf(prn, "   %s: %s\n\n", _("model"), 
		(order > 0)? ADF_model_string(j) : DF_model_string(j));
    }

    /* number of units in sample range */
    n = uN - u0 + 1;

    /* run a Dickey-Fuller test for each unit and record the
       results */

    for (i=u0; i<=uN && !err; i++) {
	adf_info ainfo = {0};

	dset->t1 = i * dset->pd;
	dset->t2 = dset->t1 + dset->pd - 1;
	err = series_adjust_sample(dset->Z[v], &dset->t1, &dset->t2);

	ainfo.v = v;
	ainfo.order = order;
	ainfo.niv = 1;
	ainfo.flags = ADF_PANEL;

	if (!err) {
	    err = real_adf_test(&ainfo, dset, opt, prn);
	    if (!err && verbose) {
		panel_unit_DF_print(&ainfo, i+1, prn);
	    }	    
	}

	if (!err) {
	    if (Ti != NULL) {
		Ti[i-u0+1] = ainfo.T;
	    }
	    if (Oi != NULL) {
		Oi[i-u0+1] = ainfo.order;
	    }	    
	    if (!na(tbar)) {
		tbar += ainfo.tau;
	    }
	    pval = ainfo.pval;
	    if (na(pval)) {
		ppv = zpv = lpv = NADBL;
	    } else if (!na(ppv)) {
		ppv += log(pval);
		zpv += normal_cdf_inverse(pval);
		lpv += log(pval / (1-pval));
	    }
	}
    }

    /* process the results as per Im-Pesaran-Shin and/or Choi */

    if (!err) {
	pprintf(prn, "%s\n\n", _("H0: all groups have unit root"));
	if (!na(tbar)) {
	    tbar /= n;
	    do_IPS_test(tbar, n, Ti, order, Oi, opt, prn);
	}
	if (!na(ppv)) {
	    pputc(prn, '\n');
	    do_choi_test(ppv, zpv, lpv, n, prn);
	}
	pputc(prn, '\n');
    }

    free(Ti);
    free(Oi);

    return err;
}

/**
 * levin_lin_test:
 * @vnum: ID number of variable to test.
 * @plist: list of ADF lag orders.
 * @dset: data information struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Levin-Lin-Chu test
 * for a unit root in panel data. 
 *
 * The list @plist should contain either a single lag order
 * to be applied to all units, or a set of unit-specific
 * orders; in the latter case the length of the list must
 * equal the number of panel units in the current sample
 * range. (This is a gretl list: the first element holds
 * a count of the number of elements following.)
 *
 * By default a test with constant is performed, but the
 * (mutually exclusive) options OPT_N and OPT_T in @opt switch to
 * the case of no constant or constant plus trend respectively.
 * The OPT_Q flag may be used to suppress printed output.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int levin_lin_test (int vnum, const int *plist,
		    DATASET *dset, gretlopt opt, 
		    PRN *prn)
{
    int (*real_levin_lin) (int, const int *, DATASET *, 
			   gretlopt, PRN *);
    int panelmode;
    int err = 0;

    panelmode = multi_unit_panel_sample(dset);

    if (!panelmode || incompatible_options(opt, OPT_N | OPT_T)) {
	return E_BADOPT;
    }

    real_levin_lin = get_plugin_function("real_levin_lin");

    if (real_levin_lin == NULL) {
	fputs(I_("Couldn't load plugin function\n"), stderr);
	err = E_FOPEN;
    } else {
	int save_t1 = dset->t1;
	int save_t2 = dset->t2;
 
	err = (*real_levin_lin) (vnum, plist, dset, opt, prn);

	dset->t1 = save_t1;
	dset->t2 = save_t2;
    }

    return err;
}

/**
 * adf_test:
 * @order: lag order for the (augmented) test.
 * @list: list of variables to test.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Augmented Dickey-Fuller 
 * test for a unit root. 
 *
 * By default two tests are performed, one for a model
 * including a constant and one including a linear trend. The 
 * deterministic components of the model can be controlled via
 * flags in @opt as follows: OPT_N, omit the constant; OPT_C,
 * run just one test using the constant; OPT_T, one test including
 * linear trend; OPT_R, one test including a quadratic trend;
 * OPT_D, include seasonal dummy variables.
 *
 * Additional flags that may be given in @opt include: 
 * OPT_V for verbose operation; OPT_F to apply first-differencing
 * before testing; OPT_G for GLS preprocessing as in Elliott, Rothenberg
 * and Stock (incompatible with OPT_N, OPT_R, OPT_D); OPT_E to
 * "test down" from a given maximum lag order (see the entry for
 * "adf" in the Gretl Command Reference for details).
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int adf_test (int order, const int *list, DATASET *dset, 
	      gretlopt opt, PRN *prn)
{
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int panelmode;
    int err;

    /* GLS incompatible with no const, quadratic trend or seasonals */
    err = incompatible_options(opt, OPT_G | OPT_N | OPT_R);
    if (!err) {
	err = incompatible_options(opt, OPT_D | OPT_G);
    }

    if (!err && (opt & OPT_G)) {
	/* under GLS, have to choose between cases */
	err = incompatible_options(opt, OPT_C | OPT_T);
    }

    panelmode = multi_unit_panel_sample(dset);

    if (panelmode) {
	err = panel_DF_test(list[1], order, dset, opt, prn);
    } else {
	/* regular time series case */
	int i, v, vlist[2] = {1, 0};
	adf_info ainfo = {0};

	ainfo.niv = 1;

	for (i=1; i<=list[0] && !err; i++) {
	    ainfo.v = vlist[1] = list[i];
	    ainfo.order = order;
	    vlist[1] = v = list[i];
	    err = list_adjust_sample(vlist, &dset->t1, &dset->t2, dset, NULL);
	    if (!err && order == -1) {
		/* default to L_{12}: see G. W. Schwert, "Tests for Unit Roots:
		   A Monte Carlo Investigation", Journal of Business and
		   Economic Statistics, 7(2), 1989, pp. 5-17. Note that at
		   some points Ng uses floor(T/100.0) in the following
		   expression, which can give a lower max order.
		*/
		int T = dset->t2 - dset->t1 + 1;

		ainfo.order = 12.0 * pow(T/100.0, 0.25);
	    }
	    if (!err) {
		err = real_adf_test(&ainfo, dset, opt, prn);
	    }

	    dset->t1 = save_t1;
	    dset->t2 = save_t2;
	}
    }

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return err;
}

/* See Peter S. Sephton, "Response surface estimates of the KPSS 
   stationarity test", Economics Letters 47 (1995) 255-261.

   The estimates below of \beta_\infty and \beta_1 (based on
   a bigger replication using Sephton's methodology) allow the
   construction of better critical values for finite samples
   than the values given in the original KPSS article.
*/

static void kpss_parms (double a, int trend, double *b)
{
    const double b0_level[] = { 0.74404, 0.46158, 0.34742 };
    const double b1_level[] = { -0.99120, 0.01642, 0.19814 };
    const double b0_trend[] = { 0.21787, 0.14797, 0.11925 };
    const double b1_trend[] = { -0.25128, 0.03270, 0.10244 };
    int i = (a == .01)? 0 : (a == .05)? 1 : 2;

    if (trend) {
	b[0] = b0_trend[i];
	b[1] = b1_trend[i];
    } else {
	b[0] = b0_level[i];
	b[1] = b1_level[i];
    }	
}

static double kpss_critval (double alpha, int T, int trend)
{
    double b[2];

    kpss_parms(alpha, trend, b);

    return b[0] + b[1]/T;
}

gretl_matrix *kpss_critvals (int T, int trend, int *err)
{
    gretl_matrix *m = NULL;

    if (T < 5) {
	*err = E_TOOFEW;
    } else {
	m = gretl_matrix_alloc(1, 3);
	if (m == NULL) {
	    *err = E_ALLOC;
	} else {
	    m->val[0] = kpss_critval(0.10, T, trend);
	    m->val[1] = kpss_critval(0.05, T, trend);
	    m->val[2] = kpss_critval(0.01, T, trend);
	}
    }

    return m;
}

#define PV_GT10 1.1
#define PV_LT01 -1.0

static double kpss_interp (double s, int T, int trend)
{
    double c10, c05, c01;
    double pv;

    c10 = kpss_critval(.10, T, trend);
    if (s < c10) {
	return PV_GT10;
    }

    c01 = kpss_critval(.01, T, trend);
    if (s > c01) {
	return PV_LT01;
    }  

    /* OK, p-value must lie between .01 and .10 */

    c05 = kpss_critval(.05, T, trend);
    if (s > c05) {
	pv = .01 + .04 * (c01 - s) / (c01 - c05);
    } else {
	pv = .05 + .05 * (c05 - s) / (c05 - c10);
    }

    return pv;
}

static int 
real_kpss_test (int order, int varno, DATASET *dset, 
		gretlopt opt, kpss_info *kinfo, 
		PRN *prn)
{
    MODEL KPSSmod;
    int *list = NULL;
    int hastrend = 0, hasseas = 0;
    double et, s2 = 0.0;
    double cumsum = 0.0, cumsum2 = 0.0;
    double teststat, pval = NADBL;
    double *autocov;
    int t1, t2, T;
    int i, t, ndum, nreg;
    int err = 0;

    /* sanity check */
    if (varno <= 0 || varno >= dset->v) {
	return E_DATA;
    }

    if (gretl_isconst(dset->t1, dset->t2, dset->Z[varno])) {
	gretl_errmsg_sprintf(_("%s is a constant"), dset->varname[varno]);
	return E_DATA;
    }

    if (opt & OPT_F) {
	/* difference the variable before testing */
	varno = diffgenr(varno, DIFF, dset);
	if (varno < 0) {
	    return E_DATA;
	}
    }

    if (opt & OPT_T) {
	hastrend = 1;
    }

    if (opt & OPT_D) {
	hasseas = 1;
    }

    ndum = hasseas ? (dset->pd - 1) : 0;
    nreg = 1 + hastrend + ndum;

    list = gretl_list_new(nreg + 1);
    if (list == NULL) {
	return E_ALLOC;
    }

    list[1] = varno;
    list[2] = 0;

    if (hastrend) {
	list[3] = gettrend(dset, 0);
	if (list[3] == 0) {
	    return E_ALLOC;
	}
    }

    if (hasseas) {
	int *slist = NULL;

	slist = seasonals_list(dset, dset->pd, 0, &err);
	if (err) {
	    free(list);
	    return err;
	} else {
	    for (i=0; i<ndum; i++) {
		list[3 + hastrend + i] = slist[i+1];
	    }
	    free(slist);
	}
    }

    /* OPT_M: reject missing values within sample range */
    KPSSmod = lsq(list, dset, OLS, OPT_A | OPT_M);
    if (KPSSmod.errcode) {
	clear_model(&KPSSmod);
	return KPSSmod.errcode;
    }

    t1 = KPSSmod.t1;
    t2 = KPSSmod.t2;
    T = KPSSmod.nobs;

    if (order < 0) {
	order = 4.0 * pow(T / 100.0, 0.25);
    }

    if (kinfo == NULL && (opt & OPT_V)) {
	KPSSmod.aux = AUX_KPSS;
	printmodel(&KPSSmod, dset, OPT_NONE, prn);
    }
  
    autocov = malloc(order * sizeof *autocov);
    if (autocov == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<order; i++) {
	autocov[i] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
	et = KPSSmod.uhat[t];
	if (na(et)) {
	    continue;
	}
	cumsum += et;
	cumsum2 += cumsum * cumsum;
	s2 += et * et;
	for (i=0; i<order; i++) {
	    int s = i + 1;

	    if (t - s >= t1) {
		autocov[i] += et * KPSSmod.uhat[t - s];
	    }
	}
#if KPSS_DEBUG
	fprintf(stderr, "%d: %#12.4g %#12.4g %#12.4g %#12.4g \n", 
		t, et, KPSSmod.uhat[t-1], s2, cumsum2);
#endif
    }

    for (i=0; i<order; i++) {
	double wt = 1.0 - ((double) (i + 1)) / (order + 1);

	s2 += 2.0 * wt * autocov[i];
    }

    s2 /= T;

    if (s2 <= 0.0) {
	teststat = pval = NADBL;
    } else {
	teststat = cumsum2 / (s2 * T * T);
	pval = kpss_interp(teststat, T, hastrend);
    }

    if (kinfo != NULL) {
	/* storing info for panel test */
	kinfo->T = T;
	kinfo->test = teststat;
	kinfo->pval = pval;
    } else {
	/* testing individual time series */
	if (opt & OPT_V) {
	    pprintf(prn, "  %s: %g\n", _("Robust estimate of variance"), s2);
	    pprintf(prn, "  %s: %g\n", _("Sum of squares of cumulated residuals"), 
		    cumsum2);
	}
    }

    if (!(opt & OPT_Q)) {
	double a[] = { 0.1, 0.05, 0.01 };
	double cv[3];

	cv[0] = kpss_critval(a[0], T, hastrend);
	cv[1] = kpss_critval(a[1], T, hastrend);
	cv[2] = kpss_critval(a[2], T, hastrend);

	pprintf(prn, _("\nKPSS test for %s"), dset->varname[varno]);
	if (hastrend) {
	    if (hasseas) {
		pputs(prn, _(" (including trend and seasonals)\n\n"));
	    } else {
		pputs(prn, _(" (including trend)\n\n"));
	    }
	} else {
	    if (hasseas) {
		pputs(prn, _(" (including seasonals)\n\n"));
	    } else {
		pputs(prn, "\n\n");
	    }
	}

	pprintf(prn, "T = %d\n", T);
	pprintf(prn, _("Lag truncation parameter = %d\n"), order);
	pprintf(prn, "%s = %g\n\n", _("Test statistic"), teststat);
	print_critical_values(a, cv, KPSS, prn);
	if (pval == PV_GT10) {
	    pprintf(prn, "%s > .10\n", _("P-value"));
	} else if (pval == PV_LT01) {
	    pprintf(prn, "%s < .01\n", _("P-value"));
	} else if (!xna(pval)) {
	    pprintf(prn, "%s %.3f\n", _("Interpolated p-value"), pval);
	}
	pputc(prn, '\n');
    }

    if (kinfo == NULL) {
	if (pval == PV_GT10 || pval == PV_LT01) {
	    /* invalidate for record_test_result */
	    pval = NADBL;
	}
	record_test_result(teststat, pval, "KPSS");
    }

    clear_model(&KPSSmod);
    free(list);
    free(autocov);

    return err;
}

static int panel_kpss_test (int order, int v, DATASET *dset, 
			    gretlopt opt, PRN *prn)
{
    kpss_info kinfo;
    int u0 = dset->t1 / dset->pd;
    int uN = dset->t2 / dset->pd;
    int n = uN - u0 + 1;
    int verbose = (opt & OPT_V);
    double ppv = 0.0, zpv = 0.0, lpv = 0.0;
    int gt_10 = 0, lt_01 = 0;
    double pval;
    int i, err = 0;

    /* run a KPSS test for each unit and record the
       results */

    pprintf(prn, _("\nKPSS test for %s %s\n"), dset->varname[v],
	    (opt & OPT_T)? _("(including trend)") : _("(without trend)"));
    pprintf(prn, _("Lag truncation parameter = %d\n"), order);
    pputc(prn, '\n');

    for (i=u0; i<=uN && !err; i++) {
	dset->t1 = i * dset->pd;
	dset->t2 = dset->t1 + dset->pd - 1;
	err = series_adjust_sample(dset->Z[v], &dset->t1, &dset->t2);
	if (!err) {
	    err = real_kpss_test(order, v, dset, opt | OPT_Q, &kinfo, prn);
	    if (!err && verbose) {
		pprintf(prn, "Unit %d, T = %d\n", i + 1, kinfo.T);
		if (na(kinfo.pval)) {
		    pputs(prn, "\n\n");
		} else {
		    pprintf(prn, "test = %g, ", kinfo.test);
		    if (kinfo.pval == PV_GT10) {
			pprintf(prn, "%s > .10\n", _("p-value"));
		    } else if (kinfo.pval == PV_LT01) {
			pprintf(prn, "%s < .01\n", _("p-value"));
		    } else {
			pprintf(prn, "%s %.3f\n", _("interpolated p-value"), 
				kinfo.pval);
		    }
		    pputc(prn, '\n');
		}
	    }
	}

	if (!err) {
	    pval = kinfo.pval;

	    if (pval == PV_GT10) {
		gt_10++;
		if (lt_01 == 0) {
		    /* record lower bound */
		    pval = .10;
		} else {
		    pval = NADBL;
		}
	    } else if (pval == PV_LT01) {
		lt_01++;
		if (gt_10 == 0) {
		    /* record upper bound */
		    pval = .01;
		} else {
		    pval = NADBL;
		}
	    }

	    if (xna(pval)) {
		ppv = zpv = lpv = NADBL;
	    } else if (!na(ppv)) {
		ppv += log(pval);
		zpv += normal_cdf_inverse(pval);
		lpv += log(pval / (1-pval));
	    }
	}
    }

    if (!err && !na(ppv)) {
	/* process the results as per Choi, as best we can */
	pprintf(prn, "%s\n\n", _("H0: all groups are stationary"));
	do_choi_test(ppv, zpv, lpv, n, prn);
	if (gt_10 > 0) {
	    pputs(prn, "   Note: these are LOWER BOUNDS "
		  "on the true p-values\n");
	    pprintf(prn, "   (Individual p-values > .10, and recorded as .10: %d)\n",
		    gt_10);
	} else if (lt_01 > 0) {
	    pputs(prn, "   Note: these are UPPER BOUNDS "
		  "on the true p-values\n");
	    pprintf(prn, "   (Individual p-values < .01, and recorded as .01: %d)\n",
		    lt_01);
	} 
	pputc(prn, '\n');
    } else {
	pprintf(prn, "Choi test: cannot be calculated\n");
    }

    return err;
}

/**
 * kpss_test:
 * @order: window size for Bartlett smoothing.
 * @list: list of variables to test.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the KPSS test for 
 * stationarity. Flags that may be given in @opt include:
 * OPT_T to include a linear trend; OPT_F to apply 
 * first-differencing before testing; OPT_V for verbose
 * operation.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int kpss_test (int order, const int *list, DATASET *dset, 
	       gretlopt opt, PRN *prn)
{
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int orig_nvars = dset->v;
    int err = 0;

    if (multi_unit_panel_sample(dset)) {
	err = panel_kpss_test(order, list[1], dset, opt, prn);
    } else {
	/* regular time series case */
	int i, v, vlist[2] = {1, 0};

	for (i=1; i<=list[0] && !err; i++) {
	    v = list[i];
	    vlist[1] = v;
	    err = list_adjust_sample(vlist, &dset->t1, &dset->t2, 
				     dset, NULL);
	    if (!err) {
		err = real_kpss_test(order, v, dset, opt, NULL, prn);
	    }
	    dset->t1 = save_t1;
	    dset->t2 = save_t2;
	}
    }

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    /* added 2012-03-22 for consistency with adf test */
    dataset_drop_last_variables(dset, dset->v - orig_nvars);

    return err;
}

static int *make_coint_list (const int *list, int detcode, int *nv, 
			     DATASET *dset, int *err)
{
    int *clist = NULL;
    int ifc = 0;
    int i, j = 1;

    /* does the incoming list contain a constant? */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) {
	    ifc = 1;
	    break;
	}
    }

    /* check for sufficient arguments */
    *nv = list[0] - ifc;
    if (*nv < 2) {
	*err = E_ARGS;
	return NULL;
    }

    /* allocate list for cointegrating regression */
    clist = gretl_list_new(*nv + detcode - 1);
    if (clist == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* transcribe original vars */
    for (i=1; i<=list[0]; i++) {
	if (list[i] != 0) {
	    clist[j++] = list[i];
	}
    }    

    /* add trend, if wanted */
    if (detcode >= UR_TREND) {
	clist[j] = gettrend(dset, 0);
	if (clist[j++] == 0) {
	    *err = E_ALLOC;
	} 
    }

    /* add trend-squared, if wanted */
    if (!*err && detcode == UR_QUAD_TREND) {
	clist[j] = gettrend(dset, 1);
	if (clist[j++] == 0) {
	    *err = E_ALLOC;
	}
    }

    /* add const, if wanted */
    if (!*err && detcode != UR_NO_CONST) {
	clist[j] = 0;
    } 

    return clist;
}

static int 
coint_check_opts (gretlopt opt, int *detcode, gretlopt *adf_opt)
{
    if (opt & OPT_N) {
	if ((opt & OPT_T) || (opt & OPT_R)) {
	    return E_BADOPT;
	} else {
	    *detcode = UR_NO_CONST;
	    *adf_opt = OPT_N;
	}
    } else if (opt & OPT_T) {
	if (opt & OPT_R) {
	    return E_BADOPT;
	} else {
	    *detcode = UR_TREND;
	    *adf_opt = OPT_T;
	}
    } else if (opt & OPT_R) {
	*detcode = UR_QUAD_TREND;
	*adf_opt = OPT_R;
    }

    if (opt & OPT_E) {
	*adf_opt |= OPT_E;
    }

    return 0;
}

/* Engle-Granger: try to ensure a uniform sample for the individual
   (A)DF tests and the test on the cointegrating regression
*/

static int coint_set_sample (const int *list, int nv, int order,
			     DATASET *dset)
{
    int anymiss;
    int i, v, t;

    for (t=dset->t1; t<dset->t2; t++) {
	anymiss = 0;
	for (i=1; i<=nv; i++) {
	    v = list[i];
	    if (na(dset->Z[v][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (!anymiss) {
	    break;
	}
    }

    dset->t1 = t + order + 1;
    
    for (t=dset->t2; t>dset->t1; t--) {
	anymiss = 0;
	for (i=1; i<=nv; i++) {
	    v = list[i];
	    if (na(dset->Z[v][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (!anymiss) {
	    break;
	}
    }

    dset->t2 = t;

    return 0;
}

#define EG_MIN_SAMPLE 0

/**
 * engle_granger_test:
 * @order: lag order for the test.
 * @list: specifies the variables to use.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out the Engle-Granger test for cointegration. 
 * Flags that may be given in @opt include: OPT_N, do
 * not an include a constant in the cointegrating regression;
 * OPT_T include constant and linear trend; OPT_R, include
 * quadratic trend; OPT_S, skip DF tests for individual variables; 
 * OPT_E, test down from maximum lag order (see the entry for
 * "adf" in the Gretl Command Reference for details); OPT_V,
 * verbose operation.
 *
 * Returns: 0 on successful completion, non-zero code
 * on error.
 */

int engle_granger_test (int order, const int *list, DATASET *dset, 
			gretlopt opt, PRN *prn)
{
#if EG_MIN_SAMPLE
    int test_t1, test_t2;
#endif
    int orig_t1 = dset->t1;
    int orig_t2 = dset->t2;
    adf_info ainfo = {0};
    gretlopt adf_opt = OPT_C;
    MODEL cmod;
    int detcode = UR_CONST;
    int i, nv, k = 0;
    int step = 1;
    int skip = 0;
    int silent = 0;
    int *clist = NULL;
    int err = 0;

    if (multi_unit_panel_sample(dset)) {
	gretl_errmsg_set("Sorry, this command is not yet available "
			 "for panel data");
	return E_DATA;
    }

    err = coint_check_opts(opt, &detcode, &adf_opt);
    if (err) {
	return err;
    }

    clist = make_coint_list(list, detcode, &nv, dset, &err);
    if (err) {
	return err;
    }

    /* backward compatibility: let a negative lag order
       indicate that we should test down */
    if (order < 0) {
	order = -order;
	adf_opt |= OPT_E;
    }

    /* verbosity? */
    if (opt & OPT_V) {
	adf_opt |= OPT_V;
    }

    /* or silence? */
    if (opt & OPT_I) {
	adf_opt |= OPT_I;
	silent = skip = 1;
    } else if (opt & OPT_S) {
	skip = 1;
    }

    gretl_model_init(&cmod, dset);

    if (!skip) {
	/* start by testing all candidate vars for unit root */
	int uniform_order = (opt & OPT_E)? 0 : order;
	
	coint_set_sample(clist, nv, uniform_order, dset);
	for (i=1; i<=nv; i++) {
	    ainfo.v = clist[i];
	    ainfo.order = order;
	    ainfo.niv = 1;
	    ainfo.flags = ADF_EG_TEST;
	    
	    if (step == 1) {
		pputc(prn, '\n');
	    }
	    pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
		    step++, dset->varname[ainfo.v]);
	    real_adf_test(&ainfo, dset, adf_opt, prn);
	}
    }

    if (!silent) {
	if (step == 1) {
	    pputc(prn, '\n');
	}
	pprintf(prn, _("Step %d: cointegrating regression\n"), step++);
    }

#if EG_MIN_SAMPLE
    test_t1 = dset->t1;
    test_t2 = dset->t2;
#endif

    dset->t1 = orig_t1;
    dset->t2 = orig_t2;

    cmod = lsq(clist, dset, OLS, OPT_NONE);
    err = cmod.errcode;
    if (err) {
	goto bailout;
    }

    if (!silent) {
	cmod.aux = AUX_COINT;
	printmodel(&cmod, dset, OPT_NONE, prn);
    }

    /* add residuals from cointegrating regression to data set */
    err = dataset_add_allocated_series(dset, cmod.uhat);
    if (err) {
	goto bailout;
    }

    k = dset->v - 1;
    strcpy(dset->varname[k], "uhat");
    cmod.uhat = NULL;

    if (!silent) {
	pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
		step, dset->varname[k]);
    }

#if EG_MIN_SAMPLE
    dset->t1 = test_t1;
    dset->t2 = test_t2;
#endif

    ainfo.v = k;
    ainfo.order = order;
    ainfo.niv = nv;
    ainfo.flags = ADF_EG_TEST | ADF_EG_RESIDS;    

    /* Run (A)DF test on the residuals */
    real_adf_test(&ainfo, dset, adf_opt, prn); 

    if (!silent) {
	pputs(prn, _("\nThere is evidence for a cointegrating relationship if:\n"
		     "(a) The unit-root hypothesis is not rejected for the individual"
		     " variables, and\n(b) the unit-root hypothesis is rejected for the "
		     "residuals (uhat) from the \n    cointegrating regression.\n"));
	pputc(prn, '\n');
    }

 bailout:
    
    clear_model(&cmod);
    free(clist);
    if (k > 0) {
	dataset_drop_variable(k, dset);
    }

    dset->t1 = orig_t1;
    dset->t2 = orig_t2;

    return err;
}
