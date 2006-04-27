#define MAX_ARMA_ORDER 6
#define MAX_ARIMA_DIFF 2

struct arma_info {
    int yno;         /* ID of dependent variable */
    ArmaFlags flags; /* specification flags */
    int ifc;         /* specification includes a constant? */
    int p;           /* non-seasonal AR order */
    int d;           /* non-seasonal difference */
    int q;           /* non-seasonal MA order */
    int P;           /* seasonal AR order */
    int D;           /* seasonal difference */
    int Q;           /* seasonal MA order */
    int maxlag;      /* longest lag in model */
    int nexo;        /* number of other regressors (ARMAX) */
    int nc;          /* total number of coefficients */
    int t1;          /* starting observation */
    int t2;          /* ending observation */
    int pd;          /* periodicity of data */
    int T;           /* full length of data series */
    double *dy;      /* differenced dependent variable */
};

#define arma_has_seasonal(a)   ((a)->flags & ARMA_SEAS)
#define arma_is_arima(a)       ((a)->flags & ARMA_DSPEC)
#define arma_by_x12a(a)        ((a)->flags & ARMA_X12A)
#define arma_exact_ml(a)       ((a)->flags & ARMA_EXACT)

#define set_arma_has_seasonal(a)  ((a)->flags |= ARMA_SEAS)
#define set_arma_is_arima(a)      ((a)->flags |= ARMA_DSPEC)
#define unset_arma_is_arima(a)    ((a)->flags &= ~ARMA_DSPEC)

static void 
arma_info_init (struct arma_info *ainfo, char flags, const DATAINFO *pdinfo)
{
    ainfo->yno = 0;
    ainfo->flags = flags;

    ainfo->p = 0;
    ainfo->d = 0;
    ainfo->q = 0;
    ainfo->P = 0;
    ainfo->D = 0;
    ainfo->Q = 0; 

    ainfo->maxlag = 0;
    ainfo->ifc = 0;
    ainfo->nexo = 0;
    ainfo->nc = 0;

    ainfo->t1 = pdinfo->t1;
    ainfo->t2 = pdinfo->t2;
    ainfo->pd = pdinfo->pd;
    ainfo->T = pdinfo->n;

    ainfo->dy = NULL;
}

static int arma_list_y_position (struct arma_info *ainfo)
{
    int ypos;

    if (arma_is_arima(ainfo)) {
	ypos = (arma_has_seasonal(ainfo))? 9 : 5;
    } else {
	ypos = (arma_has_seasonal(ainfo))? 7 : 4;
    }

    return ypos;
}

#define INT_DEBUG 0

static int arima_integrate (double *dx, const double *x,
			    int t1, int t2, int d, int D, int s)
{
    double *ix;
    int t;

#if INT_DEBUG
    fprintf(stderr, "arima_integrate: t1=%d, t2=%d, d=%d, D=%d, s=%d\n",
	    t1, t2, d, D, s);
#endif

    ix = malloc((t2 + 1) * sizeof *ix);
    if (ix == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<t1; t++) {
	ix[t] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
	ix[t] = dx[t];
	if (d > 0) {
	    ix[t] += x[t-1];
	} 
	if (d > 1) {
	    ix[t] += x[t-1];
	    ix[t] -= x[t-2];
	}
	if (D > 0) {
	    ix[t] += x[t-s];
	    if (d > 0) {
		ix[t] -= x[t-s-1];
	    }
	    if (d > 1) {
		ix[t] -= x[t-s-1];
		ix[t] += x[t-2*s];
	    }	    
	} 
	if (D > 1) {
	    ix[t] += x[t-s];
	    ix[t] -= x[t-2*s];
	    if (d > 0) {
		ix[t] += x[t-s];
		ix[t] -= x[t-s-1];
		ix[t] += x[t-2*s-1];
	    }
	    if (d > 1) {
		ix[t] -= 2 * x[t-s-1];
		ix[t] += 2 * x[t-s-2];
		ix[t] += x[t-2*s-1];
		ix[t] -= x[t-2*s-2];
	    }
	}
    }

#if INT_DEBUG
    for (t=0; t<=t2; t++) {
	fprintf(stderr, "%2d: %12.5g %12.5g %12.5g\n",
		t, x[t], dx[t], ix[t]);
    }
#endif

    /* transcribe integrated result back into "dx" */
    for (t=0; t<=t2; t++) {
	if (t < t1) {
	    dx[t] = NADBL;
	} else {
	    dx[t] = ix[t];
	}
    }

    free(ix);

    return 0;
}

static void ainfo_data_to_model (struct arma_info *ainfo, MODEL *pmod)
{
    pmod->ifc = ainfo->ifc;
    pmod->dfn = ainfo->nc - pmod->ifc;
    pmod->dfd = pmod->nobs - pmod->dfn;
    pmod->ncoeff = ainfo->nc;

    if (arma_has_seasonal(ainfo)) {
	gretl_model_set_int(pmod, "arma_P", ainfo->P);
	gretl_model_set_int(pmod, "arma_Q", ainfo->Q);
	gretl_model_set_int(pmod, "arma_pd", ainfo->pd);	
    }

    if (ainfo->d > 0 || ainfo->D > 0) {
	gretl_model_set_int(pmod, "arima_d", ainfo->d);
	gretl_model_set_int(pmod, "arima_D", ainfo->D);
    }

    if (ainfo->nexo > 0) {
	gretl_model_set_int(pmod, "armax", 1);
    }
}

/* write the various statistics from ARMA estimation into
   a gretl MODEL struct */

static void write_arma_model_stats (MODEL *pmod, const int *list, 
				    struct arma_info *ainfo,
				    const double **Z, 
				    const DATAINFO *pdinfo)
{
    const double *y = NULL;
    double mean_error;
    int t;

    pmod->ci = ARMA;

    ainfo_data_to_model(ainfo, pmod);

    free(pmod->list);
    pmod->list = gretl_list_copy(list);

    if (arma_is_arima(ainfo)) {
	y = ainfo->dy;
    } else {
	y = Z[ainfo->yno];
    }

    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, y);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, y);

    mean_error = pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(y[t]) && !na(pmod->uhat[t])) {
	    pmod->yhat[t] = y[t] - pmod->uhat[t];
	    pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	    mean_error += pmod->uhat[t];
	}
    }

    if (arma_is_arima(ainfo)) {
	arima_integrate(pmod->yhat, Z[ainfo->yno], pmod->t1, pmod->t2, 
			ainfo->d, ainfo->D, ainfo->pd);
    }

    mean_error /= pmod->nobs;
    gretl_model_set_double(pmod, "mean_error", mean_error);

    if (na(pmod->sigma)) {
	/* in X12A or native exact cases this is already done */
	pmod->sigma = sqrt(pmod->ess / pmod->nobs);
    } 

    pmod->rsq = pmod->adjrsq = pmod->fstt = NADBL;
    pmod->tss = NADBL;

    if (!arma_by_x12a(ainfo)) {
	mle_criteria(pmod, 1);
    }

    gretl_model_add_arma_varnames(pmod, pdinfo, ainfo->yno,
				  ainfo->p, ainfo->q, 
				  ainfo->P, ainfo->Q,
				  ainfo->nexo);
}

static void calc_max_lag (struct arma_info *ainfo)
{
    int pmax = ainfo->p;
    int dmax = ainfo->d;

    if (arma_has_seasonal(ainfo)) {
	pmax += ainfo->P * ainfo->pd;
	dmax += ainfo->D * ainfo->pd;
    }

    ainfo->maxlag = pmax + dmax;
}

static int 
arma_adjust_sample (const DATAINFO *pdinfo, const double **Z, const int *list,
		    struct arma_info *ainfo)
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int an, i, v, t, t1min;
    int vstart, pmax, anymiss;

    vstart = arma_list_y_position(ainfo);

    pmax = ainfo->p + ainfo->P * ainfo->pd;

    /* determine starting point for valid data, t1min */
    t1min = 0;
    for (t=0; t<=pdinfo->t2; t++) {
	anymiss = 0;
	for (i=vstart; i<=list[0]; i++) {
	    v = list[i];
	    if (na(Z[v][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (anymiss) {
	    t1min++;
        } else {
	    break;
	}
    }

    if (arma_by_x12a(ainfo) || arma_exact_ml(ainfo)) {
	/* FIXME x12a in conditional mode? */
	;
    } else {
	t1min += ainfo->maxlag;
    }

    /* if the notional starting point is before the start of
       valid data, advance it */
    if (t1 < t1min) {
	t1 = t1min;
    }

    /* trim any missing obs from the end of the specified sample
       range 
    */
    for (t=pdinfo->t2; t>=t1; t--) {
	anymiss = 0;
	for (i=vstart; i<=list[0]; i++) {
	    v = list[i];
	    if (na(Z[v][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (anymiss) {
	    t2--;
        } else {
	    break;
	}
    }

    if (arma_by_x12a(ainfo) || arma_exact_ml(ainfo)) {
	/* FIXME x12a in conditional mode? */
	t1min = t1;
    } else {    
	t1min = t1 - pmax; /* wrong? */
	if (t1min < 0) {
	    t1min = 0;
	}
    }

    /* check for missing obs within the sample range */
    for (t=t1min; t<t2; t++) {
	for (i=vstart; i<=list[0]; i++) {
	    if (t < t1 && i > vstart) {
		continue;
	    }
	    v = list[i];
	    if (na(Z[v][t])) {
		char msg[64];

		sprintf(msg, _("Missing value encountered for "
			       "variable %d, obs %d"), v, t + 1);
		gretl_errmsg_set(msg);
		return 1;
	    }
	}
    }

    an = t2 - t1 + 1;
    if (an <= ainfo->nc) {
	return 1; 
    }

    ainfo->t1 = t1;
    ainfo->t2 = t2;

    return 0;
}

#define ARIMA_DEBUG 0

/* remove the intercept from list of regressors */

static int arma_remove_const (int *list, int seasonal, int diffs,
			      const double **Z, const DATAINFO *pdinfo)
{
    int xstart, ret = 0;
    int i, j;

    if (diffs) {
	xstart = (seasonal)? 10 : 6;
    } else {
	xstart = (seasonal)? 8 : 5;
    }

    for (i=xstart; i<=list[0]; i++) {
	if (list[i] == 0 || true_const(list[i], Z, pdinfo)) {
	    for (j=i; j<list[0]; j++) {
		list[j] = list[j+1];
	    }
	    list[0] -= 1;
	    ret = 1;
	    break;
	}
    }

    return ret;
}

static int check_arma_sep (int *list, int sep1, struct arma_info *ainfo)
{
    int sep2 = (sep1 == 3)? 6 : 8;
    int i, err = 0;

    for (i=sep1+1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    if (i == sep2) {
		/* there's a second list separator in the right place:
		   we've got a seasonal specification */
		set_arma_has_seasonal(ainfo);
	    } else {
		err = 1;
	    }
	}
    }

    if (!err && sep1 == 4) {
	/* check for apparent but not "real" arima spec */
	if (arma_has_seasonal(ainfo)) {
	    if (list[2] == 0 && list[6] == 0) {
		gretl_list_delete_at_pos(list, 2);
		gretl_list_delete_at_pos(list, 5);
		unset_arma_is_arima(ainfo);
	    }
	} else {
	    if (list[2] == 0) {
		gretl_list_delete_at_pos(list, 2);
		unset_arma_is_arima(ainfo);
	    }
	}
    }

#if ARIMA_DEBUG
    fprintf(stderr, "check_arma_sep: returning %d\n", err);
#endif

    return err;
}

static int check_arma_list (int *list, gretlopt opt, 
			    const double **Z, const DATAINFO *pdinfo,
			    struct arma_info *ainfo)
{
    int armax = 0;
    int hadconst = 0;
    int err = 0;

    if (arma_has_seasonal(ainfo)) {
	armax = (list[0] > 7);
    } else {
	armax = (list[0] > 4);
    }

    if (list[1] < 0 || list[1] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[2] < 0 || list[2] > MAX_ARMA_ORDER) {
	err = 1;
    } 

    if (!err) {
	ainfo->p = list[1];
	ainfo->q = list[2];
    }

    if (!err && arma_has_seasonal(ainfo)) {
	if (list[0] < 7) {
	    err = 1;
	} else if (list[4] < 0 || list[4] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[5] < 0 || list[5] > MAX_ARMA_ORDER) {
	    err = 1;
	} 
    }

    if (!err && arma_has_seasonal(ainfo)) {
	ainfo->P = list[4];
	ainfo->Q = list[5];
    }

    /* If there's an explicit constant in the list here, we'll remove
       it, since it is added implicitly later.  But if we're supplied
       with OPT_N (meaning: no intercept) we'll flag this by
       setting ifc = 0.  Also, if the user gave an armax list
       (specifying regressors) we'll respect the absence of a constant
       from that list by setting ifc = 0.
    */

    if (!err) {
	if (armax) {
	    hadconst = arma_remove_const(list, arma_has_seasonal(ainfo),
					 0, Z, pdinfo);
	}
	if ((opt & OPT_N) || (armax && !hadconst)) {
	    ;
	} else {
	    ainfo->ifc = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command"));
    } else {
	ainfo->nexo = list[0] - ((arma_has_seasonal(ainfo))? 7 : 4);
	ainfo->nc = ainfo->p + ainfo->q + ainfo->P + ainfo->Q
	    + ainfo->nexo + ainfo->ifc;
	ainfo->yno = (arma_has_seasonal(ainfo))? list[7] : list[4];
    }

    return err;
}

static int check_arima_list (int *list, gretlopt opt, 
			     const double **Z, const DATAINFO *pdinfo,
			     struct arma_info *ainfo)
{
    int armax = 0;
    int hadconst = 0;
    int err = 0;

#if ARIMA_DEBUG
    printlist(list, "check_arima_list");
#endif

    if (arma_has_seasonal(ainfo)) {
	armax = (list[0] > 9);
    } else {
	armax = (list[0] > 5);
    }

    if (list[1] < 0 || list[1] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[2] < 0 || list[2] > MAX_ARIMA_DIFF) {
	err = 1;
    } else if (list[3] < 0 || list[3] > MAX_ARMA_ORDER) {
	err = 1;
    } 

    if (!err) {
	ainfo->p = list[1];
	ainfo->d = list[2];
	ainfo->q = list[3];
    }

    if (!err && arma_has_seasonal(ainfo)) {
	if (list[0] < 9) {
	    err = 1;
	} else if (list[5] < 0 || list[5] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[6] < 0 || list[6] > MAX_ARIMA_DIFF) {
	    err = 1;
	} else if (list[7] < 0 || list[7] > MAX_ARMA_ORDER) {
	    err = 1;
	} 
    }

    if (!err && arma_has_seasonal(ainfo)) {
	ainfo->P = list[5];
	ainfo->D = list[6];
	ainfo->Q = list[7];
    }

    /* If there's an explicit constant in the list here, we'll remove
       it, since it is added implicitly later.  But if we're supplied
       with OPT_N (meaning: no intercept) we'll flag this by
       setting ifc = 0.  Also, if the user gave an armax list
       (specifying regressors) we'll respect the absence of a constant
       from that list by setting ifc = 0.
    */

    if (!err) {
	if (armax) {
	    hadconst = arma_remove_const(list, arma_has_seasonal(ainfo),
					 1, Z, pdinfo);
	}
	if ((opt & OPT_N) || (armax && !hadconst)) {
	    ;
	} else {
	    ainfo->ifc = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command"));
    } else {
	ainfo->nexo = list[0] - ((arma_has_seasonal(ainfo))? 9 : 5);
	ainfo->nc = ainfo->p + ainfo->q + ainfo->P + ainfo->Q
	    + ainfo->nexo + ainfo->ifc;
	ainfo->yno = (arma_has_seasonal(ainfo))? list[9] : list[5];
    }

    return err;
}

static int arma_check_list (int *list, gretlopt opt,
			    const double **Z, const DATAINFO *pdinfo,
			    struct arma_info *ainfo)
{
    int sep1 = gretl_list_separator_position(list);
    int err = 0;

#if ARIMA_DEBUG
    fprintf(stderr, "arma_check_list: sep1 = %d\n", sep1);
    printlist(list, "incoming list");
#endif

    if (sep1 == 3) {
	if (list[0] < 4) {
	    err = E_PARSE;
	}
    } else if (sep1 == 4) {
	if (list[0] < 5) {
	    err = E_PARSE;
	} else {
	    set_arma_is_arima(ainfo);
	}
    } else {
	err = E_PARSE;
    }

    if (!err) {
	err = check_arma_sep(list, sep1, ainfo);
    }

    if (!err) {
	if (arma_is_arima(ainfo)) {
	    /* check for arima spec */
	    err = check_arima_list(list, opt, Z, pdinfo, ainfo);
	} else {	    
	    /* check for simple arma spec */
	    err = check_arma_list(list, opt, Z, pdinfo, ainfo);
	} 
    }

    /* catch null model */
    if (ainfo->nc == 0) {
	err = E_ARGS;
    }

#if ARIMA_DEBUG
    printlist(list, "ar(i)ma list after checking");
    fprintf(stderr, "err = %d\n", err);
#endif

    return err;
}

static int
arima_difference (const double *y, struct arma_info *ainfo)
{
    double *dy;
    int s = ainfo->pd;
    int t, t1 = 0;

#if ARMA_DEBUG
    fprintf(stderr, "doing arima_difference: d = %d, D = %d\n",
	    ainfo->d, ainfo->D);
#endif

    dy = malloc(ainfo->T * sizeof *dy);
    if (dy == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<ainfo->T; t++) {
	if (na(y[t])) {
	    t1++;
	} else {
	    break;
	}
    }

    t1 += ainfo->d + ainfo->D * s;

    for (t=0; t<t1; t++) {
	dy[t] = NADBL;
    }

    for (t=t1; t<ainfo->T; t++) {
	dy[t] = y[t];
	if (ainfo->d > 0) {
	    dy[t] -= y[t-1];
	} 
	if (ainfo->d > 1) {
	    dy[t] -= y[t-1];
	    dy[t] += y[t-2];
	}
	if (ainfo->D > 0) {
	    dy[t] -= y[t-s];
	    if (ainfo->d > 0) {
		dy[t] += y[t-s-1];
	    }
	    if (ainfo->d > 1) {
		dy[t] += y[t-s-1];
		dy[t] -= y[t-2*s];
	    }	    
	} 
	if (ainfo->D > 1) {
	    dy[t] -= y[t-s];
	    dy[t] += y[t-2*s];
	    if (ainfo->d > 0) {
		dy[t] -= y[t-s];
		dy[t] += y[t-s-1];
		dy[t] -= y[t-2*s-1];
	    }
	    if (ainfo->d > 1) {
		dy[t] += 2 * y[t-s-1];
		dy[t] -= 2 * y[t-s-2];
		dy[t] -= y[t-2*s-1];
		dy[t] += y[t-2*s-2];
	    }
	}
    }

    ainfo->dy = dy;

    return 0;
}

