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
#include "usermat.h"
#include "uservar.h"
#include "matrix_extra.h"
#include "gretl_midas.h"

#define MIDAS_DEBUG 0
#define FC_DEBUG 0

struct midas_info_ {
    char lnam0[VNAMELEN];  /* name of MIDAS list on input */
    char lname[VNAMELEN];  /* name of MIDAS list */
    char mname[VNAMELEN];  /* name of initial theta vector */
    int prelag;            /* input list already holds lags */
    int minlag;            /* minimum lag */
    int maxlag;            /* maximum lag */
    int type;              /* type of parameterization */
    int nparm;             /* number of parameters */
    int nterms;            /* number of lag terms */
    int *laglist;          /* list of lag series */
};

typedef struct midas_info_ midas_info;

static midas_info *minfo_from_array (gretl_array *A,
				     int *nmidas,
				     int *err);

/* days per month or quarter: maybe make this user-
   configurable? */

int midas_days_per_period (int days_per_week, int pd)
{
    int ret;
    
    if (days_per_week == 5) {
	ret = 22;
    } else if (days_per_week == 6) {
	ret = 26;
    } else {
	ret = 30;
    }

    return (pd == 12)? ret : 3 * ret;
}

/* Could @m be a valid frequency ratio (= number of members
   of a valid "MIDAS list"), given the periodicity of @dset?
   If so, return 1; if not, return 0.
*/

int is_valid_midas_frequency_ratio (const DATASET *dset, int m)
{
    if (dset->pd == 1) {
	/* lf = annual, hf = quarterly or monthly */
	return (m == 4 || m == 12);
    } else if (dset->pd == 4 && m == 3) {
	/* lf = quarterly, hf = monthly */
	return 1;
    } else if (dset->pd == 4 || dset->pd == 12) {
	/* lf = quarterly or monthly, hf = daily */
	if (m == midas_days_per_period(5, dset->pd)) {
	    return 1;
	} else if (m == midas_days_per_period(6, dset->pd)) {
	    return 1;
	} else if (m == midas_days_per_period(7, dset->pd)) {
	    return 1;
	}
    }

    return 0;
}

/* Infer the current month from the current @qtr along
   with the number of days per period, @ndays, and the
   current index within the month-days array, @day.
*/

static int quarter_to_month (int qtr, int ndays, int day)
{
    return qtr * 3 - 2 + (day - 1) / (ndays/3);
}

/* Construct an auxiliary dataset in which the data from
   a MIDAS list are represented at their "native"
   frequency. We use this for pretty-printing a high
   frequency series.
*/

DATASET *midas_aux_dataset (const int *list,
			    const DATASET *dset,
			    int *err)
{
    DATASET *mset = NULL;
    gretlopt opt = 0;
    int mpd, pd = dset->pd;
    int T, m = list[0];
    int yr, mon;
    int daily = 0;

    if (m < 3 || gretl_list_has_separator(list)) {
	*err = E_INVARG;
    } else if (!dataset_is_time_series(dset)) {
	*err = E_INVARG;
    } else if (pd != 1 && pd != 4 && pd != 12) {
	/* host dataset should be annual, quarterly or monthly */
	*err = E_PDWRONG;
    }

    if (*err) {
	return NULL;
    }

    if (pd == 1) {
	/* annual: midas series should be quarterly or monthly */
	if (m != 4 && m != 12) {
	    *err = E_INVARG;
	} else {
	    mpd = m;
	}
    } else if (pd == 4) {
	/* quarterly: midas series should be monthly or daily */
	if (m == 3) {
	    mpd = 12;
	} else if (m == midas_days_per_period(5, 4)) {
	    mpd = 5;
	} else if (m == midas_days_per_period(6, 4)) {
	    mpd = 6;
	} else if (m == midas_days_per_period(7, 4)) {
	    mpd = 7;
	} else {
	    *err = E_INVARG;
	}
    } else {
	/* monthly: midas series should be daily */
	if (m == midas_days_per_period(5, 12)) {
	    mpd = 5;
	} else if (m == midas_days_per_period(6, 12)) {
	    mpd = 6;
	} else if (m == midas_days_per_period(7, 12)) {
	    mpd = 7;
	} else {
	    *err = E_INVARG;
	}
    }

    if (*err) {
	return NULL;
    }    

    if (!gretl_is_midas_list(list, dset)) {
	gretl_warnmsg_set("The argument does not seem to be a MIDAS list");
    }

    T = sample_size(dset) * m;

    if (mpd >= 5 && mpd <= 7) {
	/* we'll add markers for daily dates */
	daily = 1;
	opt = OPT_M;
    }

    mset = create_auxiliary_dataset(1, T, opt);
    if (mset == NULL) {
	*err = E_ALLOC;
    }

    if (!*err) {
	char *p, obs[OBSLEN];
	int nonex, qtr = 0;
	int i, t, s, m3 = 0;

	mset->pd = mpd;
	mset->structure = TIME_SERIES;
	strcpy(mset->varname[0], dset->varname[list[1]]);
	p = strrchr(mset->varname[0], '_');
	if (p != NULL) *p = '\0';

	ntodate(obs, dset->t1, dset);

	if (mpd == 4) {
	    sprintf(mset->stobs, "%d:1", atoi(obs));
	} else if (mpd == 12) {
	    sprintf(mset->stobs, "%d:01", atoi(obs));
	}

	if (daily && pd == 4) {
	    m3 = m / 3;
	}

	/* loop across observations in low-frequency dataset */

	s = 0;
	for (t=dset->t1; t<=dset->t2; t++) {
	    if (daily) {
		ntodate(obs, t, dset);
		sscanf(obs, "%d:%d", &yr, &mon);
		if (pd == 4) {
		    qtr = mon;
		}
	    }
	    /* read data right-to-left! */
	    for (i=m; i>0; i--) {
		int vi = list[i];
		
		if (daily) {
		    if (pd == 4) {
			mon = quarter_to_month(qtr, m, m-i+1);
			nonex = daily_index_to_date(mset->S[s], yr, mon,
						    (m-i) % m3, mpd);
		    } else {
			nonex = daily_index_to_date(mset->S[s], yr, mon,
						    m-i, mpd);
		    }
		    if (nonex) {
			/* skip any non-existent daily dates */
			mset->t2 -= 1;
		    } else {
			mset->Z[0][s++] = dset->Z[vi][t];
		    }
		} else {
		    mset->Z[0][s++] = dset->Z[vi][t];
		}
	    }
	}

	if (daily) {
	    strcpy(mset->stobs, mset->S[0]);
	    strcpy(mset->endobs, mset->S[mset->t2]);
	    mset->markers = DAILY_DATE_STRINGS;
	}

	mset->sd0 = get_date_x(mset->pd, mset->stobs);
	if (!daily) {
	    ntodate(mset->endobs, mset->t2, mset);
	}
    }

    return mset;
}

static gretl_matrix *maybe_make_auto_theta (char *name, int i,
					    int ptype,
					    int m1, int m2)
{
    gretl_matrix *theta = NULL;
    int k = 0;
    
    if (!strcmp(name, "null")) {
	/* OK if we know how many parameters are needed? */
	if (ptype == MIDAS_BETA0) {
	    k = 2;
	} else if (ptype == MIDAS_BETAN) {
	    k = 3;
	} else if (ptype == MIDAS_U) {
	    k = m2 - m1 + 1;
	}
    } else if (integer_string(name)) {
	int chk = atoi(name);

	if (chk >= 1 && chk < 100) {
	    k = chk;
	}
    }

    if (k > 0) {
	theta = gretl_zero_matrix_new(k, 1);
	if (theta != NULL) {
	    if (ptype == MIDAS_BETA0) {
		theta->val[0] = 1;
		theta->val[1] = 5;
	    } else if (ptype == MIDAS_BETAN) {
		theta->val[0] = 1;
		theta->val[1] = 1;
		theta->val[2] = 0;
	    }
	    sprintf(name, "theta___%d", i+1);
	    private_matrix_add(theta, name);
	}
    }

#if MIDAS_DEBUG
    gretl_matrix_print(theta, "auto-generated theta");
#endif

    return theta;
}

static int lag_info_from_prelag_list (midas_info *m,
				      const int *list,
				      const DATASET *dset)
{
    int m1 = series_get_lag(dset, list[1]);
    int p1 = series_get_midas_period(dset, list[1]);
    int i, p, maxp = 0;

    if (p1 > 0) {
	for (i=1; i<=list[0]; i++) {
	    p = series_get_midas_period(dset, list[i]);
	    if (p > maxp) {
		maxp = p;
	    }
	}
    }

    if (is_valid_midas_frequency_ratio(dset, maxp)) {
	int hfl = maxp * m1 - (p1 - 1);
	
	m->minlag = hfl;
	m->maxlag = hfl + list[0] - 1;
    } else {
	/* oof! just report 1,... */
	m->minlag = 1;
	m->maxlag = list[0];
    }

    return 0;
}

static void midas_info_init (midas_info *m)
{
    m->lnam0[0] = '\0';
    m->lname[0] = '\0';
    m->mname[0] = '\0';
    m->prelag = 0;
    m->minlag = 0;
    m->maxlag = 0;
    m->type = 0;
    m->nparm = 0;
    m->nterms = 0;
}

/* Parse a particular entry in the incoming array of MIDAS
   specifications. Each entry should look like one of
   the following:

   (1) mds(list, minlag, maxlag, type, theta)
   (2) mds(list, minlag, maxlag, 0)
   (3) mdsl(list, type, theta)
   (4) mdsl(list, 0)

   Forms (1) and (2) imply that lags of the MIDAS terms
   in @list should be auto-generated, while (3) and (4)
   imply that the incoming @list already holds whatever
   lags are wanted.

   Variants (2) and (4), with @type set to 0, are for
   U-MIDAS terms: in this case @theta (a vector to
   initialize the coefficients) is not (always) required.
   These variants cannot be used, however, if the full
   specification mixes both U-MIDAS and other terms.
*/

static int parse_midas_info (const char *s,
			     midas_info *minfo, int i,
			     const DATASET *dset)
{
    midas_info *m = &minfo[i];
    char lname[VNAMELEN];
    char mname[VNAMELEN];
    char fmt[48];
    int n, m1, m2, p;
    int umidas = 0;
    int err = 0;

    midas_info_init(m);

    if (!strncmp(s, "mds(", 4)) {
	/* calling for auto-generated lags */
	s += 4;
	sprintf(fmt, "%%%d[^, ] , %%d , %%d , %%d, %%%d[^) ])",
		VNAMELEN-1, VNAMELEN-1);
	n = sscanf(s, fmt, lname, &m1, &m2, &p, mname);
	if (n == 4 && p == MIDAS_U) {
	    umidas = 1;
	} else if (n != 5) {
	    err = E_PARSE;
	}
    } else if (!strncmp(s, "mdsl(", 5)) {
	/* list already hold lags */
	m->prelag = 1;
	s += 5;
	sprintf(fmt, "%%%d[^, ] , %%d, %%%d[^) ])",
		VNAMELEN-1, VNAMELEN-1);
	n = sscanf(s, fmt, lname, &p, mname);
	if (n == 2 && p == MIDAS_U) {
	    umidas = 1;
	} else if (n != 3) {
	    err = E_PARSE;
	}
	m1 = m2 = 0; /* got no min/max info */
    } else {
	err = E_INVARG;
    }

    if (!err) {
	gretl_matrix *theta = NULL;
	int *list = get_list_by_name(lname);
	int k = 0;

	if (!umidas) {
	    theta = get_matrix_by_name(mname);
	    if (theta == NULL) {
		theta = maybe_make_auto_theta(mname, i, p, m1, m2);
	    }
	}

	if (m->prelag && list == NULL) {
	    err = E_INVARG;
	} else if (!m->prelag && !gretl_is_midas_list(list, dset)) {
	    gretl_errmsg_set("mds(): the first term must be a MIDAS list");
	    err = E_INVARG;
	} else if (m1 > m2) {
	    err = E_INVARG;
	} else if (p < 0 || p >= MIDAS_MAX) {
	    err = E_INVARG;
	} else if (umidas) {
	    if (m->prelag) {
		k = list[0];
	    } else {
		k = m2 - m1 + 1;
	    }
	} else {
	    k = gretl_vector_get_length(theta);
	    if (k < 1 || (p == MIDAS_BETA0 && k != 2) ||
		(p == MIDAS_BETAN && k != 3)) {
		err = E_INVARG;
	    }
	}

	if (!err) {
	    strcpy(m->lnam0, lname);
	    strcpy(m->lname, lname);
	    if (!umidas) {
		strcpy(m->mname, mname);
	    }
	    if (m->prelag) {
		/* scrounge lag info from incoming list */
		lag_info_from_prelag_list(m, list, dset);
	    } else {
		m->minlag = m1;
		m->maxlag = m2;
	    }
	    m->type = p;
	    m->nparm = k;
	}
    }

    return err;
}

/* In case we got any U-MIDAS terms in the specification,
   check for coherence. If there's a mixture of U-MIDAS
   specs and others then every spec must have an
   initializer. If all terms are U-MIDAS it's OK if there
   are no initializers (we'll just use OLS in this case).
   It's never OK to have some specs with initializers
   and others without.
*/

static int umidas_check (midas_info *m, int nmidas,
			 int nu, int *use_ols)
{
    int ntheta = 0;
    int i, err = 0;

    /* how many U-MIDAS coeff initializers do we have? */
    for (i=0; i<nmidas; i++) {
	if (m[i].type == MIDAS_U) {
	    if (m[i].mname[0] != '\0') {
		ntheta++;
	    }
	}
    }
    
    if (nu < nmidas) {
	/* mixing U-MIDAS spec(s) with others */
	if (ntheta != nu) {
	    gretl_errmsg_set("In mixed specifications, U-MIDAS terms "
			     "must have an initializer");
	    err = E_INVARG;
	}
    } else {
	/* all specs U-MIDAS */
	if (ntheta == 0) {
	    /* OK, no initializers */
	    *use_ols = 1;
	} else if (ntheta != nu) {
	    gretl_errmsg_set("U-MIDAS: initializers must be given for "
			     "either all or no terms");
	    err = E_INVARG;
	}
    }

    return err;
}

/* Parse the @spec string, which should contain one or more
   MIDAS specifications. For details on what exactly we're
   looking for, see the comment on parse_midas_info() above.
*/

static int 
parse_midas_specs (const char *spec, const DATASET *dset,
		   midas_info **pm, int *pnspec,
		   int *use_ols)
{
    midas_info *m = NULL;
    const char *s;
    int nspec = 0;
    int umidas = 0;
    int err = 0;

    /* first check: count closing parentheses */
    s = spec;
    while (*s) {
	if (*s == ')') {
	    nspec++;
	}
	s++;
    }

    if (nspec == 0) {
	/* spec is junk */
	err = E_PARSE;
    } else {
	/* allocate info structs */
	m = malloc(nspec * sizeof *m);
	if (m == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* parse and record individual MIDAS specs */
	char test[128];
	const char *p;
	int len, i = 0;

	s = spec;
	for (i=0; i<nspec && !err; i++) {
	    while (*s == ' ') s++;
	    p = s;
	    while (*p && *p != ')') p++;
	    len = p - s + 1;
	    if (len > 127) {
		err = E_PARSE;
	    } else {
		*test = '\0';
		strncat(test, s, len);
		err = parse_midas_info(test, m, i, dset);
		if (!err && m[i].type == MIDAS_U) {
		    umidas++;
		}
		s = p + 1;
	    }		
	}
    }

    if (!err && umidas > 0) {
	err = umidas_check(m, nspec, umidas, use_ols);
    }

    if (err) {
	free(m);
	*pnspec = 0;
    } else {
	*pm = m;
	*pnspec = nspec;
    }

    return err;
}

/* Extract the list of regular (low-frequency) regressors.
   While we're at it, see if any regressors are lags of the
   dependent variable and if so record this fact in @ldv.
*/

static int *make_midas_xlist (const int *list,
			      const DATASET *dset,
			      int *ldv,
			      int *err)
{
    int i, nx = list[0] - 1;
    int xi, yno = list[1];
    int *xlist = NULL;

    if (nx > 0) {
	xlist = gretl_list_new(nx);
	if (xlist == NULL) {
	    *err = E_ALLOC;
	} else {
	    for (i=1; i<=nx; i++) {
		xi = list[i+1];
		xlist[i] = xi;
		if (standard_lag_of(xi, yno, dset) > 0) {
		    *ldv = 1;
		}
	    }
	}
    }

    return xlist;
}

/* Given the name of an incoming MIDAS list plus minlag and
   maxlag values (at this point stored in the midas_info
   structure) build the list of required lags of the
   MIDAS series. Or in case the incoming list already
   includes the required lags, just take a pointer to
   the list.
*/

static int *make_midas_laglist (midas_info *m,
				DATASET *dset,
				int *err)
{
    int *list = get_list_by_name(m->lname);

    if (m->prelag) {
	/* don't copy; and don't free either! */
	if (list == NULL) {
	    fprintf(stderr, "make_midas_laglist, prelag: '%s': no list!\n",
		    m->lname);
	    *err = E_DATA;
	}
	return list;
    } else {
	/* copy, because we're going to modify the list */
	int *lcpy = gretl_list_copy(list);

	if (lcpy == NULL) {
	    *err = E_ALLOC;
	} else {
	    *err = list_laggenr(&lcpy, m->minlag, m->maxlag,
				NULL, dset, lcpy[0], OPT_L);
	}

	return lcpy;
    }
}

static int make_midas_laglists (midas_info *minfo,
				int nmidas,
				DATASET *dset)
{
    int *mlist;
    int i, err = 0;
	
    for (i=0; i<nmidas && !err; i++) {
	mlist = make_midas_laglist(&minfo[i], dset, &err);
	if (!err) {
	    if (!minfo[i].prelag) {
		sprintf(minfo[i].lname, "ML___%d", i+1);
		err = remember_list(mlist, minfo[i].lname, NULL);
	    }
	    minfo[i].nterms = mlist[0];
	    minfo[i].laglist = mlist;
	}
    }

    return err;
}

/* If @list is non-NULL, build a full list of all series
   involved: dependent variable, regular regressors, and
   all lags of MIDAS terms. We want this either for setting
   the usable sample range, or (in the case of U-MIDAS via
   OLS), as the list to pass to lsq().

   If @list is NULL, however, the returned list will just
   contain all the MIDAS lag terms: we use this variant
   when setting up for forecasting.
*/

static int *make_midas_biglist (const int *list,
				midas_info *m,
				int nmidas)
{
    int i, j, nt = 0;
    int *biglist = NULL;

    if (list != NULL) {
	nt = list[0];
    }

    for (j=0; j<nmidas; j++) {
	nt += m[j].nterms;
    }

    biglist = gretl_list_new(nt);

    if (biglist != NULL) {
	int k = 1;

	if (list != NULL) {
	    for (i=1; i<=list[0]; i++) {
		biglist[k++] = list[i];
	    }
	}
	for (j=0; j<nmidas; j++) {
	    for (i=1; i<=m[j].nterms; i++) {
		biglist[k++] = m[j].laglist[i];
	    }
	}
    }

    return biglist;
}

/* If we're doing MIDAS via nls, set the usable
   sample range first. That way we'll know how
   big the X data matrix ought to be.
*/

static int midas_set_sample (const int *list,
			     DATASET *dset,
			     midas_info *m,
			     int nmidas)
{
    int *biglist;
    int err = 0;

    biglist = make_midas_biglist(list, m, nmidas);
    if (biglist == NULL) {
	err = E_ALLOC;
    }
				 
    if (!err) {
	int t1 = dset->t1, t2 = dset->t2;

	err = list_adjust_sample(biglist, &t1, &t2, dset, NULL);
	if (!err) {
	    dset->t1 = t1;
	    dset->t2 = t2;
	}
    }

    free(biglist);

#if MIDAS_DEBUG
    fprintf(stderr, "midas_set_sample: returning %d\n", err);
#endif

    return err;
}

/* estimate U-MIDAS via OLS */

static int umidas_ols (MODEL *pmod,
		       const int *list,
		       DATASET *dset,
		       midas_info *m,
		       int nmidas,
		       gretlopt opt)
{
    int *biglist;
    int j;

    biglist = make_midas_biglist(list, m, nmidas);
    
    if (biglist == NULL) {
	return E_ALLOC;
    } else {
	*pmod = lsq(biglist, dset, OLS, opt | OPT_Z);
	free(biglist);
    }

    for (j=0; j<nmidas; j++) {
	/* we're done with these lists now */
	if (!m[j].prelag) {
	    free(m[j].laglist);
	}
	m[j].laglist = NULL;
    }

    return pmod->errcode;
}

/* For forecasting, put into uservar space a "private"
   matrix containing all the coefficients on MIDAS
   lag terms.
*/

static int push_midas_coeff_array (const MODEL *pmod,
				   const int *xlist,
				   midas_info *m,
				   int nmidas)
{
    gretl_matrix *hfb;
    int i, err;
	
    if (gretl_model_get_int(pmod, "umidas")) {
	/* original estimation was U-MIDAS OLS */
	int nx = xlist[0];
	int nt = pmod->ncoeff - nx;

	hfb = gretl_column_vector_alloc(nt);
	if (hfb == NULL) {
	    return E_ALLOC;
	}
	for (i=0; i<nt; i++) {
	    hfb->val[i] = pmod->coeff[nx + i];
	}
    } else {
	gretl_matrix *mc =
	    gretl_model_get_data(pmod, "midas_coeffs");
	int j, k, nt = 0;

	if (mc == NULL) {
	    return E_DATA;
	}
	for (i=0; i<nmidas; i++) {
	    nt += m[i].nterms;
	}
	hfb = gretl_column_vector_alloc(nt);
	if (hfb == NULL) {
	    return E_ALLOC;
	}
	k = 0;
	for (i=0; i<nmidas; i++) {	    
	    for (j=0; j<m[i].nterms; j++) {
		hfb->val[k++] = gretl_matrix_get(mc, j, i);
	    }
	}
    }

#if FC_DEBUG
    gretl_matrix_print(hfb, "hfb in fcast");
#endif

    err = private_matrix_add(hfb, "hfb___");    

    return err;
}

/* Identify parameterizations which take an additional
   leading coefficient: all but U-MIDAS and the plain
   Almon polynomial specification.
*/
#define takes_coeff(t) (t != MIDAS_U && t != MIDAS_ALMONP)

/* Prepare for generating a MIDAS forecast: this
   is called from midas_fcast() in forecast.c.
*/

int midas_forecast_setup (const MODEL *pmod,
			  DATASET *dset,
			  ForecastMethod method,
			  char **pformula)
{
    gretl_array *mA;
    midas_info *minfo = NULL;
    int *xlist = NULL;
    int *hflist = NULL;
    int nmidas = 0;
    int err = 0;

    mA = gretl_model_get_data(pmod, "midas_info");
    xlist = gretl_model_get_data(pmod, "lfxlist");

    if (mA == NULL || xlist == NULL) {
	err = E_DATA;
    } else {
	minfo = minfo_from_array(mA, &nmidas, &err);
    }

    if (!err) {
	/* reconstitute MIDAS lag-lists */
	err = make_midas_laglists(minfo, nmidas, dset);
    }

    if (!err) {
	/* build and push list of all MIDAS terms */
	hflist = make_midas_biglist(NULL, minfo, nmidas);
	if (hflist == NULL) {
	    err = E_ALLOC;
	} else {
	    /* note: remember_list() copies its argument */
	    err = remember_list(hflist, "HFLFC___", NULL);
	    user_var_privatize_by_name("HFLFC___");
	    free(hflist);
	}
    }

    if (!err) {
	/* build and push vector of all MIDAS coeffs */
	err = push_midas_coeff_array(pmod, xlist,
				     minfo, nmidas);
    }

    if (!err) {
	/* build a string that can be passed to "genr" to
	   calculate fitted values */
	char tmp[64], line[MAXLEN];
	double *b = pmod->coeff;
	int p, yno = pmod->list[1];
	int i, xi, j = 0;
	
	*line = '\0';
	gretl_push_c_numeric_locale();
    
	for (i=1; i<=xlist[0]; i++) {
	    xi = xlist[i];
	    if (xi == 0) {
		sprintf(tmp, "%+.15g", b[j++]);
	    } else {
		/* allow for dynamic formula? */
		p = (method == FC_STATIC)? 0 :
		    standard_lag_of(xi, yno, dset);
		if (p > 0) {
		    sprintf(tmp, "%+.15g*%s(-%d)", b[j++],
			    dset->varname[yno], p);
		} else {
		    sprintf(tmp, "%+.15g*%s", b[j++],
			    dset->varname[xi]);
		}
	    }
	    strcat(line, tmp);
	}
	strcat(line, "+lincomb(HFLFC___,hfb___)");

	gretl_pop_c_numeric_locale();

#if FC_DEBUG
	fprintf(stderr, "formula='%s'\n", line);
#endif
	*pformula = gretl_strdup(line);
    }

    free(minfo);

    return err;
}

static int add_midas_plot_matrix (MODEL *pmod,
				  midas_info *m,
				  const double *b)
{
    gretl_matrix *C = NULL;
    gretl_matrix *theta = NULL;
    gretl_matrix *w = NULL;
    int err = 0;
	
    C = gretl_matrix_alloc(m->nterms, 2);
    err = (C == NULL);

    if (!err && m->type != MIDAS_U) {
	theta = gretl_matrix_alloc(m->nparm, 1);
	err = (theta == NULL);
    }

    if (!err) {
	double ci, hfb = 0;
	int p = m->minlag;
	int i, k = 0;
	
	if (m->type == MIDAS_U) {
	    for (i=0; i<m->nparm; i++) {
		gretl_matrix_set(C, i, 0, b[k++]);
		gretl_matrix_set(C, i, 1, p++);
	    }
	} else {
	    hfb = takes_coeff(m->type) ? b[k++] : 1.0;
	    for (i=0; i<m->nparm; i++) {
		theta->val[i] = b[k++];
	    }		
	    w = midas_weights(m->nterms, theta, m->type, &err);
	    if (!err) {
		for (i=0; i<m->nterms; i++) {
		    ci = hfb * w->val[i];
		    gretl_matrix_set(C, i, 0, ci);
		    gretl_matrix_set(C, i, 1, p++);
		}
	    }
	    gretl_matrix_free(w);
	}
    }

    gretl_matrix_free(theta);

    if (!err) {
	char *cnames[2] = {"coeff", "lag"};
	char **S = strings_array_dup(cnames, 2);

	if (S != NULL) {
	    gretl_matrix_set_colnames(C, S);
	}
	err = gretl_model_set_matrix_as_data(pmod, "midas_coeffs", C);
    } else {
	gretl_matrix_free(C);
    }

    return err;
}

static int model_add_minfo_array (MODEL *pmod,
				  midas_info *minfo,
				  int nmidas)
{
    gretl_array *A;
    int err = 0;
    
    A = gretl_array_new(GRETL_TYPE_BUNDLES, nmidas, &err);

    if (A != NULL) {
	midas_info *m;
	gretl_bundle *b;
	int i;

	for (i=0; i<nmidas && !err; i++) {
	    b = gretl_bundle_new();
	    if (b == NULL) {
		err = E_ALLOC;
	    } else {
		m = &minfo[i];
		gretl_bundle_set_string(b, "lname",  m->lnam0);
		gretl_bundle_set_string(b, "mname",  m->mname);
		gretl_bundle_set_int(b, "prelag", m->prelag);
		gretl_bundle_set_int(b, "minlag", m->minlag);
		gretl_bundle_set_int(b, "maxlag", m->maxlag);
		gretl_bundle_set_int(b, "type",  m->type);
		gretl_bundle_set_int(b, "nparm", m->nparm);
		err = gretl_array_set_bundle(A, i, b, 0);
	    }
	}

	if (err) {
	    gretl_array_destroy(A);
	} else {
	    gretl_model_set_array_as_data(pmod, "midas_info", A);
	}
    }

    return err;
}

static midas_info *minfo_from_array (gretl_array *A,
				     int *nmidas,
				     int *err)
{
    int n = gretl_array_get_length(A);
    midas_info *m, *minfo = NULL;

    if (n == 0) {
	*err = E_DATA;
    } else {
	minfo = malloc(n * sizeof *minfo);
	if (minfo == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (minfo != NULL) {
	gretl_bundle *b;
	int i;

	for (i=0; i<n && !*err; i++) {
	    m = &minfo[i];
	    midas_info_init(m);
	    b = gretl_array_get_bundle(A, i);
	    if (b == NULL) {
		*err = E_DATA;
	    } else {
		strcpy(m->lname, gretl_bundle_get_string(b, "lname", err));
		strcpy(m->mname, gretl_bundle_get_string(b, "mname", err));
		m->prelag = gretl_bundle_get_int(b, "prelag", err);
		m->minlag = gretl_bundle_get_int(b, "minlag", err);
		m->maxlag = gretl_bundle_get_int(b, "maxlag", err);
		m->type = gretl_bundle_get_int(b, "type", err);
		m->nparm = gretl_bundle_get_int(b, "nparm", err);
	    }
	}

	if (*err) {
	    free(minfo);
	    minfo = NULL;
	}
    }

    if (!*err) {
	*nmidas = n;
    }

    return minfo;
}

/* Get the MIDAS model ready for shipping out. What
   exactly we do here depends in part on whether
   estimation was done by NLS or OLS.
*/

static int finalize_midas_model (MODEL *pmod,
				 const int *list,
				 const char *param,
				 const DATASET *dset,
				 midas_info *minfo,
				 int nmidas,
				 int *xlist,
				 int ldepvar,
				 int hfslopes)
{
    int type0, mixed = 0;
    int i, err = 0;

    gretl_model_set_string_as_data(pmod, "midas_spec",
				   gretl_strdup(param));

    if (pmod->ci == OLS) {
	/* @pmod is the result of OLS estimation */
	int vi;
	
	gretl_model_allocate_param_names(pmod, pmod->ncoeff);
	for (i=0; i<pmod->ncoeff; i++) {
	    vi = pmod->list[i+2];
	    gretl_model_set_param_name(pmod, i, dset->varname[vi]);
	}
	gretl_model_set_int(pmod, "umidas", 1);
    } else {
	/* @pmod is the result of NLS estimation */
	free(pmod->depvar);
	pmod->depvar = gretl_strdup(dset->varname[list[1]]);
	free(pmod->list);
	pmod->list = gretl_list_copy(list);
    }

    pmod->ci = MIDASREG;

    /* record list of low-frequency regressors */
    gretl_model_set_list_as_data(pmod, "lfxlist", xlist);

    if (ldepvar) {
	gretl_model_set_int(pmod, "dynamic", 1);
    }

    /* record the (common) type of MIDAS term? */
    type0 = minfo[0].type;
    if (nmidas > 1) {
	for (i=1; i<nmidas; i++) {
	    if (minfo[i].type != type0) {
		mixed = 1;
		break;
	    }
	}
    }
    if (!mixed && type0 > 0) {
	gretl_model_set_int(pmod, "midas_type", type0);
    }	
	
    if (nmidas == 1) {
	/* Record the "gross" MIDAS coefficients, to enable
	   drawing of a plot? We'll do this only if we have
	   a single MIDAS term, which is probably the most
	   most common case. Otherwise it becomes too hard
	   to get the plot right.
	*/
	int nx = list[0] - 1;
	const double *b = pmod->coeff + nx;
	
	add_midas_plot_matrix(pmod, &minfo[0], b);
    }

    if (!err) {
	err = model_add_minfo_array(pmod, minfo, nmidas);
    }

    return err;
}

static gretl_matrix *build_XZ (const gretl_matrix *X,
			       const DATASET *dset,
			       midas_info *minfo,
			       int T, int nmidas,
			       int hfslopes)
{
    gretl_matrix *XZ;
    int nx = 0;

    if (X != NULL) {
	/* allow for the unlikely case of no low-freq 
	   regressors */
	nx = X->cols;
    }
	    
    XZ = gretl_matrix_alloc(T, nx + hfslopes);
    
    if (XZ != NULL) {
	const int *zlist;
	double zti;
	int i, j, k, t;

	if (X != NULL) {
	    /* transcribe low-frequency regressors */
	    memcpy(XZ->val, X->val, T * nx * sizeof(double));
	}

	/* transcribe time-mean of MIDAS terms */
	k = nx * T;
	for (i=0; i<nmidas; i++) {
	    if (!takes_coeff(minfo[i].type)) {
		continue;
	    }
	    zlist = minfo[i].laglist;
	    for (t=0; t<T; t++) {
		zti = 0.0;
		for (j=1; j<=zlist[0]; j++) {
		    zti += dset->Z[zlist[j]][t+dset->t1];
		}
		XZ->val[k++] = zti / zlist[0];
	    }
	}		
    }

    return XZ;
}

/* Define "private" matrices to hold the regular X data
   (MX___) and the vector of coefficients on these data
   (bx___). Also add scalars, bmlc___i, to serve as the
   multipliers on the linear combinations of MIDAS terms
   (high-frequency slopes). While we're here, we'll also
   try running OLS to initialize bx___ and the hf slope
   coefficients.
*/

static int add_midas_matrices (int yno,
			       const int *xlist,
			       const DATASET *dset,
			       midas_info *minfo,
			       int nmidas,
			       int *pslopes)
{
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *c = NULL;
    int hfslopes = 0;
    int init_err = 0;
    int i, T, nx = 0;
    int err = 0;

    T = sample_size(dset);

    if (xlist != NULL) {
	nx = xlist[0];
	X = gretl_matrix_data_subset(xlist, dset,
				     dset->t1, dset->t2,
				     M_MISSING_ERROR,
				     &err);
	if (!err) {
	    err = private_matrix_add(X, "MX___");
	}
	if (!err) {
	    b = gretl_zero_matrix_new(nx, 1);
	    if (b!= NULL) {
		err = private_matrix_add(b, "bx___");
	    } else {
		err = E_ALLOC;
	    }
	}
    }

    if (!err) {
	/* for initialization only */
	y = gretl_column_vector_alloc(T);
	if (y != NULL) {
	    memcpy(y->val, dset->Z[yno] + dset->t1,
		   T * sizeof(double));
	} else {
	    init_err = 1;
	}
    }

    if (!err) {
	/* count the HF slope coeffs */
	for (i=0; i<nmidas && !err; i++) {
	    if (takes_coeff(minfo[i].type)) {
		hfslopes++;
	    }
	}
	/* "full-length" coeff vector */
	c = gretl_zero_matrix_new(nx + hfslopes, 1);
	if (c == NULL) {
	    init_err = 1;
	}
    }

    if (!err && !init_err) {
	gretl_matrix *XZ = NULL;

	if (hfslopes > 0) {
	    XZ = build_XZ(X, dset, minfo, T, nmidas, hfslopes);
	    if (XZ == NULL) {
		/* fallback, ignoring "Z" */
		c->rows = nx;
	    }
	}
	if (XZ != NULL) {
	    init_err = gretl_matrix_ols(y, XZ, c, NULL,
					NULL, NULL);
	} else {
	    init_err = gretl_matrix_ols(y, X, c, NULL,
					NULL, NULL);
	}	    
	gretl_matrix_free(XZ);
    }

#if MIDAS_DEBUG
    if (!err && !init_err) {
	gretl_matrix_print(c, "MIDAS OLS initialization");
    }
#endif

    if (!err) {
	if (!init_err) {
	    /* initialize X coeffs from OLS */
	    for (i=0; i<nx; i++) {
		b->val[i] = c->val[i];
	    }
	}
	if (hfslopes > 0) {
	    /* initialize hf slopes, with fallback to zero */
	    int use_c = !init_err && c->rows > nx;
	    char tmp[16];
	    double bzi;
	
	    for (i=0; i<nmidas && !err; i++) {
		if (takes_coeff(minfo[i].type)) {
		    sprintf(tmp, "bmlc___%d", i+1);
		    bzi = use_c ? c->val[nx+i] : 0.0;
		    err = private_scalar_add(bzi, tmp);
		}
	    }
	}
    }

    gretl_matrix_free(y);
    gretl_matrix_free(c);

    /* we're finished with the per-term laglists now */
    for (i=0; i<nmidas; i++) {
	if (!minfo[i].prelag) {
	    free(minfo[i].laglist);
	}
	minfo[i].laglist = NULL;
    }    

#if MIDAS_DEBUG
    fprintf(stderr, "add_midas_matrices: returning %d\n", err);
#endif

    *pslopes = hfslopes;

    return err;
}

static int put_midas_nls_line (char *line,
			       DATASET *dset,
			       PRN *prn)
{
    int err;
    
#if MIDAS_DEBUG
    /* display on stderr what we're passing to nls */
    fputs(line, stderr);
    fputc('\n', stderr);
#endif
    
    err = nl_parse_line(NLS, line, dset, prn);
    *line = '\0';

    return err;
}

/* Append @pname to @s, ensuring that it's preceded by
   a single space if it's not preceded by a double quote.
*/

static void append_pname (char *s, const char *pname)
{
    char c = s[strlen(s) - 1];
    
    if (c != ' ' && c != '"') {
	strncat(s, " ", 1);
    }
    strcat(s, pname);
}

static void make_pname (char *targ, midas_info *m, int i,
			const DATASET *dset)
{
    if (m->type == MIDAS_NEALMON) {
	sprintf(targ, "Almon%d", i+1);
    } else if (m->type == MIDAS_BETA0 || m->type == MIDAS_BETAN) {
	sprintf(targ, "Beta%d", i+1);
    } else if (m->type == MIDAS_ALMONP) {
	sprintf(targ, "Almon%d", i);
    } else {
	/* U-MIDAS */
	int *list = get_list_by_name(m->lname);

	if (list != NULL) {
	    strcpy(targ, dset->varname[list[i+1]]);
	} else {
	    sprintf(targ, "U-MIDAS%d", i+1);
	}
    }
}

/* Main driver function for built-in MIDAS estimation.
   The actual engine used is either NLS or OLS.
*/

MODEL midas_model (const int *list,
		   const char *param,
		   DATASET *dset,
		   gretlopt opt,
		   PRN *prn)
{
    midas_info *minfo = NULL;
    char tmp[64];
    int *xlist = NULL;
    int i, nmidas = 0;
    int ldepvar = 0;
    int hfslopes = 0;
    int origv = dset->v;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int use_ols = 0;
    MODEL mod;
    int err = 0;

    gretl_model_init(&mod, dset);

    if (param == NULL || *param == '\0') {
	err = E_DATA;
    } else {
	err = parse_midas_specs(param, dset, &minfo,
				&nmidas, &use_ols);
    }

    if (!err) {
	/* build list of regular regressors */
	xlist = make_midas_xlist(list, dset, &ldepvar, &err);
	if (xlist != NULL) {
	    err = remember_list(xlist, "XL___", NULL);
	    user_var_privatize_by_name("XL___");
	}
    }

    if (!err) {
	/* build (or just read) MIDAS lag-lists */
	err = make_midas_laglists(minfo, nmidas, dset);
    }

    if (!err && use_ols) {
	err = umidas_ols(&mod, list, dset, minfo, nmidas, opt);
#if MIDAS_DEBUG	
	fputs("*** U-MIDAS via OLS ***\n", stderr);
#endif	
	goto umidas_finish;
    }

    if (!err) {
	/* determine usable sample range */
	err = midas_set_sample(list, dset, minfo, nmidas);
    }

    if (!err) {
	/* add the required matrices */
	err = add_midas_matrices(list[1], xlist, dset, minfo,
				 nmidas, &hfslopes);
    }

#if MIDAS_DEBUG	
    fputs("*** MIDAS via NLS ***\n\n", stderr);
#endif    

    if (!err) {
	char line[MAXLEN];
	int j = 0;

	/* initial "nls" line */
	sprintf(line, "nls %s = ", dset->varname[list[1]]);
	if (xlist != NULL) {
	    strcat(line, "lincomb(XL___, bx___)");
	}
	for (i=0; i<nmidas; i++) {
	    if (takes_coeff(minfo[i].type)) {
		sprintf(tmp, " + bmlc___%d*mlc___%d", ++j, i+1);
	    } else {
		sprintf(tmp, " + mlc___%d", i+1);
	    }
	    strcat(line, tmp);
	}
	err = put_midas_nls_line(line, dset, prn);

	/* MIDAS series and gradient matrices */
	for (i=0; i<nmidas && !err; i++) {
	    if (minfo[i].type == MIDAS_U) {
		sprintf(line, "series mlc___%d = lincomb(%s, %s)",
			i+1, minfo[i].lname, minfo[i].mname);
	    } else {
		sprintf(line, "series mlc___%d = mlincomb(%s, %s, %d)",
			i+1, minfo[i].lname, minfo[i].mname,
			minfo[i].type);
	    }
	    err = put_midas_nls_line(line, dset, prn);
	    if (!err && minfo[i].type != MIDAS_U) {
		sprintf(line, "matrix mgr___%d = mgradient(%d, %s, %d)",
			i+1, minfo[i].nterms, minfo[i].mname, minfo[i].type);
		err = put_midas_nls_line(line, dset, prn);
	    }
	}

	/* derivatives */
	if (!err && xlist != NULL) {
	    strcpy(line, "deriv bx___ = MX___");
	    err = put_midas_nls_line(line, dset, prn);
	}
	for (i=0; i<nmidas && !err; i++) {
	    if (takes_coeff(minfo[i].type)) {
		sprintf(line, "deriv bmlc___%d = {mlc___%d}", i+1, i+1);
		err = put_midas_nls_line(line, dset, prn);
	    }
	    if (!err) {
		if (takes_coeff(minfo[i].type)) {
		    sprintf(line, "deriv %s = bmlc___%d * {%s} * mgr___%d",
			    minfo[i].mname, i+1, minfo[i].lname, i+1);
		} else if (minfo[i].type == MIDAS_ALMONP) {
		    sprintf(line, "deriv %s = {%s} * mgr___%d",
			    minfo[i].mname, minfo[i].lname, i+1);
		} else {
		    sprintf(line, "deriv %s = {%s}", minfo[i].mname,
			    minfo[i].lname);
		}
		err = put_midas_nls_line(line, dset, prn);
	    }
	}

	/* parameter names */
	if (!err) {
	    strcpy(line, "param_names \"");
	    if (xlist != NULL) {
		for (i=1; i<=xlist[0]; i++) {
		    strcpy(tmp, dset->varname[xlist[i]]);
		    append_pname(line, tmp);
		}
	    }
	    for (i=0; i<nmidas && !err; i++) {
		if (takes_coeff(minfo[i].type)) {
		    if (hfslopes > 1) {
			sprintf(tmp, "HF_slope%d", i+1);
		    } else {
			strcpy(tmp, "HF_slope");
		    }
		    append_pname(line, tmp);
		}
		for (j=0; j<minfo[i].nparm; j++) {
		    make_pname(tmp, &minfo[i], j, dset);
		    append_pname(line, tmp);
		}
	    }
	    strcat(line, "\"");
	    err = put_midas_nls_line(line, dset, prn);
	}
    }

#if MIDAS_DEBUG
    fputc('\n', stderr);
#endif

    if (!err) {
	mod = nl_model(dset, (opt | OPT_G | OPT_M), prn);
    }

 umidas_finish:

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    for (i=0; i<nmidas; i++) {
	if (!minfo[i].prelag) {
	    user_var_delete_by_name(minfo[i].lname, NULL);
	}
	if (minfo[i].type != MIDAS_U) {
	    sprintf(tmp, "mgr___%d", i+1);
	    user_var_delete_by_name(tmp, NULL);
	}
    }

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	finalize_midas_model(&mod, list, param, dset,
			     minfo, nmidas, xlist,
			     ldepvar, hfslopes);
    } else {
	free(xlist);
    }

    destroy_private_uvars();
    free(minfo);

    if (dset->v > origv) {
	/* or maybe not? */
	dataset_drop_last_variables(dset, dset->v - origv);
    }

    return mod;
}
