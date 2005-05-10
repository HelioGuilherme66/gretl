/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* non-parametric stats for gretl */

#include "libgretl.h"

enum {
    RUN_FIRST,
    RUN_REPEAT,
    RUN_CLEANUP
};

static int inverse_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
     
    return (*da < *db) - (*da > *db);
}

/* alpha = .001, .01, .025, .05, .1 */

static double spearman_critical[18][5] = {      /* n = */
    { 0.9643, 0.8571, 0.7450, 0.6786, 0.5357 }, /*  7 */
    { 0.9286, 0.8095, 0.7143, 0.6190, 0.5000 }, /*  8 */
    { 0.9000, 0.7667, 0.6833, 0.5833, 0.4667 }, /*  9 */
    { 0.8667, 0.7333, 0.6364, 0.5515, 0.4424 }, /* 10 */
    { 0.8364, 0.7000, 0.6091, 0.5273, 0.4182 }, /* 11 */
    { 0.8182, 0.6713, 0.5804, 0.4965, 0.3986 }, /* 12 */
    { 0.7912, 0.6429, 0.5549, 0.4780, 0.3791 }, /* 13 */
    { 0.7670, 0.6220, 0.5341, 0.4593, 0.3626 }, /* 14 */
    { 0.7464, 0.6000, 0.5179, 0.4429, 0.3500 }, /* 15 */
    { 0.7265, 0.5824, 0.5000, 0.4265, 0.3382 }, /* 16 */
    { 0.7083, 0.5637, 0.4853, 0.4118, 0.3260 }, /* 17 */
    { 0.6904, 0.5480, 0.4716, 0.3994, 0.3148 }, /* 18 */
    { 0.6737, 0.5333, 0.4579, 0.3895, 0.3070 }, /* 19 */
    { 0.6586, 0.5203, 0.4451, 0.3789, 0.2977 }, /* 20 */
    { 0.6455, 0.5078, 0.4351, 0.3688, 0.2909 }, /* 21 */
    { 0.6318, 0.4963, 0.4241, 0.3597, 0.2829 }, /* 22 */
    { 0.6186, 0.4852, 0.4150, 0.3518, 0.2767 }, /* 23 */
    { 0.6070, 0.4748, 0.4061, 0.3435, 0.2704 }  /* 24 */
};

static double spearman_signif (double rho, int n)
{
    const double *vals = spearman_critical[n - 7];

    if (rho > vals[0]) return .001;
    else if (rho > vals[1]) return .01;
    else if (rho > vals[2]) return .025;
    else if (rho > vals[3]) return .05;
    else if (rho > vals[4]) return .1;
    else return 1.0;
}

static int spearman_rho (const double *x, const double *y, int n, 
			 double *rho, double *sd, double *pval,
			 double **rxout, double **ryout, int *m,
			 int code)
{
    static double *sx = NULL;
    static double *sy = NULL;
    static double *rx = NULL;
    static double *ry = NULL;
    static int nn;

    double xx, yy, r;
    int i, j, t;

    if (code == RUN_CLEANUP) {
	free(sx); sx = NULL;
	free(sy); sy = NULL;
	free(rx); rx = NULL;
	free(ry); ry = NULL;
	nn = 0;
	return 0;
    }

    *rho = NADBL;
    *sd = NADBL;
    if (pval != NULL) {
	*pval = NADBL;
    }

    if (code == RUN_FIRST) {
	sx = malloc(n * sizeof *sx);
	sy = malloc(n * sizeof *sy);
	rx = malloc(n * sizeof *rx);
	ry = malloc(n * sizeof *ry);

	if (sx == NULL || sy == NULL || 
	    rx == NULL || ry == NULL) { 
	    return E_ALLOC;
	}
    }

    /* copy non-missing x and y into sx, sy */
    nn = 0;
    for (t=0; t<n; t++) {
	if (na(x[t]) || na(y[t])) {
	    continue;
	}
	sx[nn] = x[t];
	sy[nn] = y[t];
	nn++;
    }

    /* get sorted series */
    qsort(sx, nn, sizeof *sx, inverse_compare_doubles);
    qsort(sy, nn, sizeof *sy, inverse_compare_doubles);

    for (t=0; t<n; t++) {
	rx[t] = ry[t] = 0.0;
    }

    /* make rankings by comparing "raw" x, y with sorted */

    r = 1.0;
    for (i=0; i<nn; i++) {
	/* scan sorted x */
	int cases = 0, k = 0;

	if (i > 0 && sx[i] == sx[i-1]) {
	    continue;
	}

	for (j=0; j<n; j++) {
	    /* scan raw x for matches */
	    if (na(x[j])) {
		continue;
	    }
	    if (x[j] == sx[i]) {
		rx[k] = r;
		cases++;
	    }
	    k++;
	}
	if (cases > 1) {
	    double avg = (r + r + cases - 1.0) / 2.0;

	    for (t=0; t<nn; t++) {
		if (rx[t] == r) {
		    rx[t] = avg;
		}
	    }
	} 
	r += cases;
    }
		
    r = 1.0;
    for (i=0; i<nn; i++) {
	/* scan sorted y */
	int cases = 0, k = 0;

	if (i > 0 && sy[i] == sy[i-1]) {
	    continue;
	}

	for (j=0; j<n; j++) {
	    /* scan raw y for matches */
	    if (na(y[j])) {
		continue;
	    }
	    if (y[j] == sy[i]) {
		ry[k] = r;
		cases++;
	    }
	    k++;
	}
	if (cases > 1) {
	    double avg = (r + r + cases - 1.0) / 2.0;

	    for (t=0; t<nn; t++) {
		if (ry[t] == r) {
		    ry[t] = avg;
		}
	    }
	} 
	r += cases;
    }	    
	    
    /* calculate rho and standard error */
    xx = 0.0;
    for (i=0; i<nn; i++) { 
	xx += (rx[i] - ry[i]) * (rx[i] - ry[i]);
    }
    yy = 1.0 - 6.0 * xx / (nn * (nn * nn - 1));
    xx = sqrt(1.0 / (nn - 1.0));

    *rho = yy;
    *sd = xx;

    if (nn >= 20 && pval != NULL) {
	*pval = normal_pvalue_1(fabs(yy / xx)); 
    } 

    /* save the ranks, if wanted */
    if (rxout != NULL) {
	*rxout = rx;
    } 
    if (ryout != NULL) {
	*ryout = ry;
    } 

    if (m != NULL) {
	*m = nn;
    }

    return 0;
}

/**
 * spearman:
 * @list: list of (two) variables to process.
 * @Z: data matrix.
 * @pdinfo: information on the data set.
 * @opt: if non-zero, print both the "raw" and the ranked data.
 * @prn: gretl printing struct.
 *
 * Calculates and prints Spearman's rank correlation coefficient for the two
 * variables specified in the @list.
 * 
 * Returns: 0 on successful completion, 1 on error.
 */

int spearman (const int *list, const double **Z, const DATAINFO *pdinfo, 
	      gretlopt opt, PRN *prn)
{
    double *rx = NULL, *ry = NULL;
    double rho, sd, pval;
    int vx, vy, t, m;

    if (list[0] != 2) {
	pputs(prn, _("spearman command requires two variables\n"));
	return 1;
    }

    vx = list[1];
    vy = list[2];

    if (spearman_rho(Z[vx] + pdinfo->t1, 
		     Z[vy] + pdinfo->t1, 
		     pdinfo->t2 - pdinfo->t1 + 1, 
		     &rho, &sd, &pval,
		     (opt)? &rx : NULL,
		     (opt)? &ry : NULL,
		     &m, RUN_FIRST)) {
	return E_ALLOC;
    }

    pprintf(prn, _("\nFor the variables '%s' and '%s'\n"), pdinfo->varname[vx],
	    pdinfo->varname[vy]);
    pprintf(prn, _("Spearman's rank correlation coefficient (rho) = %f\n"), rho);

    if (!na(pval)) {
	pprintf(prn, _("Under the null hypothesis of no correlation, rho "
		       "follows N(0, %f)\n"), sd);
	pprintf(prn, _("z-score = %f, with one-tailed p-value %f\n"), rho / sd,
		pval);
    } else if (m >= 7) {
	pval = spearman_signif(m, fabs(rho));
	if (pval < 1.0) {
	    pprintf(prn, _("significant at the %g%% level (one-tailed)\n"), 100.0 * pval);
	} else {
	    /* xgettext:no-c-format */
	    pputs(prn, _("not significant at the 10% level\n"));
	}
    } else {
	pputs(prn, _("Sample is too small to calculate a p-value based on "
		"the normal distribution\n"));
    }

    if (opt) { /* print raw and ranked data */
	int i = 0;

	pprintf(prn, "\n     %s ", _("Obs"));
	pprintf(prn, "%13s%13s%13s%13s\n\n", pdinfo->varname[vx], _("rank"),
	       pdinfo->varname[vy], _("rank"));
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    print_obs_marker(t, pdinfo, prn);
	    if (!(na(Z[vx][t])) && !(na(Z[vy][t]))) {
		gretl_printxn(Z[vx][t], 15, prn);
		gretl_printxn(rx[i], 15, prn);
		gretl_printxn(Z[vy][t], 15, prn);
		gretl_printxn(ry[i], 15, prn);
		i++;
	    }
	    pputc(prn, '\n');
	}
    }

    spearman_rho(NULL, NULL, 0, NULL, NULL, NULL, 
		 NULL, NULL, NULL, RUN_CLEANUP);

    return 0;
}

#undef LOCKE_DEBUG

static int randomize_doubles (const void *a, const void *b)
{
    return gretl_rand_int_max(8096) - 4097;
}

/* put sample into random order */

static double *locke_shuffle (const double *x, int *n, int code)
{
    static double *sx;
    static int m;
    int i;

    if (code == RUN_CLEANUP) {
	free(sx);
	sx = NULL;
	m = 0;
	return NULL;
    }	

    if (code == RUN_FIRST) {
	m = *n;
	
	sx = malloc(m * sizeof *sx);
	if (sx == NULL) {
	    return NULL;
	}

	/* check for missing and negative values as we go */

	m = 0;
	for (i=0; i<*n; i++) {
	    if (na(x[i])) {
		continue;
	    } else if (x[i] < 0.0) {
		m = 0;
		break;
	    } else {
		sx[m++] = x[i];
	    }
	}

	if (m == 0) {
	    free(sx);
	    return NULL;
	}

	m = 2 * m / 2;

	*n = m;
    } 

    qsort(sx, m, sizeof *sx, randomize_doubles);

    return sx;
}

#define NREPEAT 100

/* Charles Locke's nonparametric test for whether an empirical
   distribution is gamma.  See C. Locke, "A Test for the Composite
   Hypothesis that a Population has a Gamma Distribution,"
   Commun. Statis.-Theor. Meth. A5(4), 351-364 (1976).

   See also Shapiro and Chen, Journal of Quality Technology 33(1),
   Jan 2001.
*/

double lockes_test (const double *x, int t1, int t2)
{
    double rho, sd, z;
    double *sx, *u = NULL, *v = NULL;
    int i, j, t, m = t2 - t1 + 1;

    sx = locke_shuffle(x + t1, &m, RUN_FIRST);

    if (sx == NULL) {
	return NADBL;
    }

    m /= 2;
	
    u = malloc(m * sizeof *u);
    v = malloc(m * sizeof *v);

    if (u == NULL || v == NULL) {
	free(u);
	free(v);
	locke_shuffle(NULL, NULL, RUN_CLEANUP);
	return NADBL;
    }

    z = 0.0;

    /* repeat the shuffling of the series NREPEAT times, since the
       test statistic is sensitive to the ordering under the null
    */

    for (j=0; j<NREPEAT; j++) {
	t = 0;
	for (i=0; i<m; i++) {
	    u[i] = sx[t] + sx[t+1];
	    v[i] = sx[t] / sx[t+1];
	    if (sx[t+1] / sx[t] > v[i]) {
		v[i] = sx[t+1] / sx[t];
	    }
	    t += 2;
	}
	spearman_rho(u, v, m, &rho, &sd, NULL, NULL, NULL, NULL,
		     (j == 0)? RUN_FIRST : RUN_REPEAT);
	z += rho / sd;
#if LOCKE_DEBUG
	printf("z[%d] = %g/%g = %g\n", j, rho, sd, rho / sd);
#endif
	locke_shuffle(NULL, NULL, RUN_REPEAT);
    }   

    free(u);
    free(v);

    z /= (double) NREPEAT;

#if LOCKE_DEBUG
    fprintf(stderr, "Spearman's rho: z = %g\n", z);
#endif 

    locke_shuffle(NULL, NULL, RUN_CLEANUP);
    spearman_rho(NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL, NULL,
		 RUN_CLEANUP);

    return z;
}

/**
 * runs_test:
 * @varno: ID number of the variable to process.
 * @Z: data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 *
 * Performs, and prints the results of, the runs test for randomness
 * for the variable specified by @varno.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int runs_test (int varno, const double **Z, const DATAINFO *pdinfo, 
	       PRN *prn)
{
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2, n = pdinfo->n;
    int nn, runs = 1;
    double xx, *x, mean, sd;
    double z, pval;

    nn = t2 - t1 + 1;
    x = malloc(nn * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    nn = 0;
    for (t=t1; t<=t2; t++) {
	xx = Z[varno][t];
	if (na(xx)) {
	    continue;
	}
	x[nn++] = xx;
    } 

    if (nn <= 1) {
	pputs(prn, _("\nInsufficient data for runs test\n"));
	free(x);
	return 1;
    }

    for (t=1; t<nn; t++) {
	if ((x[t] > 0 && x[t-1] < 0) || (x[t] < 0 && x[t-1] > 0)) { 
	    runs++;
	}
    }

    mean = (1 + nn / 2.0);
    sd = sqrt((double) n - 1) / 2.0;
    z = fabs((runs - mean) / sd);
    pval = normal_pvalue_2(z);

    pprintf(prn, _("\nNumber of runs (R) in the variable '%s' = %d\n"), 
	    pdinfo->varname[varno], runs);
    pprintf(prn, _("Under the null hypothesis of randomness, R "
	    "follows N(%f, %f)\n"), mean, sd);
    pprintf(prn, _("z-score = %f, with two-tailed p-value %f\n"), z, pval);

    record_test_result(z, pval, "runs");
  
    free(x);

    return 0;
}



