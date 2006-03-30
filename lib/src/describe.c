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

/* describe.c - gretl descriptive statistics */

#include "libgretl.h"
#include "gretl_matrix.h"
#include "libset.h"
#include "compat.h"

#include <unistd.h>

#ifdef WIN32
# include <windows.h>
#endif

#if defined(ENABLE_NLS) && defined(USE_GTK2)
# include <glib.h>
#endif

static char gretl_tmp_str[MAXLEN];

/* return 1 if series x has any missing observations, 0 if it does
   not */

static int missvals (double *x, int n)
{
    int t, ret = 0;
    
    for (t=0; t<n; t++) {
	if (na(x[t])) {
	    ret = 1;
	    break;
	}
    }
	    
    return ret;
}

/* return the number of "good" (non-missing) observations in series x;
   also (optionally) write to x0 the first non-missing value */

static int good_obs (const double *x, int n, double *x0)
{
    int t, good = n;

    if (x0 != NULL) {
	*x0 = NADBL;
    }

    for (t=0; t<n; t++) {
	if (na(x[t])) {
	    good--;
	} else if (x0 != NULL && na(*x0)) {
	    *x0 = x[t];
	}
    }

    return good;
}

/**
 * gretl_minmax:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @min: pointer to receive minimum value.
 * @max: pointer to receive maximum value.
 *
 * Puts the minimum and maximum values of the series @x,
 * from obs @t1 to obs @t2, into the variables @min and @max.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_minmax (int t1, int t2, const double *x, 
		  double *min, double *max)
{
    int t;

    while (na(x[t1]) && t1 <= t2) {
	t1++;
    }

    if (t1 >= t2) {
        *min = *max = NADBL;
        return 1;
    }

    *min = x[t1];
    *max = x[t1];

    for (t=t1; t<=t2; t++) {
	if (!(na(x[t]))) {
	    *max = x[t] > *max ? x[t] : *max;
	    *min = x[t] < *min ? x[t] : *min;
	}
    }

    return 0;
}

/**
 * gretl_mean:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the arithmetic mean of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL 
 * in case there are no valid observations.
 */

double gretl_mean (int t1, int t2, const double *x)
{
    int n;
    register int t;
    double xbar, sum = 0.0;

    n = t2 - t1 + 1;
    if (n <= 0) {
	return NADBL;
    }

    for (t=t1; t<=t2; t++) {
	if (!(na(x[t]))) {
	    sum += x[t];
	} else {
	    n--;
	}
    }

    if (n == 0) {
	return NADBL;
    }

    xbar = sum / n;
    sum = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!(na(x[t]))) {
	    sum += (x[t] - xbar); 
	}
    }

    return xbar + sum / n;
}

/**
 * eval_ytest:
 * @y: reference value.
 * @op: operator.
 * @test: test value.
 *
 * Returns: 1 if the expression @y @yop @test (for example
 * "y = 2" or "y <= 45") evaluates as true, else 0.
 */

int eval_ytest (double y, GretlOp op, double test)
{
    int ret = 0;

    switch (op) {
    case OP_EQ:
	ret = (y == test);
	break;
    case OP_GT:
	ret = (y > test);
	break;
    case OP_LT:
	ret = (y < test);
	break;
    case OP_NEQ:
	ret = (y != test);
	break;
    case OP_GTE:
	ret = (y <= test);
	break;
    case OP_LTE:
	ret = (y <= test);
	break;
    default:
	break;
    }

    return ret;
}

/**
 * gretl_restricted_mean:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: criterion series.
 * @yop: criterion operator.
 * @yval: criterion value.
 *
 * Returns: the arithmetic mean of the series @x in the
 * range @t1 to @t2 (inclusive), but including only
 * observations where the criterion variable @y bears the
 * relationship @yop to the value @yval -- or #NADBL in case 
 * there are no observations that satisfy the restriction.
 */

double gretl_restricted_mean (int t1, int t2, const double *x,
			      const double *y, GretlOp yop, 
			      double yval)
{
    int n;
    register int t;
    double xbar, sum = 0.0;

    n = t2 - t1 + 1;
    if (n <= 0) {
	return NADBL;
    }

    for (t=t1; t<=t2; t++) {
	if (!na(x[t]) && eval_ytest(y[t], yop, yval)) {
	    sum += x[t];
	} else {
	    n--;
	}
    }

    if (n == 0) {
	return NADBL;
    }

    xbar = sum / n;
    sum = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t]) && eval_ytest(y[t], yop, yval)) {
	    sum += (x[t] - xbar); 
	}
    }

    return xbar + sum / n;
}

/**
 * gretl_median:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the median value of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL 
 * on failure.
 */

double gretl_median (int t1, int t2, const double *x)
{
    int m = t2 - t1 + 1;
    double *sx, med;
    int t, n, n2p;

    sx = malloc(m * sizeof *sx);

    if (sx == NULL) {
	return NADBL;
    }

    n = 0;
    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    sx[n++] = x[t];
	}
    }

    qsort(sx, n, sizeof *sx, gretl_compare_doubles); 

    n2p = (m = n / 2) + 1;
    med = (n % 2)? sx[n2p - 1] : 0.5 * (sx[m - 1] + sx[n2p - 1]);

    free(sx);

    return med;
}

/**
 * gretl_sst:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the sum of squared deviations from the mean for
 * the series @x from obs @t1 to obs @t2, skipping any missing 
 * values, or #NADBL on failure.
 */

double gretl_sst (int t1, int t2, const double *x)
{
    int t, n = t2 - t1 + 1;
    double sumsq, xx, xbar;

    if (n == 0) {
	return NADBL;
    }

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
	return NADBL;
    }

    sumsq = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    xx = x[t] - xbar;
	    sumsq += xx * xx;
	} 
    }

    return sumsq;
}

/**
 * gretl_variance:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the variance of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL 
 * on failure.
 */

double gretl_variance (int t1, int t2, const double *x)
{
    int t, n = t2 - t1 + 1;
    double sumsq, xx, xbar;

    if (n == 0) {
	return NADBL;
    }

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
	return NADBL;
    }

    sumsq = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    xx = x[t] - xbar;
	    sumsq += xx * xx;
	} else {
	    n--;
	}
    }

    sumsq = (n > 1)? sumsq / (n - 1) : 0.0;

    return (sumsq >= 0)? sumsq : NADBL;
}

/**
 * gretl_restricted_variance:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: criterion series.
 * @yop: criterion operator. 
 * @yval: criterion value.
 *
 * Returns: the variance of the series @x from obs
 * @t1 to obs @t2, skipping any missing values and
 * observations where the series @y does not bear the
 * relationship @yop to the value @yval, or #NADBL on 
 * failure.
 */

double gretl_restricted_variance (int t1, int t2, const double *x,
				  const double *y, GretlOp yop, 
				  double yval)
{
    int t, n = t2 - t1 + 1;
    double sumsq, xx, xbar;

    if (n == 0) {
	return NADBL;
    }

    xbar = gretl_restricted_mean(t1, t2, x, y, yop, yval);
    if (na(xbar)) {
	return NADBL;
    }

    sumsq = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t]) && eval_ytest(y[t], yop, yval)) {
	    xx = x[t] - xbar;
	    sumsq += xx * xx;
	} else {
	    n--;
	}
    }

    sumsq = (n > 1)? sumsq / (n - 1) : 0.0;

    return (sumsq >= 0)? sumsq : NADBL;
}

/**
 * gretl_stddev:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the standard deviation of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL 
 * on failure.
 */

double gretl_stddev (int t1, int t2, const double *x)
{
    double xx = gretl_variance(t1, t2, x);

    return (na(xx))? xx : sqrt(xx);
}

/**
 * gretl_restricted_stddev:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: criterion series.
 * @yop: criterion operator.
 * @yval: criterion value.
 *
 * Returns: the standard deviation of the series @x from obs
 * @t1 to obs @t2, skipping any missing values and observations
 * where the series @y does not bear the relationship @yop to
 * the value @yval, or #NADBL on failure.
 */

double gretl_restricted_stddev (int t1, int t2, const double *x,
				const double *y, GretlOp yop, 
				double yval)
{
    double xx = gretl_restricted_variance(t1, t2, x, y, yop, yval);

    return (na(xx))? xx : sqrt(xx);
}

/**
* gretl_long_run_variance:
* @t1: starting observation.
* @t2: ending observation.
* @x: data series.
* @m: bandwidth.
*
* Returns: the long-run variance of the series @x from obs
* @t1 to obs @t2, using Bartlett kernel weights, or #NADBL 
* on failure (which includes encountering missing values). 
*/

double gretl_long_run_variance (int t1, int t2, const double *x, int m)
{
    double zt, wt, xbar, s2 = 0.0;
    double *autocov;
    int i, t, n;

    if (array_adjust_t1t2(x, &t1, &t2)) {
	return NADBL;
    }

    n = t2 - t1 + 1;

    if (n < 2) {
	return NADBL;
    }

    xbar = gretl_mean(t1, t2, x);

    autocov = malloc(m * sizeof *autocov);
    if (autocov == NULL) {
	return NADBL;
    }
  
    for (i=0; i<m; i++) {
	autocov[i] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
	zt = x[t] - xbar;
	s2 += zt * zt;
	for (i=0; i<m; i++) {
	    int s = i + 1;

	    if (t - s >= t1) {
		autocov[i] += zt * (x[t - s] - xbar);
	    }
	}
    }

    for (i=0; i<m; i++) {
	wt = 1.0 - ((double) (i + 1)) / (m + 1.0);
	s2 += 2.0 * wt * autocov[i];
    }

    s2 /= n;

    free(autocov);

    return s2;
}

/**
 * gretl_covar:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: data series.
 *
 * Returns: the covariance of the series @x and @y from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL 
 * on failure.
 */

double gretl_covar (int t1, int t2, const double *x, const double *y)
{
    int t, nn, n = t2 - t1 + 1;
    double sx, sy, sxy, xt, yt, xbar, ybar;

    if (n == 0) {
	return NADBL;
    }

    nn = n;
    sx = sy = 0.0;

    for (t=t1; t<=t2; t++) {
        xt = x[t];
        yt = y[t];
        if (na(xt) || na(yt)) {
            nn--;
            continue;
        }
        sx += xt;
        sy += yt;
    }

    if (nn == 0) {
	return NADBL;
    }

    xbar = sx / nn;
    ybar = sy / nn;
    sxy = 0.0;

    for (t=t1; t<=t2; t++) {
        xt = x[t];
        yt = y[t];
        if (na(xt) || na(yt)) {
	    continue;
	}
        sx = xt - xbar;
        sy = yt - ybar;
        sxy = sxy + (sx * sy);
    }

    return sxy / (nn - 1);
}

/**
 * gretl_corr:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: data series.
 * @missing: location to receive information on the number
 * of missing observations that were skipped, or %NULL.
 *
 * Returns: the correlation coefficient for the series @x and @y 
 * from obs @t1 to obs @t2, skipping any missing values, or #NADBL 
 * on failure.
 */

double gretl_corr (int t1, int t2, const double *x, const double *y,
		   int *missing)
{
    int t, nn, n = t2 - t1 + 1;
    double sx, sy, sxx, syy, sxy, den, xbar, ybar;
    double cval = 0.0;

    if (n == 0) {
	return NADBL;
    }

    if (gretl_isconst(t1, t2, x) || gretl_isconst(t1, t2, y)) {
	return NADBL;
    }

    nn = n;
    sx = sy = 0.0;
    for (t=t1; t<=t2; t++) {
        if (na(x[t]) || na(y[t])) {
            nn--;
            continue;
        }
        sx += x[t];
        sy += y[t];
    }

    if (nn == 0) {
	return NADBL;
    }

    xbar = sx / nn;
    ybar = sy / nn;
    sxx = syy = sxy = 0.0;

    for (t=t1; t<=t2; t++) {
        if (na(x[t]) || na(y[t])) {
	    continue;
	}
        sx = x[t] - xbar;
        sy = y[t] - ybar;
	sxx += sx * sx;
	syy += sy * sy;
	sxy += sx * sy;
    }

    if (sxy != 0.0) {
        den = sxx * syy;
        if (den > 0.0) {
	    cval = sxy / sqrt(den);
        } else {
	    cval = NADBL;
	}
    }

    if (missing != NULL) {
	*missing = n - nn;
    }

    return cval;
}

/**
 * gretl_corr_rsq:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: data series.
 *
 * Returns: the square of the correlation coefficient for the series 
 * @x and @y from obs @t1 to obs @t2, skipping any missing values, 
 * or #NADBL on failure.  Used as alternative value for R^2 in a
 * regression without an intercept.
 */

double gretl_corr_rsq (int t1, int t2, const double *x, const double *y)
{
    double r = gretl_corr(t1, t2, x, y, NULL);

    if (na(r)) {
	return NADBL;
    } else {
	return r * r;
    }
}

/* we're supposing a variance smaller than this is just noise */
#define TINYVAR 1.0e-16 

/**
 * gretl_moments:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @xbar: pointer to receive mean.
 * @sd: pointer to receive standard deviation.
 * @skew: pointer to receive skewness.
 * @kurt: pointer to receive excess kurtosis.
 * @k: degrees of freedom loss (generally 1).
 *
 * Calculates sample moments for series @x from obs @t1 to obs
 * @t2.  
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_moments (int t1, int t2, const double *x, 
		   double *xbar, double *std, 
		   double *skew, double *kurt, int k)
{
    int t, n;
    double dev, var;
    double s, s2, s3, s4;
    int allstats = 1;

    if (skew == NULL && kurt == NULL) {
	allstats = 0;
    }

    while (na(x[t1]) && t1 <= t2) {
	t1++;
    }

    if (gretl_isconst(t1, t2, x)) {
	*xbar = x[t1];
	*std = 0.0;
	if (allstats) {
	    *skew = *kurt = NADBL;
	}
	return 1;
    }

    s = 0.0;
    n = 0;
    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    s += x[t];
	    n++;
	}
    }

    if (n == 0) {
	*xbar = *std = NADBL;
	if (allstats) {
	    *skew = *kurt = 0.0;
	}	
	return 1;
    }

    *xbar = s / n;
    var = 0.0;
    if (allstats) {
	*skew = *kurt = 0.0;
    }

    s2 = s3 = s4 = 0.0;
    for (t=t1; t<=t2; t++) {
	if (na(x[t])) {
	    continue;
	}
	dev = x[t] - *xbar;
	s2 += dev * dev;
	if (allstats) {
	    s3 += pow(dev, 3);
	    s4 += pow(dev, 4);
	}
    }

    var = s2 / (n - k);

    if (var < 0.0) {
	*std = NADBL;
	if (allstats) {
	    *skew = *kurt = NADBL;
	}
	return 1;
    }

    if (var > TINYVAR) {
	*std = sqrt(var);
    } else {
	*std = 0.0;
    }

    if (allstats) {
	if (var > TINYVAR) {
	    /* if variance is effectively zero, these should be undef'd */
	    *skew = (s3 / n) / pow(s2 / n, 1.5);
	    *kurt = ((s4 / n) / pow(s2 / n, 2)) - 3.0; /* excess kurtosis */
	} else {
	    *skew = *kurt = NADBL;
	}
    }

    return 0;
}

/**
 * free_freq:
 * @freq: gretl frequency distribution struct
 *
 * Frees all malloced elements of the struct.
 *
 */

void free_freq (FreqDist *freq)
{
    if (freq == NULL) {
	return;
    }

    free(freq->midpt);
    free(freq->endpt);
    free(freq->f);

    free(freq);
}

static FreqDist *freq_new (void)
{
    FreqDist *freq;

    freq = malloc(sizeof *freq);
    if (freq == NULL) return NULL;

    freq->midpt = NULL;
    freq->endpt = NULL;
    freq->f = NULL;

    freq->dist = 0;
    freq->test = NADBL;

    return freq;
}

/**
 * dh_root_b1_to_z1:
 * @rb1: square root b1, skewness.
 * @n: number of observations.
 *
 * Performs the transformation from skewness, root b1, to the 
 * normal score, z1, as set out in Doornik and Hansen, "An Omnibus 
 * Test for Normality", 1994.  The transformation is originally 
 * due to D'Agostino (Biometrika, 1970).
 *
 * Returns: the z1 value.
 */

double dh_root_b1_to_z1 (double rb1, double n)
{
    double b, w2, d, y, z1;

    b = 3.0 * (n*n + 27*n - 70) * (n+1) * (n+3) /
	((n-2) * (n+5) * (n+7) * (n+9));

    w2 = -1.0 + sqrt(2 * (b-1));
    
    d = 1.0 / sqrt(log(sqrt(w2)));

    y = rb1 * sqrt(((w2-1.0)/2.0) * ((n+1.0)*(n+3.0))/(6.0*(n-2)));

    z1 = d * log(y + sqrt(y*y + 1));

    return z1;
}

/**
 * dh_b2_to_z2:
 * @b1: skewness.
 * @b2: kurtosis.
 * @n: number of observations.
 *
 * Performs the transformation from kurtosis, b2, to the 
 * normal score, z2, as set out in Doornik and Hansen, "An Omnibus 
 * Test for Normality", 1994.
 *
 * Returns: the z2 value.
 */

double dh_b2_to_z2 (double b1, double b2, double n)
{
    double d, a, c, k, alpha, chi, z2;
    double n2 = n * n;

    d = (n-3) * (n+1) * (n2 + 15*n - 4.0);

    a = ((n-2) * (n+5) * (n+7) * (n2 + 27*n - 70.0)) / (6.0 * d);

    c = ((n-7) * (n+5) * (n+7) * (n2 + 2*n - 5.0)) / (6.0 * d);

    k = ((n+5) * (n+7) * (n2*n + 37*n2 + 11*n - 313.0)) / (12.0 * d);

    alpha = a + b1 * c;

    chi = (b2 - 1.0 - b1) * 2.0 * k;

    z2 = (pow(chi/(2*alpha), 1.0/3.0) - 1.0 + (1.0 / (9.0*alpha))) *
	sqrt(9.0*alpha);

    return z2;
}

/**
 * doornik_chisq:
 * @skew: skewness.
 * @xkurt: excess kurtosis.
 * @n: number of observations.
 *
 * Calculates the Chi-square test for normality as set out by
 * Doornik and Hansen, "An Omnibus Test for Normality", 1994.
 * This is a modified version of the test proposed by Bowman and
 * Shenton (Biometrika, 1975).
 *
 * Returns: the Chi-square value, which has 2 degrees of freedom.
 */

double doornik_chisq (double skew, double xkurt, int n)
{
    double rb1, b1, b2, z1, z2;

    rb1 = skew;
    b1 = skew * skew;
    b2 = xkurt + 3.0; /* Note: convert from "excess" to regular */

    z1 = dh_root_b1_to_z1(rb1, (double) n);
    z2 = dh_b2_to_z2(b1, b2, (double) n);

    return z1*z1 + z2*z2;
}

static int
get_moments (const gretl_matrix *M, int row, double *skew, double *kurt)
{
    int j, n = gretl_matrix_cols(M);
    double xi, xbar, dev, var;
    double s = 0.0;
    double s2 = 0.0;
    double s3 = 0.0;
    double s4 = 0.0;
    int err = 0;
    
    for (j=0; j<n; j++) {
	s += gretl_matrix_get(M, row, j);
    }

    xbar = s / n;

    for (j=0; j<n; j++) {
	xi = gretl_matrix_get(M, row, j);
	dev = xi - xbar;
	s2 += dev * dev;
	s3 += pow(dev, 3);
	s4 += pow(dev, 4);
    }

    var = s2 / n;

    if (var > 0.0) {
	/* if variance is effectively zero, these should be undef'd */
	*skew = (s3 / n) / pow(s2 / n, 1.5);
	*kurt = ((s4 / n) / pow(s2 / n, 2));
    } else {
	*skew = *kurt = NADBL;
	err = 1;
    }

    return err;
}

int 
gretl_system_normality_test (const gretl_matrix *E, const gretl_matrix *Sigma, 
			     PRN *prn)
{
    gretl_matrix *S = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *C = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *tmp = NULL;

    /* convenience pointers: do not free! */
    gretl_matrix *H;
    gretl_vector *Z1;
    gretl_vector *Z2;

    double *evals = NULL;
    double x, skew, kurt;
    double X2 = NADBL;
    int n, p;
    int i, j;
    int err = 0;

    if (E == NULL || Sigma == NULL) {
	err = 1;
	goto bailout;
    }

    p = gretl_matrix_cols(E);
    n = gretl_matrix_rows(E);

    S = gretl_matrix_copy(Sigma);
    V = gretl_vector_alloc(p);
    C = gretl_matrix_alloc(p, p);
    X = gretl_matrix_copy_transpose(E);
    R = gretl_matrix_alloc(p, n);
    tmp = gretl_matrix_alloc(p, p);

    if (S == NULL || V == NULL || C == NULL || X == NULL || 
	R == NULL || tmp == NULL) {
	err = 1;
	goto bailout;
    }

    for (i=0; i<p; i++) {
	x = gretl_matrix_get(S, i, i);
	gretl_vector_set(V, i, 1.0 / sqrt(x));
    }

    err = gretl_matrix_diagonal_sandwich(V, S, C);
    if (err) {
	goto bailout;
    }

    gretl_matrix_print_to_prn(C, "\nResidual correlation matrix, C", prn);

    evals = gretl_symmetric_matrix_eigenvals(C, 1, &err);
    if (err) {
	goto bailout;
    }

    pputs(prn, "Eigenvalues of the correlation matrix:\n\n");
    for (i=0; i<p; i++) {
	pprintf(prn, " %10g\n", evals[i]);
    }
    pputc(prn, '\n');

    /* C should now contain eigenvectors of the original C:
       relabel as 'H' for perspicuity */
    H = C;
#if 0
    gretl_matrix_print_to_prn(H, "Eigenvectors, H", prn);
#endif
    gretl_matrix_copy_values(tmp, H);

    /* make "tmp" into $H \Lambda^{-1/2}$ */
    for (i=0; i<p; i++) {
	for (j=0; j<p; j++) {
	    x = gretl_matrix_get(tmp, i, j);
	    x *= 1.0 / sqrt(evals[j]);
	    gretl_matrix_set(tmp, i, j, x);
	}
    }

    /* make S into $H \Lambda^{-1/2} H'$ */
    gretl_matrix_multiply_mod(tmp, GRETL_MOD_NONE,
			      H, GRETL_MOD_TRANSPOSE,
			      S);

#if 1
    gretl_matrix_demean_by_row(X);
#endif

    /* compute VX', in X (don't need to subtract means, because these
       are OLS residuals)
    */
    for (i=0; i<p; i++) {
	for (j=0; j<n; j++) {
	    x = gretl_matrix_get(X, i, j);
	    x *= gretl_vector_get(V, i);
	    gretl_matrix_set(X, i, j, x);
	}
    }

    /* finally, compute $R' = H  \Lambda^{-1/2} H' V X'$ 
       Doornik and Hansen, 1994, section 3 */
    gretl_matrix_multiply(S, X, R);

    /* Z_1 and Z_2 are p-vectors: use existing storage */
    Z1 = V;
    Z2 = gretl_matrix_reuse(tmp, p, 1);

    for (i=0; i<p && !err; i++) {
	get_moments(R, i, &skew, &kurt);
	if (na(skew) || na(kurt)) {
	    err = 1;
	} else {
	    double z1i = dh_root_b1_to_z1(skew, n);
	    double z2i = dh_b2_to_z2(skew * skew, kurt, n);

	    gretl_vector_set(Z1, i, z1i);
	    gretl_vector_set(Z2, i, z2i);
	}
    }
	
    if (!err) {
	X2 = gretl_vector_dot_product(Z1, Z1, &err);
	X2 += gretl_vector_dot_product(Z2, Z2, &err);
    }

    if (na(X2)) {
	pputs(prn, "Calculation of test statistic failed\n");
    } else {
	pputs(prn, "Test for multivariate normality of residuals\n");
	pprintf(prn, "Doornik-Hansen Chi-square(%d) = %g, ", 2 * p, X2);
	pprintf(prn, "with p-value = %g\n", chisq(X2, 2 * p));
    }

 bailout:

    gretl_matrix_free(S);
    gretl_matrix_free(V);
    gretl_matrix_free(C);
    gretl_matrix_free(X);
    gretl_matrix_free(R);
    gretl_matrix_free(tmp);

    free(evals);

    return err;
}

static int freq_add_arrays (FreqDist *freq, int nbins)
{
    int err = 0;
	
    freq->endpt = malloc((nbins + 1) * sizeof *freq->endpt);
    freq->midpt = malloc(nbins * sizeof *freq->midpt);
    freq->f = malloc(nbins * sizeof *freq->f);
	
    if (freq->endpt == NULL || freq->midpt == NULL || freq->f == NULL) {
	err = E_ALLOC;
	strcpy(gretl_errmsg, _("Out of memory for frequency distribution"));
    } else {
	freq->numbins = nbins;
    }

    return err;
}

static int count_distinct_vals (int *x, int n)
{
    int i, c = 1;

    for (i=1; i<n; i++) {
	if (x[i] != x[i-1]) c++;
    }

    return c;
}

static int compare_ints (const void *a, const void *b)
{
    const int *ia = (const int *) a;
    const int *ib = (const int *) b;
     
    return *ia - *ib;
}

/**
 * get_freq:
 * @varno: ID number of variable to process.
 * @Z: data array.
 * @pdinfo: information on the data set.
 * @params: degrees of freedom loss (generally = 1 unless we're dealing
 * with the residual from a regression)
 * @opt: if & OPT_O, compare with gamma distribution, not normal
 *
 * Calculates the frequency distribution for the specified variable.
 *
 * Returns: pointer to struct containing the distribution.
 */

FreqDist *get_freq (int varno, const double **Z, const DATAINFO *pdinfo, 
		    int params, gretlopt opt)
{
    FreqDist *freq;
    const double *x = Z[varno];
    double xx, xmin, xmax;
    double skew, kurt;
    double range, binwidth;
    int nbins;
    int t, k, n;

    freq = freq_new();
    if (freq == NULL) {
	return NULL;
    }

    gretl_errno = 0;
    gretl_errmsg[0] = '\0';

    n = pdinfo->t2 - pdinfo->t1 + 1;
    n = good_obs(x + pdinfo->t1, n, NULL);

    if (n < 8) {
	gretl_errno = E_DATA;
	sprintf(gretl_errmsg, _("Insufficient data to build frequency "
		"distribution for variable %s"), pdinfo->varname[varno]);
	free_freq(freq);
	return NULL;
    }

    freq->t1 = pdinfo->t1; 
    freq->t2 = pdinfo->t2;
    freq->n = n;

    strcpy(freq->varname, pdinfo->varname[varno]);

    if (gretl_isconst(pdinfo->t1, pdinfo->t2, x)) {
	gretl_errno = 1;
	sprintf(gretl_errmsg, _("%s is a constant"), freq->varname);
	free_freq(freq);
	return NULL;
    } 

    freq->discrete = gretl_isdiscrete(pdinfo->t1, pdinfo->t2, x);

    gretl_moments(pdinfo->t1, pdinfo->t2, x, 
		  &freq->xbar, &freq->sdx, 
		  &skew, &kurt, params);

    gretl_minmax(pdinfo->t1, pdinfo->t2, x, &xmin, &xmax);
    range = xmax - xmin;

    /* calculate test stat for distribution */
    if (freq->n > 7) {
	if (opt & OPT_O) {
	    freq->test = lockes_test(x, pdinfo->t1, pdinfo->t2);
	    freq->dist = DIST_GAMMA;
	} else {
	    freq->test = doornik_chisq(skew, kurt, freq->n); 
	    freq->dist = DIST_NORMAL;
	}
    } else {
	freq->test = NADBL;
	freq->dist = 0;
    }

    /* if the histogram is not wanted, we're done */
    if (opt & OPT_Q) {
	freq->numbins = 0;
	return freq;
    }
    
    if (freq->discrete) {
	int *ifreq = NULL, *ivals = NULL;
	int *sorted = NULL;
	int i, last, err = 0;

	sorted = malloc(n * sizeof *sorted);
	if (sorted == NULL) {
	    free_freq(freq);
	    return NULL;
	}
	
	n = 0;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (!na(x[t])) {
		sorted[n++] = (int) x[t];
	    }
	}
	
	qsort(sorted, n, sizeof *sorted, compare_ints); 
	nbins = count_distinct_vals(sorted, n);

	ifreq = malloc(nbins * sizeof *ifreq);
	ivals = malloc(nbins * sizeof *ivals);
	if (ifreq == NULL || ivals == NULL) {
	    err = E_ALLOC;
	    goto discrete_end;
	}

	ivals[0] = last = sorted[0];
	ifreq[0] = i = 1;

	for (t=1; t<n; t++) {
	    if (sorted[t] != last) {
		last = sorted[t];
		ifreq[i] = 1;
		ivals[i++] = last;
	    } else {
		ifreq[i-1] += 1;
	    }
	}

	if (freq_add_arrays(freq, nbins)) {
	    err = E_ALLOC;
	} else {
	    for (k=0; k<nbins; k++) {
		freq->endpt[k] = freq->midpt[k] = ivals[k];
		freq->f[k] = ifreq[k];
	    }
	    freq->endpt[nbins] = xmax;
	}

    discrete_end:

	free(sorted);
	free(ivals);
	free(ifreq);
	
	if (err) {
	    gretl_errno = E_ALLOC;
	    strcpy(gretl_errmsg, _("Out of memory for frequency distribution"));
	    free_freq(freq);
	    freq = NULL;
	}
    } else {
	/* not discrete (integer) data */
	if (n < 16) {
	    nbins = 5; 
	} else if (n < 50) {
	    nbins = 7;
	} else if (n > 850) {
	    nbins = 29;
	} else {
	    nbins = (int) sqrt((double) n);
	    if (nbins % 2 == 0) {
		nbins++;
	    }
	}
	
	if (freq_add_arrays(freq, nbins)) {
	    gretl_errno = E_ALLOC;
	    strcpy(gretl_errmsg, _("Out of memory for frequency distribution"));
	    return freq;
	}

	binwidth = range / (nbins - 1);
	
	freq->endpt[0] = xmin - .5 * binwidth;
	
	if (xmin > 0.0 && freq->endpt[0] < 0.0) {
	    double rshift;
	    
	    freq->endpt[0] = 0.0;
	    rshift = 1.0 - xmin / binwidth;
	    freq->endpt[freq->numbins] = xmax + rshift * binwidth;
	} else {
	    freq->endpt[freq->numbins] = xmax + .5 * binwidth;
	}
	
	for (k=0; k<freq->numbins; k++) {
	    freq->f[k] = 0;
	    if (k > 0) {
		freq->endpt[k] = freq->endpt[k-1] + binwidth;
	    }
	    freq->midpt[k] = freq->endpt[k] + .5 * binwidth;
	}
	
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    xx = x[t];
	    if (na(xx)) {
		continue;
	    }
	    if (xx < freq->endpt[1]) {
		freq->f[0] += 1;
		continue;
	    }
	    if (xx >= freq->endpt[freq->numbins]) {
		freq->f[freq->numbins-1] += 1;
		continue;
	    }
	    for (k=1; k<freq->numbins; k++) {
		if (freq->endpt[k] <= xx && xx < freq->endpt[k+1]) {
		    freq->f[k] += 1;
		    break;
		}
	    }
	}
    }
	
    return freq;
}

int freqdist (int varno, const double **Z, const DATAINFO *pdinfo,
	      int graph, gretlopt opt, PRN *prn)
{
    FreqDist *freq;

    freq = get_freq(varno, Z, pdinfo, 1, opt); 

    if (freq == NULL) {
	return E_ALLOC;
    }

    print_freq(freq, prn); 

    if (graph && !(opt & OPT_Q)) {
	if (plot_freq(freq, (opt)? DIST_GAMMA : DIST_NORMAL)) {
	    pputs(prn, _("gnuplot command failed\n"));
	}
    }

    free_freq(freq);

    return 0;
}

/**
 * free_xtab:
 * @xtab: gretl crosstab struct
 *
 * Frees all malloced elements of the struct, and then
 * the pointer itself.
 */

void free_xtab (Xtab *tab)
{
    int i;

    free(tab->rtotal);
    free(tab->ctotal);
    free(tab->rval);
    free(tab->cval);

    for (i=0; i<tab->rows; i++) {
	free(tab->f[i]);
    }
    free(tab->f);

    free(tab);
}

static Xtab *xtab_new (int n, int t1, int t2)
{
    Xtab *tab = malloc(sizeof *tab);

    if (tab == NULL) return NULL;

    tab->rtotal = NULL;
    tab->ctotal = NULL;
    tab->rval = NULL;
    tab->cval = NULL;
    tab->f = NULL;

    tab->n = n;
    tab->t1 = t1;
    tab->t2 = t2;

    return tab;
}

static int xtab_allocate_arrays (Xtab *tab, int rows, int cols)
{
    int i, j;

    tab->rows = rows;
    tab->cols = cols;

    tab->rval = malloc(rows * sizeof tab->rval);
    tab->rtotal = malloc(rows * sizeof tab->rtotal);

    tab->cval = malloc(cols * sizeof tab->cval);
    tab->ctotal = malloc(cols * sizeof tab->ctotal);

    tab->f = malloc(rows * sizeof tab->f);

    if (tab->rval == NULL || tab->rtotal == NULL ||
	tab->cval == NULL || tab->ctotal == NULL ||
	tab->f == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<rows; i++) {
	tab->f[i] = NULL;
    }

    for (i=0; i<rows; i++) {
	tab->f[i] = malloc(cols * sizeof tab->f[i]);
	if (tab->f[i] == NULL) {
	    return E_ALLOC;
	} else {
	    for (j=0; j<cols; j++) {
		tab->f[i][j] = 0;
	    }
	}
    }

    return 0;
}

int compare_xtab_rows (const void *a, const void *b) 
{
    const int **da = (const int **) a;
    const int **db = (const int **) b;
    int ret;
    
    ret = da[0][0] - db[0][0];

    if (ret == 0) {
	ret = da[0][1] - db[0][1];
    }
    
    return ret;
}

/**
  crosstab struct creation function
*/

Xtab *get_xtab (int rvarno, int cvarno, const double **Z, 
		const DATAINFO *pdinfo)
{
    Xtab *tab = NULL;
    int **X = NULL;
    int rx, cx;
    int rows, cols;

    FreqDist *rowfreq, *colfreq;

    int t, i, nr, nc;
    int t1 = pdinfo->t1;
    int t2 = pdinfo->t2;
    int n = 0;
    int err = 0;

    /* check if both variables are discrete */

    if (!gretl_isdiscrete(t1, t2, Z[rvarno]) ||
	!gretl_isdiscrete(t1, t2, Z[cvarno])) {
	fprintf(stderr, "One of the variables is not discrete!\n");
	return NULL;
    }

    /* count non-missing values */

    for (t=t1; t<=t2; t++) {
	if (!(na(Z[rvarno][t]) || na(Z[cvarno][t]))) {
	    n++;
	}
    }

    if (n == 0) {
	fprintf(stderr, "All values invalid!\n");
	return NULL;
    }

    /* Put in some info we already know */

    tab = xtab_new(n, t1, t2);
    if (tab == NULL) {
	return NULL;
    }

    tab->missing = (t2 - t1 + 1) - n;
    strcpy(tab->rvarname, pdinfo->varname[rvarno]);
    strcpy(tab->cvarname, pdinfo->varname[cvarno]);

    /* 
       start allocating stuff; we use temporary FreqDists for rows 
       and columns to retrieve values with non-zero frequencies
       and get dimensions for the cross table
     */

    rowfreq = get_freq(rvarno, Z, pdinfo, 0, OPT_NONE); 
    colfreq = get_freq(cvarno, Z, pdinfo, 0, OPT_NONE); 
    X = malloc(n * sizeof *X);

    if (rowfreq == NULL || colfreq == NULL || X == NULL) {
	free_freq(rowfreq);
	free_freq(colfreq);
	free(X);
	return NULL;
    }

    rows = rowfreq->numbins;
    cols = colfreq->numbins;

    if (xtab_allocate_arrays(tab, rows, cols)) {
	free_xtab(tab);
	return NULL;
    }

    for (i=0; i<rows; i++) {
	tab->rval[i] = rowfreq->midpt[i];
	tab->rtotal[i] = 0;
    }

    for (i=0; i<cols; i++) {
	tab->cval[i] = colfreq->midpt[i];
	tab->ctotal[i] = 0;
    }

    /* done with FreqDists */
    free_freq(rowfreq);
    free_freq(colfreq);

    /* 
       build the matrix X holding the values to be sorted
    */

    i = 0;
    for (t=t1; t<=t2 && i<n; t++) {
	X[i] = malloc(2 * sizeof **X);
	if (X[i] == NULL) {
	    err = E_ALLOC;
	} else {
	    if (!(na(Z[rvarno][t]) || na(Z[cvarno][t]))) { 
		X[i][0] = (int) Z[rvarno][t];
		X[i][1] = (int) Z[cvarno][t];
 		i++;
	    }
	}
    }

    qsort(X, n, sizeof *X, compare_xtab_rows);

    rx = cx = 0;
    nr = tab->rval[0];
    nc = tab->cval[0];

    /* compute frequencies by going through sorted X */

    for (i=0; i<n; i++) {
	while (X[i][0] > nr) { /* skip row */
	    nr = tab->rval[++rx];
	    cx = 0;
	    nc = tab->cval[0];
	}
	while (X[i][1] > nc) { /* skip column */
	    nc = tab->cval[++cx];
	}
#if XTAB_DEBUG
	fprintf(stderr,"%d: (%d,%d) [%d,%d] %d,%d\n",
		i, rx, cx, nr, nc, X[i][0], X[i][1]);
#endif
	tab->f[rx][cx] += 1;
	tab->rtotal[rx] += 1;
	tab->ctotal[cx] += 1;
    }

    /* we're done; free stuff and quit */

    for (i=0; i<n; i++) {
	free(X[i]);
    }
    free(X);

    return tab;
}

int crosstab (const int *list, const double **Z, 
	      const DATAINFO *pdinfo, gretlopt opt, 
	      PRN *prn)
{
    Xtab *tab;
    int rowvar;
    int colvar;

    if (list[0] != 2) {
	return E_PARSE;
    }

    rowvar = list[1];
    colvar = list[2];

    if (rowvar >= pdinfo->v || colvar >= pdinfo->v ||
	!pdinfo->vector[rowvar] || !pdinfo->vector[colvar]) {
	/* paranoia */
	return E_DATA;
    }

    tab = get_xtab(rowvar, colvar, Z, pdinfo); 

    if (tab == NULL) {
	return E_ALLOC;
    }

    print_xtab(tab, opt, prn); 

    free_xtab(tab);

    return 0;
}

/* ---------- end crosstab stuff -------------------- */

int model_error_dist (const MODEL *pmod, double ***pZ,
		      DATAINFO *pdinfo, PRN *prn)
{
    FreqDist *freq = NULL;
    int err = 0;

    if (genr_fit_resid(pmod, pZ, pdinfo, GENR_RESID, 1)) {
	pputs(prn, _("Out of memory attempting to add variable\n"));
	err = E_ALLOC;
    }

    if (!err) {
	freq = get_freq(pdinfo->v - 1, (const double **) *pZ, pdinfo, 
			pmod->ncoeff, OPT_NONE);
	if (freq == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	print_freq(freq, prn); 
	free_freq(freq);
    }

    dataset_drop_last_variables(1, pZ, pdinfo);

    return err;

}

/* PACF via Durbin-Levinson algorithm */

static int get_pacf (double *pacf, const double *acf, int m)
{
    int i, j;
    gretl_matrix *phi;
    double x, num, den;

    phi = gretl_matrix_alloc(m, m);
    if (phi == NULL) return 1;

    pacf[0] = acf[0];
    gretl_matrix_set(phi, 0, 0, acf[0]);

    for (i=1; i<m; i++) {
	num = acf[i];
	for (j=0; j<i; j++) {
	    num -= gretl_matrix_get(phi, i-1, j) * acf[i-j-1];
	}
	den = 1.0;
	for (j=0; j<i; j++) {
	    den -= gretl_matrix_get(phi, i-1, j) * acf[j];
	}
	pacf[i] = num / den;
	gretl_matrix_set(phi, i, i, pacf[i]);
	for (j=0; j<i; j++) {
	    x = gretl_matrix_get(phi, i-1, j);
	    x -= pacf[i] * gretl_matrix_get(phi, i-1, i-j-1);
	    gretl_matrix_set(phi, i, j, x);
	}
    }
	
    gretl_matrix_free(phi);

    return 0;
}

int auto_acf_order (int pd, int nobs)
{
    int m;

    switch (pd) {
    case 4: 
	m = (nobs <= 20)? nobs - 5 : 14; 
	break;
    case 12: 
    case 52: 
	m = (nobs <= 40)? nobs - 13 : 28;
	break;
    case 24: 
	m = (nobs <= 100)? nobs - 25 : 96;
	break;
    default:  
	m = (nobs <= 18)? nobs - 5 : 14;
	break;
    }

    if (m > nobs / 5) {
	/* restrict to 20 percent of data (Tadeusz) */
	m = nobs / 5;
    }

    return m;
}

static double gretl_acf (int n, int k, const double *y)
{
    int t;
    double z, ybar, num = 0.0, den = 0.0;

    if (n == 0 || gretl_isconst(0, n - 1, y)) { 
	return NADBL;
    }

    ybar = gretl_mean(0, n-1, y);

    for (t=0; t<n; t++) {
	z = y[t] - ybar;
	den += z * z;
	if (t >= k) {
	    num += z * (y[t-k] - ybar);
	}
    }

    return num / den;
}

static char *corrgm_crit_string (void)
{
    if (get_local_decpoint() == ',') {
	return "1,96/T^0,5";
    } else {
	return "1.96/T^0.5";
    }
}

/**
 * corrgram:
 * @varno: ID number of variable to process.
 * @order: integer order for autocorrelation function.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * @opt: if includes %OPT_A, use ASCII graphics; if includes
 * %OPT_R, variable in question is a model residual generated
 * "on the fly".
 *
 * Computes the autocorrelation function and plots the correlogram for
 * the variable specified by @varno.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int corrgram (int varno, int order, double ***pZ, 
	      DATAINFO *pdinfo, PRN *prn, gretlopt opt)
{
    double *acf, box, pm;
    double *pacf = NULL;
    int k, l, acf_m, pacf_m, nobs; 
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int list[2];
    FILE *fq = NULL;
    int err = 0, pacf_err = 0;

    list[0] = 1;
    list[1] = varno;
    varlist_adjust_sample(list, &t1, &t2, (const double **) *pZ);
    nobs = t2 - t1 + 1;

    if (missvals(&(*pZ)[varno][t1], nobs)) {
	pprintf(prn, "\n%s",
		_("Missing values within sample -- can't do correlogram"));
	return 1;
    }

    if (nobs < 4) {
	pputs(prn, _("\nInsufficient observations for correlogram"));
	return 1;
    }

    if (gretl_isconst(t1, t2, &(*pZ)[varno][0])) {
	sprintf(gretl_tmp_str, _("%s is a constant"), 
		pdinfo->varname[varno]);
	pprintf(prn, "\n%s\n", gretl_tmp_str);
	return 1;
    }

    acf_m = order;
    if (acf_m == 0) {
	acf_m = auto_acf_order(pdinfo->pd, nobs);
    } else if (acf_m > nobs - pdinfo->pd) {
	acf_m = nobs - 1; /* ?? */
    }

    acf = malloc(acf_m * sizeof *acf);
    if (acf == NULL) {
	return E_ALLOC;    
    }

    /* calculate acf, with lag order m */
    for (l=1; l<=acf_m; l++) {
	acf[l-1] = gretl_acf(nobs, l, &(*pZ)[varno][t1]);
    }

    if (opt & OPT_R) {
	pprintf(prn, "\n%s\n\n", _("Residual autocorrelation function"));
    } else {
	sprintf(gretl_tmp_str, _("Autocorrelation function for %s"), 
		pdinfo->varname[varno]);
	pprintf(prn, "\n%s\n\n", gretl_tmp_str);
    }

    /* add Ljung-Box statistic */
    box = 0;
    for (t=0; t<acf_m; t++) { 
	box += acf[t] * acf[t] / (nobs - t + 1);
    }
    box *= nobs * (nobs + 2.0);

    pprintf(prn, "Ljung-Box Q' = %.4f\n", box);
    pprintf(prn, _("Degrees of freedom = %d, p-value = %.4f\n\n"),
	    acf_m, chisq(box, acf_m));

    /* print acf */
    for (t=0; t<acf_m; t++) {
	pprintf(prn, "%5d)%8.4f", t + 1, acf[t]);
	if ((t + 1) % 5 == 0) {
	    pputc(prn, '\n');
	}
    }
    pputc(prn, '\n');

    if (opt & OPT_A) { 
	/* use ASCII graphics, not gnuplot */
	double *xl = malloc(acf_m * sizeof *xl);

	if (xl == NULL) {
	    err = E_ALLOC;
	    goto acf_getout;
	}
	for (l=0; l<acf_m; l++) {
	    xl[l] = l + 1.0;
	}
        pprintf(prn, "\n\n%s\n\n", _("Correlogram"));
	graphyzx(NULL, acf, NULL, xl, acf_m, pdinfo->varname[varno], 
		 _("lag"), NULL, 0, prn);
	free(xl);
    } 

    /* determine lag order for pacf (may have to be shorter than acf_m) */
    if (acf_m > nobs / 2 - 1) {
	pacf_m = nobs / 2 - 1;
    } else {
	pacf_m = acf_m;
    }

    /* generate (and if not in batch mode) plot partial 
       autocorrelation function */

    pacf = malloc(pacf_m * sizeof *pacf); 
    if (pacf == NULL) {
	err = E_ALLOC;
	goto acf_getout;
    }

    /* for confidence bands */
    pm = 1.96 / sqrt((double) nobs);

    err = pacf_err = get_pacf(pacf, acf, pacf_m);

    if (!err) {
	pacf[0] = acf[0];
	pprintf(prn, "\n%s", _("Partial autocorrelations"));
	if (pacf_m < acf_m) {
	    pprintf(prn, " (%s %d):\n\n", _("to lag"), pacf_m);
	} else {
	    pputs(prn, ":\n\n");
	}
	for (k=0; k<pacf_m; k++) {
	    pprintf(prn, "%5d)%8.4f", k+1, pacf[k]);
	    if ((k + 1) % 5 == 0) {
		pputc(prn, '\n');
	    }
	}
    }
    pputc(prn, '\n');
    if (pacf_m % 5 > 0) {
	pputc(prn, '\n');
    }

    pprintf(prn, "%s: %s = %g\n", 
	    /* xgettext:no-c-format */
	    _("5% critical value"),
	    corrgm_crit_string(),
	    pm);

    if (opt & OPT_A) {
	goto acf_getout;
    } else if (gnuplot_init(PLOT_CORRELOGRAM, &fq)) {
	err = E_FOPEN;
	goto acf_getout;
    }

    gretl_push_c_numeric_locale();

    /* create two separate plots, if both are OK */
    if (!pacf_err) {
	fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.48\n", fq);
    }
    fputs("set xzeroaxis\n", fq);
    fputs("set key top right\n", fq); 
    fprintf(fq, "set xlabel '%s'\n", I_("lag"));
    fputs("set yrange [-1.1:1.1]\n", fq);

    /* upper plot: Autocorrelation Function or ACF */
    if (!pacf_err) {
	fputs("set origin 0.0,0.50\n", fq);
    }
    if (opt & OPT_R) {
	fprintf(fq, "set title '%s'\n", I_("Residual ACF"));
    } else {
	fprintf(fq, "set title '%s %s'\n", I_("ACF for"), 
		pdinfo->varname[varno]);
    }
    fprintf(fq, "set xrange [0:%d]\n", acf_m + 1);
    fprintf(fq, "plot \\\n"
	    "'-' using 1:2 notitle w impulses, \\\n"
	    "%g title '+- %s' lt 2, \\\n"
	    "%g notitle lt 2\n", pm, corrgm_crit_string(), -pm);
    for (k=0; k<acf_m; k++) {
	fprintf(fq, "%d %g\n", k + 1, acf[k]);
    }
    fputs("e\n", fq);

    if (!pacf_err) {
	/* lower plot: Partial Autocorrelation Function or PACF */
	fputs("set origin 0.0,0.0\n", fq);
	if (opt & OPT_R) {
	    fprintf(fq, "set title '%s'\n", I_("Residual PACF"));
	} else {
	    fprintf(fq, "set title '%s %s'\n", I_("PACF for"), 
		    pdinfo->varname[varno]);
	}
	fprintf(fq, "set xrange [0:%d]\n", pacf_m + 1);
	fprintf(fq, "plot \\\n"
		"'-' using 1:2 notitle w impulses, \\\n"
		"%g title '+- %s' lt 2, \\\n"
		"%g notitle lt 2\n", pm, corrgm_crit_string(), -pm);
	for (k=0; k<pacf_m; k++) {
	    fprintf(fq, "%d %g\n", k + 1, pacf[k]);
	}
	fputs("e\n", fq);
    }

    if (!pacf_err) {
	fputs("set nomultiplot\n", fq);
    }

    gretl_pop_c_numeric_locale();

    fclose(fq);

    err = gnuplot_make_graph();

 acf_getout:

    free(acf);
    free(pacf);

    return err;
}

static int roundup_mod (int i, double x)
{
    return (int) ceil((double) x * i);
}

static int fract_int_GPH (int n, double *hhat, double *omega, PRN *prn)
{
    double xx, tstat, **tmpZ = NULL;
    DATAINFO *tmpdinfo;
    MODEL tmp;
    int t, err = 0, list[4];

    tmpdinfo = datainfo_new();
    if (tmpdinfo == NULL) {
	return 1;
    }

    tmpdinfo->n = n;
    tmpdinfo->v = 3;
    if (start_new_Z(&tmpZ, tmpdinfo, 1)) {
	return 1;
    }

    /* Test from Geweke and Porter-Hudak, as set out in
       Greene, Econometric Analysis 4e, p. 787 */
    
    for (t=0; t<n; t++) {
	tmpZ[0][t] = 1.0;
	tmpZ[1][t] = log(hhat[t]);
	xx = sin(omega[t] / 2);
	tmpZ[2][t] = log(4 * xx * xx);
    }

    list[0] = 3;
    list[1] = 1;
    list[2] = 0;
    list[3] = 2;

    tmp = lsq(list, &tmpZ, tmpdinfo, OLS, OPT_A);

    if (!tmp.errcode) {
	tstat = -tmp.coeff[1] / tmp.sderr[1];
	pprintf(prn, "\n%s (Geweke, Porter-Hudak)\n"
		"  %s = %g (%g)\n"
		"  %s: t(%d) = %g, %s %.4f\n",
		_("Test for fractional integration"),
		_("Estimated degree of integration"), -tmp.coeff[1], tmp.sderr[1],
		_("test statistic"), tmp.dfd, tstat, 
		_("with p-value"), t_pvalue_2(tstat, tmp.dfd));
    } else {
	err = tmp.errcode;
    }

    clear_model(&tmp);

    destroy_dataset(tmpZ, tmpdinfo);

    return err;
}

/* Fractional integration via Local Whittle Estimator */

gretl_matrix *gretl_matrix_periodogram (const gretl_matrix *x, int m)
{
    gretl_matrix *p;
    double *autocov;
    double xbar, varx;
    double xx, yy;
    int k, T, t; 

    T = gretl_vector_get_length(x);

    p = gretl_column_vector_alloc(m);
    if (p == NULL) {
	return NULL;
    }

    autocov = malloc(T * sizeof *autocov);
    if (autocov == NULL) {
	gretl_matrix_free(p);
	return NULL;
    }

    xbar = gretl_vector_mean(x);
    varx = gretl_vector_variance(x);

#if LWE_DEBUG
    fprintf(stderr, "gretl_matrix_periodogram: T = %d, m = %d\n"
	    "  xbar = %g, varx = %g\n", T, m, xbar, varx);
#endif

    /* find autocovariances */
    for (k=1; k<=T-1; k++) {
	double xt, xtk;

	autocov[k] = 0.0;
	for (t=k; t<T; t++) {
	    xt = gretl_vector_get(x, t);
	    xtk = gretl_vector_get(x, t - k);
	    autocov[k] += (xt - xbar) * (xtk - xbar);
	}
	autocov[k] /= T;

#if LWE_DEBUG > 1
	fprintf(stderr, "autocov[%d] = %g\n", k, autocov[k]);
#endif
    }

    for (t=1; t<=m; t++) {
	yy = 2 * M_PI * t / (double) T;
	xx = varx; 
	for (k=1; k<=T-1; k++) {
	    xx += 2.0 * autocov[k] * cos(yy * k);
	}
	xx /= 2 * M_PI;
	xx *= T;
	gretl_vector_set(p, t-1, xx);
#if LWE_DEBUG
	fprintf(stderr, "periodogram[%d] = %g\n", t, xx);
#endif
    }

    free(autocov);

    return p;
}

gretl_matrix *LWE_lambda (const gretl_matrix *I, int n, double *lcm)
{
    int m = gretl_vector_get_length(I);
    gretl_matrix *lambda, *llambda;
    int i;

    lambda = gretl_column_vector_alloc(m);

    for (i=0; i<m; i++) {
	gretl_vector_set(lambda, i, (2.0 * M_PI / n) * (i + 1));
#if LWE_DEBUG
	fprintf(stderr, "LWE_obj_func: lambda[%d] = %g\n",
		i, gretl_vector_get(lambda, i));
#endif
    }

    llambda = gretl_matrix_copy(lambda);
    gretl_matrix_transform_elements(llambda, T_LOG);
    *lcm = gretl_vector_mean(llambda);

#if LWE_DEBUG
    fprintf(stderr, "LWE_lambda: col mean of log lambda = %g\n", *lcm);
#endif

    gretl_matrix_free(llambda);

    return lambda;
}

double LWE_obj_func (const gretl_matrix *I, double d,
		     const gretl_matrix *lambda, double lcm)
{
    gretl_matrix *lambda2, *Itmp;
    double dd = 2.0 * d;
    double ret;
    int err = 0;

    lambda2 = gretl_matrix_copy(lambda);
    if (lambda2 == NULL) {
	return NADBL;
    }

    gretl_matrix_dot_pow(lambda2, dd);

    Itmp = gretl_matrix_dot_multiply(I, lambda2, &err);
    if (Itmp == NULL) {
	gretl_matrix_free(lambda2);
	return NADBL;
    }

    ret = -(log(gretl_vector_mean(Itmp)) - dd * lcm);
    
    gretl_matrix_free(lambda2);
    gretl_matrix_free(Itmp);

    return ret;
}

double LWE (const gretl_matrix *X, int m)
{
    gretl_matrix *I;
    gretl_matrix *lambda;
    int n = gretl_matrix_rows(X);
    double d = 0, ret;
    double lcm;

    double dd = 1.0, f, incr, incl, deriv, h, eps = 1.0e-05;
    int iter = 0;
    const int MAX_ITER = 100;
      
    I = gretl_matrix_periodogram(X, m);
    if (I == NULL) {
	return NADBL;
    }

    lambda = LWE_lambda(I, n, &lcm);
    if (lambda == NULL) {
	gretl_matrix_free(I);
	return NADBL;
    }    

    while (fabs(dd) > 1.0e-06 && iter < MAX_ITER) {
	f = LWE_obj_func(I, d, lambda, lcm);
	incr = LWE_obj_func(I, d + eps, lambda, lcm) / eps;
	incl = LWE_obj_func(I, d - eps, lambda, lcm) / eps;
	deriv = (incr - incl) / 2.0;
	h = (0.5 * (incr + incl) - f / eps) / eps;

	if (h >= 0) {
	    dd = deriv;
	} else {
	    dd = -deriv / h;
	}

	if (fabs(dd) > 1) {
	    dd = (dd > 0) ? 1 : -1;
	}
	
#if LWE_DEBUG
	fprintf(stderr, "d = %g, f(d) = %g, deriv = %g, h = %g, dd = %g\n",
		d, f, deriv, h, dd);
#endif
	d += 0.5 * dd;
	iter++;
    }

    if (iter == MAX_ITER) {
	fprintf(stderr, "Maximum number of iterations reached\n");
	ret = NADBL;
    } else {
        ret = d;
    }

    gretl_matrix_free(I);
    gretl_matrix_free(lambda);

    return ret;
}

int fract_int_LWE (const double **Z, int varno, int t1, int t2,
		   PRN *prn)
{
    gretl_matrix *X;
    double m1, m2;
    double d, se, z;
    int T, m;

    X = gretl_data_series_to_vector(Z, varno, t1, t2);
    if (X == NULL) {
	return 1;
    }

    T = gretl_vector_get_length(X);

    m1 = floor((double) T / 2.0);
    m2 = floor(pow((double) T, 0.6));

    if (m1 < m2) {
	m = m1;
    } else {
	m = m2;
    }

    m--;

    d = LWE(X, m);
    if (na(d)) {
	gretl_matrix_free(X);
	return 1;
    }

    se = 1 / (2.0 * sqrt((double) m));
    z = d / se;

    pprintf(prn, "\n%s (T = %d, m = %d)\n"
	    "  %s = %g (%g)\n"
	    "  %s: z = %g, %s %.4f\n\n",
	    _("Local Whittle Estimator"), T, m,
	    _("Estimated degree of integration"), d, se,
	    _("test statistic"), z, 
	    _("with p-value"), normal_pvalue_2(z));    

    gretl_matrix_free(X);

    return 0;
}

/**
 * periodogram:
 * @varno: ID number of variable to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if includes %OPT_O, use Bartlett lag window for periodogram;
 * if includes %OPT_N, don't display gnuplot graph; if includes
 * %OPT_R, the variable is a model residual.
 * @prn: gretl printing struct.
 *
 * Computes and displays the periodogram for the variable specified 
 * by @varno.
 *
 * Returns: 0 on successful completion, error code on error.
 *
 */

int periodogram (int varno, double ***pZ, const DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    double *autocov = NULL;
    double *omega = NULL;
    double *hhat = NULL;
    double *savexx = NULL;
    double *stdy = NULL;
    double xx, yy, varx, stdx, w;
    int k, xmax, L, nT, nobs; 
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int list[2];
    int do_graph = !(opt & OPT_N);
    int window = (opt & OPT_O);
    FILE *fq = NULL;
    int err = 0;

    *gretl_errmsg = 0;

    list[0] = 1;
    list[1] = varno;
    varlist_adjust_sample(list, &t1, &t2, (const double **) *pZ);
    nobs = t2 - t1 + 1;

    if (missvals(&(*pZ)[varno][t1], nobs)) {
	strcpy(gretl_errmsg, 
	       _("Missing values within sample -- can't do periodogram"));
	return 1;
    }    

    if (nobs < 12) {
	strcpy(gretl_errmsg,
	       _("Insufficient observations for periodogram"));
	return 1;
    }

    if (gretl_isconst(t1, t2, &(*pZ)[varno][0])) {
	sprintf(gretl_tmp_str, _("'%s' is a constant"), pdinfo->varname[varno]);
	pprintf(prn, "\n%s\n", gretl_tmp_str);
	return 1;
    }

    /* Chatfield (1996); Greene 4ed, p. 772 */
    if (window) {
	L = (int) 2.0 * sqrt((double) nobs);
    } else {
	L = nobs - 1; 
    }

    /* prepare for fractional integration test */
    xx = sqrt((double) nobs);
    nT = (int) xx;
    if ((double) nT < xx) {
	nT++;
    }
    
    autocov = malloc((L + 1) * sizeof *autocov);
    omega = malloc(nT * sizeof *omega);
    hhat = malloc(nT * sizeof *hhat);
    stdy = malloc(nobs * sizeof *stdy);

    if (autocov == NULL || omega == NULL || hhat == NULL || stdy == NULL) {
	free(autocov);
	free(omega);
	free(hhat);
	free(stdy);
	return E_ALLOC;
    }

    xx = gretl_mean(t1, t2, (*pZ)[varno]);
    varx = gretl_variance(t1, t2, &(*pZ)[varno][0]);
    varx *= (double) (nobs - 1) / nobs;
    stdx = sqrt(varx);
    for (t=t1; t<=t2; t++) {
	stdy[t-t1] = ((*pZ)[varno][t] - xx) / stdx;
    }

    /* find autocovariances */
    for (k=1; k<=L; k++) {
	autocov[k] = 0.0;
	for (t=k; t<nobs; t++) {
	    autocov[k] += stdy[t] * stdy[t-k];
	}
	autocov[k] /= nobs;
    }

    xmax = roundup_mod(nobs, 2.0);

    if (do_graph && gnuplot_init(PLOT_PERIODOGRAM, &fq) == 0) {
	char titlestr[80];

	fputs("# literal lines = 4\n", fq);
	fputs("set xtics nomirror\n", fq); 

	if (pdinfo->pd == 4) {
	    fprintf(fq, "set x2label '%s'\n", I_("quarters"));
	} else if (pdinfo->pd == 12) {
	    fprintf(fq, "set x2label '%s'\n", I_("months"));
	} else if (pdinfo->pd == 1 && pdinfo->structure == TIME_SERIES) {
	    fprintf(fq, "set x2label '%s'\n", I_("years"));
	} else {
	    fprintf(fq, "set x2label '%s'\n", I_("periods"));
	}

	fprintf(fq, "set x2range [0:%d]\n", xmax);
	fputs("set x2tics(", fq);
	k = (nobs / 2) / 6;
	for (t = 1; t <= nobs/2; t += k) {
	    fprintf(fq, "\"%.1f\" %d, ", (double) nobs / t, 4 * t);
	}
	fprintf(fq, "\"\" %d)\n", 2 * nobs);
	fprintf(fq, "set xlabel '%s'\n", I_("scaled frequency"));
	fputs("set xzeroaxis\n", fq);
	fputs("set nokey\n", fq);

	if (opt & OPT_R) {
	    strcpy(titlestr, I_("Residual spectrum"));
	} else {
	    sprintf(titlestr, I_("Spectrum of %s"), pdinfo->varname[varno]);
	}
	fprintf(fq, "set title '%s", titlestr);

	if (window) {
	    sprintf(titlestr, I_("Bartlett window, length %d"), L);
	    fprintf(fq, " (%s)'\n", titlestr);
	} else {
	    fputs("'\n", fq);
	}

	fprintf(fq, "set xrange [0:%d]\n", roundup_mod(nobs, 0.5));
	fputs("plot '-' using 1:2 w lines\n", fq);
    }

    if (do_graph && fq == NULL) {
	do_graph = 0;
	err = 1;
    }

    if (opt & OPT_R) {
	pprintf(prn, "\n%s\n", _("Residual periodogram"));
    } else {
	pprintf(prn, _("\nPeriodogram for %s\n"), pdinfo->varname[varno]);
    }

    pprintf(prn, _("Number of observations = %d\n"), nobs);
    if (window) {
	pprintf(prn, _("Using Bartlett lag window, length %d\n\n"), L);
    }
    pputs(prn, _(" omega  scaled frequency  periods  spectral density\n\n"));

    if (do_graph) { 
	savexx = malloc((1 + nobs/2) * sizeof *savexx);
	if (savexx == NULL) {
	    err = 1;
	    fclose(fq);
	    do_graph = 0;
	}
    }

    for (t=1; t<=nobs/2; t++) {
	yy = 2 * M_PI * t / (double) nobs;
	xx = 1.0; 
	for (k=1; k<=L; k++) {
	    if (window) {
		w = 1 - (double) k/(L + 1);
	    } else {
		w = 1.0;
	    }
	    xx += 2.0 * w * autocov[k] * cos(yy * k);
	}
	xx *= varx /(2 * M_PI);
	pprintf(prn, " %.4f%9d%16.2f%16.5f\n", yy, t, 
		(double) nobs / t, xx);
	if (savexx != NULL) {
	    savexx[t] = xx;
	}
	if (t <= nT) {
	    omega[t-1] = yy;
	    hhat[t-1] = xx;
	}
    }

    pputc(prn, '\n');

    if (do_graph) {
	gretl_push_c_numeric_locale();

	for (t=1; t<=nobs/2; t++) {
	    fprintf(fq, "%d %g\n", t, savexx[t]);
	}

	gretl_pop_c_numeric_locale();

	fputs("e\n", fq);

	fclose(fq);
	free(savexx);
	err = gnuplot_make_graph();
    }

    if (!window) {
	if (fract_int_GPH(nT, hhat, omega, prn)) {
	    pprintf(prn, "\n%s\n", _("Fractional integration test failed"));
	}
	fract_int_LWE((const double **) *pZ, varno, t1, t2, prn);
    }

    free(autocov);
    free(omega);
    free(hhat);
    free(stdy);

    return err;
}

static void printf15 (double zz, PRN *prn)
{
    if (na(zz)) {
	pprintf(prn, "%*s", UTF_WIDTH(_("undefined"), 15), 
		_("undefined"));
    } else {
	pputc(prn, ' ');
	gretl_print_fullwidth_double(zz, 5, prn);	
    }
}

#define LINEWID 78

static void center_line (char *str, PRN *prn, int dblspc)
{
    size_t len = strlen(str);

    if (LINEWID > len) {
	size_t i, pad = (LINEWID - len) / 2;
	char cstr[84];

	for (i=0; i<pad; i++) cstr[i] = ' ';
	strcpy(cstr + i, str);
	if (dblspc) {
	    strcat(cstr, "\n");
	}
	pprintf(prn, "%s\n", cstr);
    } else {
	if (dblspc) {
	    strcat(str, "\n");
	}
	pprintf(prn, "%s\n", str);
    }
}

static void prhdr (const char *str, const DATAINFO *pdinfo, 
		   int ci, int missing, PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN], tmp[96];

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    pputc(prn, '\n');

    sprintf(tmp, _("%s, using the observations %s - %s"), str, date1, date2);
    center_line(tmp, prn, 0);

    if (missing) {
	strcpy(tmp, _("(missing values were skipped)"));
	center_line(tmp, prn, 1);
    }
}

static void print_summary_single (const Summary *summ, int j,
				  const DATAINFO *pdinfo,
				  PRN *prn)
{
    char obs1[OBSLEN], obs2[OBSLEN], tmp[128];
    double vals[8];
    const char *labels[] = {
	N_("Mean"),
	N_("Median"),
	N_("Minimum"),
	N_("Maximum"),
	N_("Standard deviation"),
	N_("C.V."),
	N_("Skewness"),
	N_("Ex. kurtosis")
    };
    int slen = 0, i = 0;

    ntodate(obs1, pdinfo->t1, pdinfo);
    ntodate(obs2, pdinfo->t2, pdinfo);

    prhdr(_("Summary Statistics"), pdinfo, SUMMARY, 0, prn);
    sprintf(tmp, _("for the variable '%s' (%d valid observations)"), 
	    pdinfo->varname[summ->list[j+1]], summ->n);
    center_line(tmp, prn, 1);

    vals[0] = summ->mean[j];
    vals[1] = summ->median[j];
    vals[2] = summ->low[j];
    vals[3] = summ->high[j];
    vals[4] = summ->sd[j];
    vals[5] = summ->cv[j];
    vals[6] = summ->skew[j];
    vals[7] = summ->xkurt[j];

    for (i=0; i<8; i++) {
	if (strlen(_(labels[i])) > slen) {
#if defined(ENABLE_NLS) && defined(USE_GTK2)
	    slen = g_utf8_strlen(_(labels[i]), -1);	    
#else
	    slen = strlen(_(labels[i]));
#endif
	}
    }
    slen++;

    for (i=0; i<8; i++) {
	pprintf(prn, "  %-*s", UTF_WIDTH(_(labels[i]), slen), _(labels[i]));
	printf15(vals[i], prn);
	pputc(prn, '\n');
    }

    pputs(prn, "\n\n");    
}

/**
 * print_summary:
 * @summ: gretl summary statistics struct.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 *
 * Print the summary statistics for a given variable.
 *
 */

void print_summary (const Summary *summ,
		    const DATAINFO *pdinfo,
		    PRN *prn)
{
    int pause = gretl_get_text_pause();
    int len, maxlen = 0;
    int i, vi, lineno;

    if (summ->list[0] == 1) {
	print_summary_single(summ, 0, pdinfo, prn);
	return;
    }

    for (i=1; i<=summ->list[0]; i++) {
	vi = summ->list[i];
	len = strlen(pdinfo->varname[vi]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    len = (maxlen <= 8)? 10 : maxlen + 1;

    if (len > 13) {
	for (i=0; i<summ->list[0]; i++) {
	    print_summary_single(summ, i, pdinfo, prn);
	}
	return;
    }

    prhdr(_("Summary Statistics"), pdinfo, SUMMARY, summ->missing, prn);

    pprintf(prn, "\n%s  ", _("Variable"));
    pputs(prn, _("      MEAN           MEDIAN           MIN"
            "             MAX\n\n"));

    lineno = 1;
    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (pause && (lineno % PAGELINES == 0)) {
	    scroll_pause();
	    lineno = 1;
	}
	pprintf(prn, "%-*s", len, pdinfo->varname[vi]);
	printf15(summ->mean[i], prn);
	printf15(summ->median[i], prn);
	printf15(summ->low[i], prn);
	printf15(summ->high[i], prn);
	pputc(prn, '\n');
	lineno++;
    }
    pputc(prn, '\n');

    if (pause) {
	scroll_pause();
    }

    pprintf(prn, "\n%s  ", _("Variable"));
    pputs(prn, _("      S.D.            C.V.           "
	 " SKEW          EXCSKURT\n\n"));

    lineno = 1;
    for (i=0; i<summ->list[0]; i++) {
	double cv;

	vi = summ->list[i + 1];

	if (pause && (lineno % PAGELINES == 0)) {
	    scroll_pause();
	    lineno = 1;
	}

	pprintf(prn, "%-*s", len, pdinfo->varname[vi]);

	if (floateq(summ->mean[i], 0.0)) {
	    cv = NADBL;
	} else if (floateq(summ->sd[i], 0.0)) {
	    cv = 0.0;
	} else {
	    cv = fabs(summ->sd[i] / summ->mean[i]);
	} 

	printf15(summ->sd[i], prn);
	printf15(cv, prn);
	printf15(summ->skew[i], prn);
	printf15(summ->xkurt[i], prn);
	pputc(prn, '\n');
	lineno++;
    }
    pputc(prn, '\n');
}

/**
 * free_summary:
 * @summ: gretl summary statistics struct
 *
 * Frees all malloced elements of the struct.
 *
 */

void free_summary (Summary *summ)
{
    free(summ->list);
    free(summ->stats);

    free(summ);
}

static Summary *summary_new (const int *list)
{
    Summary *summ;
    int nv = list[0];

    summ = malloc(sizeof *summ);
    if (summ == NULL) {
	return NULL;
    }

    summ->list = gretl_list_copy(list);
    if (summ->list == NULL) {
	free(summ);
	return NULL;
    }

    summ->n = 0;
    summ->missing = 0;

    summ->stats = malloc(8 * nv * sizeof *summ->stats);
    if (summ->stats == NULL) {
	free_summary(summ);
	return NULL;
    }

    summ->mean = summ->stats;
    summ->median = summ->mean + nv;
    summ->sd = summ->median + nv;
    summ->skew = summ->sd + nv;
    summ->xkurt = summ->skew + nv;
    summ->low = summ->xkurt + nv;
    summ->high = summ->low + nv;
    summ->cv = summ->high + nv;

    return summ;
}

/**
 * summary:
 * @list: list of variables to process.
 * @Z: data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 *
 * Calculates descriptive summary statistics for the specified variables.
 *
 * Returns: struct containing the summary statistics.
 *
 */

Summary *summary (const int *list, 
		       const double **Z, const DATAINFO *pdinfo,
		       PRN *prn) 
{
    Summary *summ;
    int i, vi, sn, gn;

    summ = summary_new(list);
    if (summ == NULL) {
	return NULL;
    }

    for (i=0; i<summ->list[0]; i++)  {
	double x0;

	vi = summ->list[i + 1];

	sn = pdinfo->t2 - pdinfo->t1 + 1;
	gn = good_obs(Z[vi] + pdinfo->t1, sn, &x0);

	if (gn < sn) {
	    summ->missing = 1;
	}

	if (gn > summ->n) {
	    summ->n = gn;
	}

	if (sn < 2) { 
	    /* zero or one observations */
	    if (summ->n == 0) {
		pprintf(prn, _("Dropping %s: sample range contains no valid "
			"observations\n"), pdinfo->varname[vi]);
	    } else {
		pprintf(prn, _("Dropping %s: sample range has only one "
			"obs, namely %g\n"), pdinfo->varname[vi], x0);
	    }
	    gretl_list_delete_at_pos(summ->list, i + 1);
	    if (summ->list[0] == 0) {
		free_summary(summ);
		return NULL;
	    } else {
		i--;
		continue;
	    }
	}

	gretl_minmax(pdinfo->t1, pdinfo->t2, Z[vi], 
		     &summ->low[i], 
		     &summ->high[i]);
	
	gretl_moments(pdinfo->t1, pdinfo->t2, Z[vi], 
		      &summ->mean[i], 
		      &summ->sd[i], 
		      &summ->skew[i], 
		      &summ->xkurt[i], 1);

	if (!floateq(summ->mean[i], 0.0)) {
	    summ->cv[i] = fabs(summ->sd[i] / summ->mean[i]);
	} else {
	    summ->cv[i] = NADBL;
	}

	summ->median[i] = gretl_median(pdinfo->t1, pdinfo->t2, Z[vi]);
    } 

    return summ;
}

/**
 * vmatrix_new:
 *
 * Returns: an allocated and initialized #VMatrix, or
 * %NULL on failure.
 */

VMatrix *vmatrix_new (void)
{
    VMatrix *vmat = malloc(sizeof *vmat);

    if (vmat != NULL) {
	vmat->vec = NULL;
	vmat->list = NULL;
	vmat->names = NULL;

	vmat->ci = 0;
	vmat->dim = 0;
	vmat->t1 = 0;
	vmat->t2 = 0;
	vmat->missing = 0;
    }

    return vmat;
}

/**
 * free_vmatrix:
 * @vmat: gretl correlation matrix struct
 *
 * Frees all malloced elements of the struct.
 */

void free_vmatrix (VMatrix *vmat)
{
    if (vmat != NULL) {
	free_strings_array(vmat->names, vmat->dim);
	free(vmat->vec);
	if (vmat->list != NULL) {
	    free(vmat->list);
	}
	free(vmat);
    }
}

/**
 * corrlist:
 * @list: list of variables to process, by ID number.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 *
 * Computes pairwise correlation coefficients for the variables
 * specified in @list, skipping any constants.
 *
 * Returns: gretl correlation matrix struct, or %NULL on failure.
 */

VMatrix *corrlist (int *list, const double **Z, const DATAINFO *pdinfo)
{
    VMatrix *corrmat;
    int i, j, lo, nij, mm;
    int t1 = pdinfo->t1, t2 = pdinfo->t2; 
    int missing = 0;

    corrmat = vmatrix_new();
    if (corrmat == NULL) {
	return NULL;
    }

    /* drop any constants from list */
    for (i=1; i<=list[0]; i++) {
	if (gretl_isconst(t1, t2, Z[list[i]])) {
	    gretl_list_delete_at_pos(list, i);
	    i--;
	}
    }

    corrmat->dim = lo = list[0];  
    mm = (lo * (lo + 1)) / 2;

    corrmat->names = create_strings_array(lo);
    if (corrmat->names == NULL) {
	free(corrmat);
	return NULL;
    }

    corrmat->vec = malloc(mm * sizeof *corrmat->vec);
    if (corrmat->vec == NULL) {
	free_vmatrix(corrmat);
	return NULL;
    }

    for (i=0; i<lo; i++) {  
	int vi = list[i+1];

	corrmat->names[i] = gretl_strdup(pdinfo->varname[vi]);
	if (corrmat->names[i] == NULL) {
	    free_vmatrix(corrmat);
	    return NULL;
	}

	for (j=0; j<lo; j++)  {
	    int vj = list[j+1];

	    nij = ijton(i, j, lo);
	    if (i == j) {
		corrmat->vec[nij] = 1.0;
		continue;
	    }
	    corrmat->vec[nij] = gretl_corr(t1, t2, Z[vi], Z[vj], &missing);
	    if (missing > 0) {
		corrmat->missing = 1;
	    }
	}
    }

    corrmat->list = gretl_list_copy(list);
    corrmat->ci = CORR;
    corrmat->t1 = t1;
    corrmat->t2 = t2;

    return corrmat;
}

/**
 * matrix_print_corr:
 * @corr: gretl correlation matrix
 * @pdinfo: dataset information.
 * @prn: gretl printing struct.
 *
 * Prints a gretl correlation matrix to @prn.
 */

void matrix_print_corr (VMatrix *corr, const DATAINFO *pdinfo, PRN *prn)
{
    char tmp[96];
    int n = corr->t2 - corr->t1 + 1;

    prhdr(_("Correlation Coefficients"), pdinfo, CORR, corr->missing, 
	  prn);
    sprintf(tmp, _("5%% critical value (two-tailed) = "
	    "%.4f for n = %d"), rhocrit95(n), n);
    center_line(tmp, prn, 1);
    text_print_vmatrix(corr, prn);
}

/**
 * gretl_corrmx:
 * @list: gives the ID numbers of the variables to process.
 * @Z: data array.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Computes and prints the correlation matrix for the specified list
 * of variables.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int gretl_corrmx (int *list, const double **Z, const DATAINFO *pdinfo, 
		  PRN *prn)
{
    VMatrix *corr;

    corr = corrlist(list, Z, pdinfo);
    if (corr == NULL) {
	return 1;
    }

    matrix_print_corr(corr, pdinfo, prn);
    free_vmatrix(corr);

    return 0;
}

/**
 * means_test:
 * @list: gives the ID numbers of the variables to compare.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @opt: if OPT_O, assume population variances are different.
 * @prn: gretl printing struct.
 *
 * Carries out test of the null hypothesis that the means of two
 * variables are equal.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int means_test (const int *list, const double **Z, const DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn)
{
    double m1, m2, s1, s2, skew, kurt, se, mdiff, t, pval;
    double v1, v2;
    double *x = NULL, *y = NULL;
    int vardiff = (opt & OPT_O);
    int df, n1, n2, n = pdinfo->n;

    if (list[0] < 2) {
	return E_ARGS;
    }

    x = malloc(n * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    y = malloc(n * sizeof *y);
    if (y == NULL) {
	free(x);
	return E_ALLOC;
    }    

    n1 = ztox(list[1], x, Z, pdinfo);
    n2 = ztox(list[2], y, Z, pdinfo);

    if (n1 == 0 || n2 == 0) {
	pputs(prn, _("Sample range has no valid observations."));
	free(x); free(y);
	return 1;
    }

    if (n1 == 1 || n2 == 1) {
	pputs(prn, _("Sample range has only one observation."));
	free(x); free(y);
	return 1;
    }

    df = n1 + n2 - 2;

    gretl_moments(0, n1-1, x, &m1, &s1, &skew, &kurt, 1);
    gretl_moments(0, n2-1, y, &m2, &s2, &skew, &kurt, 1);
    mdiff = m1 - m2;

    v1 = s1 * s1;
    v2 = s2 * s2;

    if (vardiff) {
	se = sqrt((v1 / n1) + (v2 / n2));
    } else {
	/* form pooled estimate of variance */
	double sp2;

	sp2 = ((n1-1) * v1 + (n2-1) * v2) / df;
	se = sqrt(sp2 / n1 + sp2 / n2);
    }

    t = mdiff / se;
    pval = t_pvalue_2(t, df);

    pprintf(prn, _("\nEquality of means test "
	    "(assuming %s variances)\n\n"), (vardiff)? _("unequal") : _("equal"));
    pprintf(prn, "   %s: ", pdinfo->varname[list[1]]);
    pprintf(prn, _("Number of observations = %d\n"), n1);
    pprintf(prn, "   %s: ", pdinfo->varname[list[2]]);
    pprintf(prn, _("Number of observations = %d\n"), n2);
    pprintf(prn, _("   Difference between sample means = %g - %g = %g\n"), 
	    m1, m2, mdiff);
    pputs(prn, _("   Null hypothesis: The two population means are the same.\n"));
    pprintf(prn, _("   Estimated standard error = %g\n"), se);
    pprintf(prn, _("   Test statistic: t(%d) = %g\n"), df, t);
    pprintf(prn, _("   p-value (two-tailed) = %g\n\n"), pval);
    if (pval > .10)
	pputs(prn, _("   The difference is not statistically significant.\n\n"));

    record_test_result(t, pval, _("difference of means"));

    free(x);
    free(y);

    return 0;
}

/**
 * vars_test:
 * @list: gives the ID numbers of the variables to compare.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Carries out test of the null hypothesis that the variances of two
 * variables are equal.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int vars_test (const int *list, const double **Z, const DATAINFO *pdinfo, 
	       PRN *prn)
{
    double m, skew, kurt, s1, s2, var1, var2;
    double F, pval;
    double *x = NULL, *y = NULL;
    int dfn, dfd, n1, n2, n = pdinfo->n;

    if (list[0] < 2) return E_ARGS;

    if ((x = malloc(n * sizeof *x)) == NULL) return E_ALLOC;
    if ((y = malloc(n * sizeof *y)) == NULL) return E_ALLOC;

    n1 = ztox(list[1], x, Z, pdinfo);
    n2 = ztox(list[2], y, Z, pdinfo);

    if (n1 == 0 || n2 == 0) {
	pputs(prn, _("Sample range has no valid observations."));
	free(x); free(y);
	return 1;
    }

    if (n1 == 1 || n2 == 1) {
	pputs(prn, _("Sample range has only one observation."));
	free(x); free(y);
	return 1;
    }
    
    gretl_moments(0, n1-1, x, &m, &s1, &skew, &kurt, 1);
    gretl_moments(0, n2-1, y, &m, &s2, &skew, &kurt, 1);

    var1 = s1*s1;
    var2 = s2*s2;
    if (var1 > var2) { 
	F = var1/var2;
	dfn = n1 - 1;
	dfd = n2 - 1;
    } else {
	F = var2/var1;
	dfn = n2 - 1;
	dfd = n1 - 1;
    }

    pval = fdist(F, dfn, dfd);

    pputs(prn, _("\nEquality of variances test\n\n"));
    pprintf(prn, "   %s: ", pdinfo->varname[list[1]]);
    pprintf(prn, _("Number of observations = %d\n"), n1);
    pprintf(prn, "   %s: ", pdinfo->varname[list[2]]);
    pprintf(prn, _("Number of observations = %d\n"), n2);
    pprintf(prn, _("   Ratio of sample variances = %g\n"), F);
    pprintf(prn, "   %s: %s\n", _("Null hypothesis"), 
	    _("The two population variances are equal"));
    pprintf(prn, "   %s: F(%d,%d) = %g\n", _("Test statistic"), dfn, dfd, F);
    pprintf(prn, _("   p-value (two-tailed) = %g\n\n"), pval);
    if (fdist(F, dfn, dfd) > .10)
	pputs(prn, _("   The difference is not statistically significant.\n\n"));

    record_test_result(F, pval, _("difference of variances"));

    free(x);
    free(y);

    return 0;
}

struct MahalDist_ {
    int *list;
    int n;
    double *d;
};

const double *mahal_dist_get_distances (const MahalDist *md)
{
    return md->d;
}

int mahal_dist_get_n (const MahalDist *md)
{
    return md->n;
}

const int *mahal_dist_get_varlist(const MahalDist *md)
{
    return md->list;
}

void free_mahal_dist (MahalDist *md)
{
    free(md->list);
    free(md->d);
    free(md);
}

static MahalDist *mahal_dist_new (const int *list, int n)
{
    MahalDist *md = malloc(sizeof *md);
    
    if (md != NULL) {
	md->d = malloc(n * sizeof *md->d);
	if (md->d == NULL) {
	    free(md);
	    md = NULL;
	} else {
	    md->list = gretl_list_copy(list);
	    if (md->list == NULL) {
		free(md->d);
		free(md);
		md = NULL;
	    } else {
		md->n = n;
	    }
	}
    }

    if (md != NULL) {
	int t;

	for (t=0; t<n; t++) {
	    md->d[t] = NADBL;
	}
    }

    return md;
}

static int mdist_saver (double ***pZ, DATAINFO *pdinfo)
{
    int sv = 0;
    int err;

    err = dataset_add_series(1, pZ, pdinfo);

    if (!err) {
	int t;

	sv = pdinfo->v - 1;
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[sv][t] = NADBL;
	}

	strcpy(pdinfo->varname[sv], "mdist");
	make_varname_unique(pdinfo->varname[sv], sv, pdinfo);

	strcpy(VARLABEL(pdinfo, sv), _("Mahalanobis distances"));	
    }
		
    return sv;
}

static int 
real_mahalanobis_distance (const int *list, double ***pZ,
			   DATAINFO *pdinfo, gretlopt opt,
			   MahalDist *md, PRN *prn)
{
    gretl_matrix *S = NULL;
    gretl_vector *means = NULL;
    gretl_vector *xdiff;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    int n, err = 0;

    varlist_adjust_sample(list, &pdinfo->t1, &pdinfo->t2, 
			  (const double **) *pZ);

    n = pdinfo->t2 - pdinfo->t1 + 1;
    if (n < 2) {
	pdinfo->t1 = orig_t1;
	pdinfo->t2 = orig_t2;
	return E_DATA;
    }

    xdiff = gretl_column_vector_alloc(list[0]);
    if (xdiff == NULL) {
	pdinfo->t1 = orig_t1;
	pdinfo->t2 = orig_t2;
	return E_ALLOC;
    }

    S = gretl_covariance_matrix_from_varlist(list, 
					     (const double **) *pZ, 
					     pdinfo, 
					     &means,
					     &err);

    if (!err) {
	if (opt & OPT_V) {
	    gretl_matrix_print_to_prn(S, _("Covariance matrix"), prn);
	}
	err = gretl_invert_symmetric_matrix(S);
	if (err) {
	    fprintf(stderr, "error inverting covariance matrix\n");
	} else if (opt & OPT_V) {
	    gretl_matrix_print_to_prn(S, _("Inverse of covariance matrix"), prn);
	}
    }

    if (!err) {
	int k = gretl_vector_get_length(means);
	char obs_string[OBSLEN];
	int savevar = 0;
	double m;
	int i, t;

	if (opt & OPT_S) {
	    /* save the results to a data series */
	    savevar = mdist_saver(pZ, pdinfo);
	}

	pprintf(prn, "%s\n", _("Mahalanobis distances from the centroid"));
	pprintf(prn, "%s\n", _("using the variables:"));
	for (i=1; i<=list[0]; i++) {
	    pprintf(prn, " %s\n", pdinfo->varname[list[i]]);
	}
	pputc(prn, '\n');

	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {

	    /* write vector of deviations from centroid for
	       observation t */
	    for (i=0; i<k; i++) {
		int vi = list[i+1];
		double xbar = gretl_vector_get(means, i);

		gretl_vector_set(xdiff, i, (*pZ)[vi][t] - xbar);
	    }

	    m = gretl_scalar_b_X_b(xdiff, GRETL_MOD_TRANSPOSE, S, &err);

	    get_obs_string(obs_string, t, pdinfo);
	    pprintf(prn, "%8s ", obs_string);

	    if (err) {
		pprintf(prn, "NA\n");
	    } else {
		m = sqrt(m);
		pprintf(prn, "%9.6f\n", m);
		if (savevar > 0) {
		    (*pZ)[savevar][t] = m;
		} else if (md != NULL) {
		    md->d[t] = m;
		}
	    }
	}

	if (savevar > 0) {
	    pputc(prn, '\n');
	    pprintf(prn, _("Distances saved as '%s'"), 
		    pdinfo->varname[savevar]);
	    pputc(prn, '\n');
	}

	pputc(prn, '\n');
    }

    gretl_matrix_free(xdiff);
    gretl_matrix_free(means);
    gretl_matrix_free(S);

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    return err;
}

int mahalanobis_distance (const int *list, double ***pZ,
			  DATAINFO *pdinfo, gretlopt opt, 
			  PRN *prn)
{
    return real_mahalanobis_distance(list, pZ, pdinfo, opt, 
				     NULL, prn);
}

MahalDist *get_mahal_distances (const int *list, double ***pZ,
			     DATAINFO *pdinfo, gretlopt opt,
			     PRN *prn)
{
    MahalDist *md = mahal_dist_new(list, pdinfo->n);
    int err;

    if (md != NULL) {
	err = real_mahalanobis_distance(list, pZ, pdinfo, opt, 
					md, prn);
	if (err) {
	    free_mahal_dist(md);
	    md = NULL;
	}
    }

    return md;
}

static double gini_coeff (const double *x, int t1, int t2, double **plz,
			  int *pn, int *err)
{
    int m = t2 - t1 + 1;
    double *sx = NULL;
    double csx = 0.0, sumx = 0.0, sisx = 0.0;
    double idx, gini;
    int t, n = 0;

    *gretl_errmsg = '\0';

    sx = malloc(m * sizeof *sx);

    if (sx == NULL) {
	*err = E_ALLOC;
	return NADBL;
    }

    n = 0;
    for (t=t1; t<=t2; t++) {
	if (na(x[t])) {
	    continue;
	} else if (x[t] < 0.0) {
	    strcpy(gretl_errmsg, _("illegal negative value"));
	    *err = E_DATA;
	    break;
	} else {
	    sx[n++] = x[t];
	    sumx += x[t];
	}
    }

    if (*err) {
	free(sx);
	return NADBL;
    } else if (n == 0) {
	free(sx);
	*err = E_DATA;
	return NADBL;
    }

    qsort(sx, n, sizeof *sx, gretl_compare_doubles); 

    for (t=0; t<n; t++) {
	csx += sx[t];
	idx = t + 1;
	sisx += idx * sx[t];
	sx[t] = csx / sumx; /* now = Lorenz curve */
    }

    gini = 2.0 * sisx / (n * sumx) - ((double) n + 1) / n;

    if (plz != NULL) {
	*plz = sx;
	*pn = n;
    } else {
	free(sx);
    }

    return gini;
}

/**
 * gretl_gini:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the Gini coefficient for the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL 
 * on failure.
 */

double gretl_gini (int t1, int t2, const double *x)
{
    double g;
    int err = 0;

    g = gini_coeff(x, t1, t2, NULL, NULL, &err);

    if (err) {
	g = NADBL;
    }

    return g;
}

static int lorenz_graph (const char *vname, double *lz, int n)
{
    FILE *fp;
    double idx;
    int t, err = 0;

    if (gnuplot_init(PLOT_REGULAR, &fp)) {
	return E_FOPEN;
    }

    fputs("set key top left\n", fp);

    fprintf(fp, "set title '%s'\n", vname);
    fprintf(fp, "plot \\\n"
	    "'-' using 1:2 title '%s' w lines, \\\n"
	    "'-' using 1:2 notitle w lines\n",
	    I_("Lorenz curve"));

    gretl_push_c_numeric_locale();

    for (t=0; t<n; t++) {
	idx = t + 1;
	fprintf(fp, "%g %g\n", idx / n, lz[t]);
    }    
    fputs("e\n", fp);

    for (t=0; t<n; t++) {
	idx = ((double) t + 1) / n;
	fprintf(fp, "%g %g\n", idx, idx);
    }    
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    err = gnuplot_make_graph();

    return err;
}

/**
 * gini:
 * @vnum: ID number of variable to examine.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: not used yet.
 * @prn: gretl printing struct.
 *
 * Graphs the Lorenz curve for variable @vnum and prints the
 * Gini coefficient.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gini (int vnum, double ***pZ, DATAINFO *pdinfo, 
	  gretlopt opt, PRN *prn)
{
    double *lz = NULL;
    double gini;
    int n, fulln;
    int err = 0;

    gini = gini_coeff((*pZ)[vnum], pdinfo->t1, pdinfo->t2, 
		      &lz, &n, &err);

    if (err) {
	return err;
    }

    fulln = pdinfo->t2 - pdinfo->t1 - 1;
    pprintf(prn, "%s\n", pdinfo->varname[vnum], n);
    pprintf(prn, _("Number of observations = %d\n"), n);

    if (n < fulln) {
	pputs(prn, _("Warning: there were missing values\n"));
    } 

    pputc(prn, '\n');

    pprintf(prn, "%s = %g\n", _("Sample Gini coefficient"), gini);
    pprintf(prn, "%s = %g\n", _("Estimate of population value"), 
	    gini * (double) n / (n - 1));

    err = lorenz_graph(pdinfo->varname[vnum], lz, n);

    free(lz);

    return err;
}

