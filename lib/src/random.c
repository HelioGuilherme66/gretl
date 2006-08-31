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

/* random.c for gretl */

#include "libgretl.h"

#include <time.h>

#ifdef USE_GLIB2
# include <glib.h>
#else
# include "mt19937ar.c"
#endif

#ifdef USE_GLIB2
static GRand *gretl_rand;
#endif

static unsigned int useed;

unsigned int get_gretl_random_seed (void)
{
    return useed;
}

/**
 * gretl_rand_init:
 *
 * Initialize gretl's PRNG, using the system time as seed.
 */

void gretl_rand_init (void)
{
    useed = time(NULL);
#ifdef USE_GLIB2
    gretl_rand = g_rand_new();
    gretl_rand_set_seed((guint32) useed);
#else
    init_genrand(useed);
#endif
}

/**
 * gretl_rand_set_seed:
 * @seed: the chosen seed value.
 *
 * Set a specific (and hence reproducible) seed for gretl's PRNG.
 * But if the value 0 is given for @seed, set the seed using
 * the system time (which is the default when libgretl is
 * initialized).
 */

void gretl_rand_set_seed (unsigned int seed)
{
    if (seed == 0) {
	useed = time(NULL);
    } else {
	useed = seed;
    }
#ifdef USE_GLIB2
    g_rand_set_seed(gretl_rand, useed);
#else
    init_genrand(useed);
#endif
}

/**
 * gretl_uniform_dist_minmax:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 * @min: lower bound of range.
 * @min: upper bound of range.
 *
 * Fill the selected subset of array @a with pseudo-random drawings
 * from the uniform distribution on @min-@max, using the Mersenne
 * Twister.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_uniform_dist_minmax (double *a, int t1, int t2,
			       double min, double max) 
{
    int t;

    if (na(min) && na(max)) {
	min = 0.0;
	max = 1.0;
    } else if (na(min) || na(max) || max <= min) {
	return E_INVARG;
    }

    for (t=t1; t<=t2; t++) {
#ifdef USE_GLIB2
	a[t] = g_rand_double_range(gretl_rand, min, max);
#else
	a[t] = min + genrand_int32() * (max / 4294967296.0);
#endif 
    }

    return 0;
}

/**
 * gretl_uniform_dist:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the uniform distribution on 0-1, using the Mersenne
 * Twister.
 */

void gretl_uniform_dist (double *a, int t1, int t2) 
{
    int t;

    for (t=t1; t<=t2; t++) {
#ifdef USE_GLIB2
	a[t] = g_rand_double_range(gretl_rand, 0.0, 1.0);
#else
	a[t] = genrand_int32() * (1.0 / 4294967296.0);
#endif 
    }
}

/**
 * gretl_normal_dist:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the standard normal distribution, using the Mersenne Twister
 * for uniform input and the Box-Muller method for converting to the
 * normal distribution.
 */

void gretl_normal_dist (double *a, int t1, int t2) 
{
    double x, y, z;
    int t;

    for (t=t1; t<=t2; t++) {
    tryagain:
#ifdef USE_GLIB2
	x = g_rand_double(gretl_rand);
	y = g_rand_double(gretl_rand);
#else
	x = genrand_real2();
	y = genrand_real2();
#endif
	z = sqrt(-2. * log(x));
	if (isnan(z) || isinf(z)) {
	    goto tryagain;
	}
	a[t] = z * cos(2. * M_PI * y);
    }
}

/**
 * gretl_normal_dist_with_params:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 * @mean: mean of the distribution
 * @sd: standard deviation
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the normal distribution with the given mean and standard
 * deviation, using the Mersenne Twister for uniform input and the 
 * Box-Muller method for converting to the normal distribution.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int
gretl_normal_dist_with_params (double *a, int t1, int t2,
			       double mean, double sd) 
{
    int t;

    if (na(mean) && na(sd)) {
	mean = 0.0;
	sd = 1.0;
    } else if (na(mean) || na(sd) || sd <= 0.0) {
	return E_INVARG;
    }

    gretl_normal_dist(a, t1, t2);

    if (mean != 0.0 || sd != 1.0) {
	for (t=t1; t<=t2; t++) {
	    a[t] = mean + a[t] * sd;
	}
    }
    
    return 0;
}

/**
 * gretl_chisq_dist:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @v: degrees of freedom.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Chi-Squared distribution with @v degrees of freedom, 
 * using the Mersenne Twister for uniform input and the Box-Muller 
 * method for converting to the normal distribution.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_chisq_dist (double *a, int t1, int t2, int v) 
{
    double x, y, z;
    int i, t;

    if (v < 1) {
	return E_INVARG;
    }

    for (t=t1; t<=t2; t++) {
	a[t] = 0.0;
	for (i=0; i<v; i++) {
	tryagain:
#ifdef USE_GLIB2
	    x = g_rand_double(gretl_rand);
	    y = g_rand_double(gretl_rand);
#else
	    x = genrand_real2();
	    y = genrand_real2();
#endif
	    z = sqrt(-2. * log(x));
	    if (isnan(z) || isinf(z)) {
		goto tryagain;
	    }
	    z *= cos(2. * M_PI * y);
	    a[t] += z * z;
	}
    }

    return 0;
}

/**
 * gretl_t_dist:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @v: degrees of freedom.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Student t distribution with @v degrees of freedom, 
 * using the Mersenne Twister for uniform input and the Box-Muller 
 * method for converting to the normal distribution.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_t_dist (double *a, int t1, int t2, int v) 
{
    double *X2 = NULL;
    int T = t2 - t1 + 1;
    int t;

    if (v < 1) {
	return E_INVARG;
    }

    X2 = malloc(T * sizeof *X2);
    if (X2 == NULL) {
	return E_ALLOC;
    }

    gretl_normal_dist(a, t1, t2);
    gretl_chisq_dist(X2, 0, T-1, v);

    for (t=0; t<T; t++) {
	a[t + t1] /= sqrt(X2[t] / v);
    }

    free(X2);

    return 0;
}

/**
 * gretl_rand_int_max:
 * @max: the maximum value (open)
 *
 * Returns: a pseudo-random unsigned int in the interval [0, @max)
 * using the Mersenne Twister.
 */

unsigned int gretl_rand_int_max (unsigned int max)
{
#ifdef USE_GLIB2
    return g_rand_int_range(gretl_rand, 0, max);
#else
    return genrand_int32() * (max / 4294967296.0);
#endif
}

/**
 * gretl_rand_int:
 *
 * Returns: a pseudo-random unsigned int on the interval
 * [0,0xffffffff] using the Mersenne Twister.
 */

unsigned int gretl_rand_int (void)
{
#ifdef USE_GLIB2
    return g_rand_int(gretl_rand);
#else
    return genrand_int32();
#endif
}

/**
 * gretl_rand_free:
 *
 * Free the gretl_rand structure (may be called at program exit).
 */

void gretl_rand_free (void)
{
#ifdef USE_GLIB2
    g_rand_free(gretl_rand);
#else
    return;
#endif
}


