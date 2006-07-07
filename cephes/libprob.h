#ifndef LIBPROB_H
#define LIBPROB_H

/* area under the binomial pdf, from 0 to k, of the
 * binomial distribution with success probability p
 * and n trials.
*/

double bdtr (int k, int n, double p);

/* area under the binomial pdf, from k+1 to n, of the
 * binomial distribution with success probability p
 * and n trials. 
*/

double bdtrc (int k, int n, double p);

/* area under the left hand tail (from 0 to x)
 * of the Chi square probability density function with
 * v degrees of freedom.
 */

double chdtr (double v, double x);

/* area under the right hand tail (from x to infinity)
 * of the Chi square probability density function with
 * v degrees of freedom.
 */

double chdtrc (double v, double x);

/*
 * Finds the Chi-square argument x such that the integral
 * from x to infinity of the Chi-square density is equal
 * to the given cumulative probability y.
 */

double chdtri (double df, double y);

/*
 * Returns the area from x to infinity under the F density
 * function (also known as Snedcor's density or the
 * variance ratio density).
 */

double fdtrc (int ia, int ib, double x);

/*
 * Finds the F density argument x such that the integral
 * from x to infinity of the F density is equal to the
 * given probability p.
 */

double fdtri (int ia, int ib, double y);

/*
 * Computes the integral from minus infinity to t of the Student
 * t distribution with integer k > 0 degrees of freedom.
 */

double stdtr (int k, double t);

/*
 * Given probability p, finds the argument t such that stdtr(k,t)
 * is equal to p.
 */

double stdtri (int k, double p);

/*
 * Returns the area under the Gaussian probability density
 * function, integrated from minus infinity to x.
 */

double ndtr (double a);

/*
 * Returns the argument, x, for which the area under the
 * Gaussian probability density function (integrated from
 * minus infinity to x) is equal to y.
 */

double ndtri (double y0);

/*
 * Returns gamma function of the argument.  The result is
 * correctly signed.
 */

double cephes_gamma (double x); /* cephes' gamma(), renamed */

/*
 * Returns the base e (2.718...) logarithm of the absolute
 * value of the gamma function of the argument.
 * The sign (+1 or -1) of the gamma function is set in a
 * global variable named cephes_sgngam.
 */

double cephes_lgamma (double x); /* alias for cephes' lgam() */

/*
 * Returns the current value of cephes_sgngam.
 */

int get_cephes_sgngam (void);

/* Accessor for cephes error code */

int get_cephes_errno (void);

#endif /* LIBPROB_H */





