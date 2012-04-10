#include "minpack.h"
#include <math.h>

/* given an m by n matrix a, an n by n nonsingular diagonal
   matrix d, an m-vector b, and a positive number delta,
   the problem is to determine a value for the parameter
   par such that if x solves the system

   a*x = b ,     sqrt(par)*d*x = 0 ,

   in the least squares sense, and dxnorm is the euclidean
   norm of d*x, then either par is zero and

   (dxnorm-delta) .le. 0.1*delta ,

   or par is positive and

   abs(dxnorm-delta) .le. 0.1*delta .

   this subroutine completes the solution of the problem
   if it is provided with the necessary information from the
   qr factorization, with column pivoting, of a. that is, if
   a*p = q*r, where p is a permutation matrix, q has orthogonal
   columns, and r is an upper triangular matrix with diagonal
   elements of nonincreasing magnitude, then lmpar expects
   the full upper triangle of r, the permutation matrix p,
   and the first n components of (q transpose)*b. on output
   lmpar also provides an upper triangular matrix s such that

   t   t                   t
   p *(a *a + par*d*d)*p = s *s .

   s is employed within lmpar and may be of separate interest.

   only a few iterations are generally needed for convergence
   of the algorithm. if, however, the limit of 10 iterations
   is reached, then the output par will contain the best
   value obtained so far.

   n is a positive integer input variable set to the order of r.

   r is an n by n array. on input the full upper triangle
   must contain the full upper triangle of the matrix r.
   on output the full upper triangle is unaltered, and the
   strict lower triangle contains the strict upper triangle
   (transposed) of the upper triangular matrix s.

   ldr is a positive integer input variable not less than n
   which specifies the leading dimension of the array r.

   ipvt is an integer input array of length n which defines the
   permutation matrix p such that a*p = q*r. column j of p
   is column ipvt(j) of the identity matrix.

   diag is an input array of length n which must contain the
   diagonal elements of the matrix d.

   qtb is an input array of length n which must contain the first
   n elements of the vector (q transpose)*b.

   delta is a positive input variable which specifies an upper
   bound on the euclidean norm of d*x.

   par is a nonnegative variable. on input par contains an
   initial estimate of the levenberg-marquardt parameter.
   on output par contains the final estimate.

   x is an output array of length n which contains the least
   squares solution of the system a*x = b, sqrt(par)*d*x = 0,
   for the output par.

   sdiag is an output array of length n which contains the
   diagonal elements of the upper triangular matrix s.

   wa1 and wa2 are work arrays of length n.

   argonne national laboratory. minpack project. march 1980.
   burton s. garbow, kenneth e. hillstrom, jorge j. more
*/

int lmpar_(int n, double *r, int ldr, 
	   int *ipvt, double *diag, double *qtb, double *delta, 
	   double *par, double *x, double *sdiag, double *wa1, 
	   double *wa2)
{
    const double p1 = .1;
    const double p001 = .001;

    int r_dim1, r_offset;
    double d1, d2;

    int i, j, k, l;
    double fp;
    int jm1, jp1;
    double sum, parc, parl;
    int iter;
    double temp, paru, dwarf;
    int nsing;
    double gnorm;
    double dxnorm;

    --wa2;
    --wa1;
    --sdiag;
    --x;
    --qtb;
    --diag;
    --ipvt;
    r_dim1 = ldr;
    r_offset = 1 + r_dim1;
    r -= r_offset;

    /* dwarf is the smallest positive magnitude */
    dwarf = dpmpar_(2);

    /* compute and store in x the gauss-newton direction. if the */
    /* jacobian is rank-deficient, obtain a least squares solution. */

    nsing = n;
    for (j = 1; j <= n; ++j) {
	wa1[j] = qtb[j];
	if (r[j + j * r_dim1] == 0.0 && nsing == n) {
	    nsing = j - 1;
	}
	if (nsing < n) {
	    wa1[j] = 0.0;
	}
    }

    if (nsing >= 1) {
	for (k = 1; k <= nsing; ++k) {
	    j = nsing - k + 1;
	    wa1[j] /= r[j + j * r_dim1];
	    temp = wa1[j];
	    jm1 = j - 1;
	    if (jm1 >= 1) {
		for (i = 1; i <= jm1; ++i) {
		    wa1[i] -= r[i + j * r_dim1] * temp;
		}
	    }
	}
    }

    for (j = 1; j <= n; ++j) {
	l = ipvt[j];
	x[l] = wa1[j];
    }

    /* initialize the iteration counter. */
    /* evaluate the function at the origin, and test */
    /* for acceptance of the gauss-newton direction. */

    iter = 0;
    for (j = 1; j <= n; ++j) {
	wa2[j] = diag[j] * x[j];
    }
    dxnorm = enorm_(n, &wa2[1]);
    fp = dxnorm - *delta;
    if (fp <= p1 * *delta) {
	*par = 0.0;
	return 0;
    }

    /* if the jacobian is not rank deficient, the newton */
    /* step provides a lower bound, parl, for the zero of */
    /* the function. otherwise set this bound to zero. */

    parl = 0.0;
    if (nsing >= n) {
	for (j = 1; j <= n; ++j) {
	    l = ipvt[j];
	    wa1[j] = diag[l] * (wa2[l] / dxnorm);
	}
	for (j = 1; j <= n; ++j) {
	    sum = 0.0;
	    jm1 = j - 1;
	    if (jm1 >= 1) {
		for (i = 1; i <= jm1; ++i) {
		    sum += r[i + j * r_dim1] * wa1[i];
		}
	    }
	    wa1[j] = (wa1[j] - sum) / r[j + j * r_dim1];
	}
	temp = enorm_(n, &wa1[1]);
	parl = fp / *delta / temp / temp;
    }

    /* calculate an upper bound, paru, for the zero of the function. */

    for (j = 1; j <= n; ++j) {
	sum = 0.0;
	for (i = 1; i <= j; ++i) {
	    sum += r[i + j * r_dim1] * qtb[i];
	}
	l = ipvt[j];
	wa1[j] = sum / diag[l];
    }
    gnorm = enorm_(n, &wa1[1]);
    paru = gnorm / *delta;
    if (paru == 0.0) {
	paru = dwarf / min(*delta,p1);
    }

    /* if the input par lies outside of the interval (parl, paru), */
    /* set par to the closer endpoint. */

    *par = max(*par, parl);
    *par = min(*par, paru);
    if (*par == 0.0) {
	*par = gnorm / dxnorm;
    }

    /* beginning of iteration */

    while (1) {
	++iter;

	/* evaluate the function at the current value of par */

	if (*par == 0.0) {
	    d1 = dwarf, d2 = p001 * paru;
	    *par = max(d1, d2);
	}
	temp = sqrt(*par);
	for (j = 1; j <= n; ++j) {
	    wa1[j] = temp * diag[j];
	}

	qrsolv_(n, &r[r_offset], ldr, &ipvt[1], &wa1[1], &qtb[1], &x[1], 
		&sdiag[1], &wa2[1]);

	for (j = 1; j <= n; ++j) {
	    wa2[j] = diag[j] * x[j];
	}

	dxnorm = enorm_(n, &wa2[1]);
	temp = fp;
	fp = dxnorm - *delta;

	/* If the function is small enough, accept the current value
	   of par. Also test for the exceptional cases where parl
	   is zero or the number of iterations has reached 10. 
	*/

	if (fabs(fp) <= p1 * *delta || 
	    (parl == 0.0 && fp <= temp && temp < 0.0) ||
	    iter == 10) {
	    break;
	}

	/* compute the newton correction */

	for (j = 1; j <= n; ++j) {
	    l = ipvt[j];
	    wa1[j] = diag[l] * (wa2[l] / dxnorm);
	}
	for (j = 1; j <= n; ++j) {
	    wa1[j] /= sdiag[j];
	    temp = wa1[j];
	    jp1 = j + 1;
	    if (n >= jp1) {
		for (i = jp1; i <= n; ++i) {
		    wa1[i] -= r[i + j * r_dim1] * temp;
		}
	    }
	}

	temp = enorm_(n, &wa1[1]);
	parc = fp / *delta / temp / temp;

	/* depending on the sign of the function, update parl or paru */
	if (fp > 0.0) {
	    parl = max(parl, *par);
	}
	if (fp < 0.0) {
	    paru = min(paru, *par);
	}

	/* compute an improved estimate for par */
	d1 = parl, d2 = *par + parc;
	*par = max(d1, d2);
    }

    if (iter == 0) {
	*par = 0.0;
    }

    return 0;
}

