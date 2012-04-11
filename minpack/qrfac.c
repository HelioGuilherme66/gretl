#include "minpack.h"
#include <math.h>
#include <float.h>

/*
c     subroutine qrfac
c
c     this subroutine uses householder transformations with column
c     pivoting (optional) to compute a qr factorization of the
c     m by n matrix a. that is, qrfac determines an orthogonal
c     matrix q, a permutation matrix p, and an upper trapezoidal
c     matrix r with diagonal elements of nonincreasing magnitude,
c     such that a*p = q*r. the householder transformation for
c     column k, k = 1,2,...,min(m,n), is of the form
c
c                           t
c           i - (1/u(k))*u*u
c
c     where u has zeros in the first k-1 positions. the form of
c     this transformation and the method of pivoting first
c     appeared in the corresponding linpack subroutine.
c
c     the subroutine statement is
c
c       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of a.
c
c       n is a positive integer input variable set to the number
c         of columns of a.
c
c       a is an m by n array. on input a contains the matrix for
c         which the qr factorization is to be computed. on output
c         the strict upper trapezoidal part of a contains the strict
c         upper trapezoidal part of r, and the lower trapezoidal
c         part of a contains a factored form of q (the non-trivial
c         elements of the u vectors described above).
c
c       lda is a positive integer input variable not less than m
c         which specifies the leading dimension of the array a.
c
c       pivot is a logical input variable. if pivot is set true,
c         then column pivoting is enforced. if pivot is set false,
c         then no column pivoting is done.
c
c       ipvt is an integer output array of length lipvt. ipvt
c         defines the permutation matrix p such that a*p = q*r.
c         column j of p is column ipvt(j) of the identity matrix.
c         if pivot is false, ipvt is not referenced.
c
c       lipvt is a positive integer input variable. if pivot is false,
c         then lipvt may be as small as 1. if pivot is true, then
c         lipvt must be at least n.
c
c       rdiag is an output array of length n which contains the
c         diagonal elements of r.
c
c       acnorm is an output array of length n which contains the
c         norms of the corresponding columns of the input matrix a.
c         if this information is not needed, then acnorm can coincide
c         with rdiag.
c
c       wa is a work array of length n. if pivot is false, then wa
c         can coincide with rdiag.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar,enorm
c
c       fortran-supplied ... dmax1,dsqrt,min0
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
*/

int qrfac_(int m, int n, double *a, int lda, 
	   int *ipvt, double *rdiag, double *acnorm, 
	   double *wa)
{
    const double p05 = .05;

    int a_dim1, a_offset, i2, i3;
    double d1, d2, d3;

    int i, j, k, jp1;
    double sum;
    int kmax;
    double temp;
    int minmn;
    double epsmch;
    double ajnorm;

    --wa;
    --acnorm;
    --rdiag;
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;

    epsmch = DBL_EPSILON;

    /* compute the initial column norms and initialize several arrays */
    for (j = 1; j <= n; ++j) {
	acnorm[j] = enorm_(m, &a[j * a_dim1 + 1]);
	rdiag[j] = acnorm[j];
	wa[j] = rdiag[j];
	ipvt[j] = j;
    }

    /* reduce a to r with Householder transformations. */

    minmn = min(m,n);
    for (j = 1; j <= minmn; ++j) {
	/* bring the column of largest norm into the pivot position */
	kmax = j;
	for (k = j; k <= n; ++k) {
	    if (rdiag[k] > rdiag[kmax]) {
		kmax = k;
	    }
	}
	if (kmax != j) {
	    for (i = 1; i <= m; ++i) {
		temp = a[i + j * a_dim1];
		a[i + j * a_dim1] = a[i + kmax * a_dim1];
		a[i + kmax * a_dim1] = temp;
	    }
	    rdiag[kmax] = rdiag[j];
	    wa[kmax] = wa[j];
	    k = ipvt[j];
	    ipvt[j] = ipvt[kmax];
	    ipvt[kmax] = k;
	}

	/* compute the householder transformation to reduce the
	   j-th column of a to a multiple of the j-th unit vector 
	*/

	i2 = m - j + 1;
	ajnorm = enorm_(i2, &a[j + j * a_dim1]);
	if (ajnorm == 0.0) {
	    goto L100;
	}
	if (a[j + j * a_dim1] < 0.0) {
	    ajnorm = -ajnorm;
	}
	for (i = j; i <= m; ++i) {
	    a[i + j * a_dim1] /= ajnorm;
	}
	a[j + j * a_dim1] += 1.0;

	/* apply the transformation to the remaining columns
	   and update the norms 
	*/

	jp1 = j + 1;
	if (n < jp1) {
	    goto L100;
	}
	for (k = jp1; k <= n; ++k) {
	    sum = 0.0;
	    for (i = j; i <= m; ++i) {
		sum += a[i + j * a_dim1] * a[i + k * a_dim1];
	    }
	    temp = sum / a[j + j * a_dim1];
	    for (i = j; i <= m; ++i) {
		a[i + k * a_dim1] -= temp * a[i + j * a_dim1];
	    }
	    if (rdiag[k] != 0.0) {
		temp = a[j + k * a_dim1] / rdiag[k];
		d3 = temp;
		d1 = 0.0, d2 = 1.0 - d3 * d3;
		rdiag[k] *= sqrt((max(d1,d2)));
		d1 = rdiag[k] / wa[k];
		if (p05 * (d1 * d1) <= epsmch) {
		    i3 = m - j;
		    rdiag[k] = enorm_(i3, &a[jp1 + k * a_dim1]);
		    wa[k] = rdiag[k];
		}
	    }
	}
    L100:
	rdiag[j] = -ajnorm;
    }

    return 0;
}

