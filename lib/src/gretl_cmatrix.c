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
#include "gretl_f2c.h"
#include "clapack_complex.h"
#include "gretl_cmatrix.h"

#include <fftw3.h>
#include <complex.h>

/* helper function for fftw-based real FFT functions */

static int fft_allocate (double **px, gretl_matrix **pm,
			 fftw_complex **pc, int r, int c)
{
    *pm = gretl_matrix_alloc(r, c);
    if (*pm == NULL) {
	return E_ALLOC;
    }

    *px = fftw_malloc(r * sizeof **px);
    if (*px == NULL) {
	gretl_matrix_free(*pm);
	return E_ALLOC;
    }

    *pc = fftw_malloc((r/2 + 1 + r % 2) * sizeof **pc);
    if (*pc == NULL) {
	gretl_matrix_free(*pm);
	free(*px);
	return E_ALLOC;
    }

    return 0;
}

/* start fftw-based real FFT functions */

/**
 * gretl_matrix_fft:
 * @y: input matrix.
 * @err: location to receive error code.
 *
 * Add description here.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_fft (const gretl_matrix *y, int *err)
{
    gretl_matrix *ft = NULL;
    fftw_plan p = NULL;
    double *tmp = NULL;
    fftw_complex *out;
    int r = gretl_matrix_rows(y);
    int c, m, odd, cr, ci;
    int i, j;

    if (r < 2) {
	*err = E_DATA;
	return NULL;
    }

    c = gretl_matrix_cols(y);
    m = r / 2;
    odd = r % 2;
    cr = 0;
    ci = 1;

    *err = fft_allocate(&tmp, &ft, &out, r, 2 * c);
    if (*err) {
	return NULL;
    }

    for (j=0; j<c; j++) {
	/* load the data */
	for (i=0; i<r; i++) {
	    tmp[i] = gretl_matrix_get(y, i, j);
	}

	if (j == 0) {
	    /* make the plan just once */
	    p = fftw_plan_dft_r2c_1d(r, tmp, out, FFTW_ESTIMATE);
	}

	/* run the transform */
	fftw_execute(p);

	/* transcribe the result */
	for (i=0; i<=m+odd; i++) {
	    gretl_matrix_set(ft, i, cr, out[i][0]);
	    gretl_matrix_set(ft, i, ci, out[i][1]);
	}
	for (i=m; i>0; i--) {
	    gretl_matrix_set(ft, r-i, cr,  out[i][0]);
	    gretl_matrix_set(ft, r-i, ci, -out[i][1]);
	}
	cr += 2;
	ci += 2;
    }

    fftw_destroy_plan(p);
    fftw_free(out);
    fftw_free(tmp);

    return ft;
}

/**
 * gretl_matrix_ffti:
 * @y: input matrix.
 * @err: location to receive error code.
 *
 * Add description here.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_ffti (const gretl_matrix *y, int *err)
{
    gretl_matrix *ft = NULL;
    fftw_plan p = NULL;
    double *tmp = NULL;
    fftw_complex *in;
    int c, r = gretl_matrix_rows(y);
    int m, odd, cr, ci;
    int i, j;

    if (r < 2) {
	*err = E_DATA;
	return NULL;
    }

    c = gretl_matrix_cols(y) / 2;
    m = r / 2;
    odd = r % 2;

    if (c == 0) {
	*err = E_NONCONF;
	return NULL;
    }

    *err = fft_allocate(&tmp, &ft, &in, r, c);
    if (*err) {
	return NULL;
    }

    cr = 0;
    ci = 1;

    for (j=0; j<c; j++) {
	/* load the data */
	for (i=0; i<=m+odd; i++) {
	    in[i][0] = gretl_matrix_get(y, i, cr);
	    in[i][1] = gretl_matrix_get(y, i, ci);
	}

	if (j == 0) {
	    /* make the plan just once */
	    p = fftw_plan_dft_c2r_1d(r, in, tmp, FFTW_ESTIMATE);
	}

	/* run the transform */
	fftw_execute(p);

	/* transcribe result */
	for (i=0; i<r; i++) {
	    gretl_matrix_set(ft, i, j, tmp[i] / r);
	}
	cr += 2;
	ci += 2;
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(tmp);

    return ft;
}

/* end fftw-based real FFT functions */

/* Multiplication of complex matrices via BLAS zgemm() */

gretl_matrix *gretl_zgemm (const gretl_matrix *A, gretl_matrix *B,
			   int *err)
{
    gretl_matrix *C;
    cmplx *a, *b, *c;
    cmplx alpha = {1, 0};
    cmplx beta = {0, 0};
    char transa = 'N';
    char transb = 'N';
    integer lda, ldb, ldc;
    integer m, n, k;

    if (gretl_is_null_matrix(A) || A->rows % 2 != 0 ||
	gretl_is_null_matrix(B) || B->rows % 2 != 0) {
	*err = E_INVARG;
	return NULL;
    } else if (A->cols != B->rows / 2) {
	*err = E_NONCONF;
	return NULL;
    }	

    m = A->rows / 2;
    n = B->cols;
    k = A->cols;

    C = gretl_matrix_alloc(A->rows, B->cols);
    if (C == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    a = (cmplx *) A->val;
    b = (cmplx *) B->val;
    c = (cmplx *) C->val;

    lda = A->rows / 2;
    ldb = B->rows / 2;
    ldc = C->rows / 2;

    zgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda,
	   b, &ldb, &beta, c, &ldc);

    return C;
}

/* Eigen decomposition of complex (Hermitian) matrix using
   LAPACK's zheev() */

gretl_matrix *gretl_zheev (const gretl_matrix *A, gretl_matrix *V,
			   int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *Acpy = NULL;
    integer n, info, lwork;
    double *w = NULL;
    double *rwork = NULL;
    cmplx *a = NULL;
    cmplx *work = NULL;
    cmplx wsz;
    char jobz = V != NULL ? 'V' : 'N';
    char uplo = 'U';

    if (gretl_is_null_matrix(A) || A->rows != 2 * A->cols) {
	*err = E_INVARG;
	return NULL;
    }

    if (V != NULL && V->rows * V->cols == A->rows * A->cols) {
	/* we can use V->val as workspace */
	gretl_matrix_destroy_info(V);
	V->rows = A->rows;
	V->cols = A->cols;
	gretl_matrix_copy_values(V, A);
	a = (cmplx *) V->val;
    } else {
	/* we need a copy for workspace */
	Acpy = gretl_matrix_copy(A);
	if (Acpy == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	}
	a = (cmplx *) Acpy->val;
    }

    n = A->rows / 2;

    ret = gretl_matrix_alloc(n, 1);
    if (ret == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    w = ret->val;

    /* get optimal workspace size */
    lwork = -1;
    zheev_(&jobz, &uplo, &n, a, &n, w, &wsz, &lwork, rwork, &info);

    lwork = (integer) wsz.r;
    work = malloc(lwork * sizeof *work);
    rwork = malloc((3 * n - 2) * sizeof *rwork);
    if (work == NULL || rwork == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    /* do the actual decomposition */
    zheev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, rwork, &info);
    if (info != 0) {
	fprintf(stderr, "zheev: info = %d\n", info);
	*err = E_DATA;
    } else if (V != NULL && Acpy != NULL) {
	gretl_matrix_replace_content(V, Acpy);
    }

 bailout:

    free(rwork);
    free(work);
    if (Acpy != NULL) {
	gretl_matrix_free(Acpy);
    }

    if (*err) {
	gretl_matrix_free(ret);
	ret = NULL;
    }

    return ret;
}

/* Inverse of a complex matrix via LU decomposition using the
   LAPACK functions zgetrf() and zgetri()
*/

gretl_matrix *gretl_zgetri (const gretl_matrix *A, int *err)
{
    gretl_matrix *Ainv = NULL;
    integer lwork = -1;
    integer *ipiv;
    cmplx *a, *work = NULL;
    integer n, info;

    if (gretl_is_null_matrix(A) || A->rows != 2 * A->cols) {
	*err = E_INVARG;
	return NULL;
    }

    Ainv = gretl_matrix_copy(A);
    if (Ainv == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    n = A->cols;

    ipiv = malloc(2 * n * sizeof *ipiv);
    if (ipiv == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    a = (cmplx *) Ainv->val;

    zgetrf_(&n, &n, a, &n, ipiv, &info);
    if (info != 0) {
	printf("zgetrf: info = %d\n", info);
	*err = E_DATA;
    }

    if (!*err) {
	/* workspace size query */
	cmplx wsz;

	zgetri_(&n, a, &n, ipiv, &wsz, &lwork, &info);
	lwork = (integer) wsz.r;
	work = malloc(lwork * sizeof *work);
	if (work == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	/* actual inversion */
	zgetri_(&n, a, &n, ipiv, work, &lwork, &info);
	if (info != 0) {
	    printf("zgetri: info = %d\n", info);
	    *err = E_DATA;
	}
    }

 bailout:
    
    free(work);
    free(ipiv);

    if (*err && Ainv != NULL) {
	gretl_matrix_free(Ainv);
	Ainv = NULL;
    }

    return Ainv;
}

/* Complex FFT (or inverse) via fftw */

gretl_matrix *gretl_complex_fft (const gretl_matrix *A, int inverse,
				 int *err)
{
    gretl_matrix *B = NULL;
    fftw_complex *tmp, *ptr;
    fftw_plan p;
    int sign;
    int r, c, j;

    if (gretl_is_null_matrix(A) || A->rows % 2 != 0) {
	*err = E_INVARG;
	return NULL;
    }    

    B = gretl_matrix_copy(A);
    if (B == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    r = A->rows / 2;
    c = A->cols;

    tmp = (fftw_complex *) B->val;
    sign = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    ptr = tmp;
    for (j=0; j<c; j++) {
	p = fftw_plan_dft_1d(r, ptr, ptr, sign, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	/* advance pointer to next column */
	ptr += r;
    }

    if (inverse) {
	/* "FFTW computes an unnormalized transform: computing a
	    forward followed by a backward transform (or vice versa)
	    will result in the original data multiplied by the size of
	    the transform (the product of the dimensions)."
	    So should we do the following?
	*/
	for (j=0; j<r*c; j++) {
	    tmp[j][0] = tmp[j][0] / r;
	    tmp[j][1] = tmp[j][1] / r;
	}
    }

    return B;
}

/* Hadamard product for complex matrices that match by both
   row and column dimension, or we have a row match and one
   matrix has a single column, or we have a column match and
   one matrix has a single (complex) row.
*/

gretl_matrix *gretl_complex_hprod (const gretl_matrix *A,
				   const gretl_matrix *B,
				   int *err)
{
    gretl_matrix *C = NULL;
    const gretl_matrix *L = A;
    const gretl_matrix *R = B;
    complex double *a;
    complex double *b;
    complex double *c;
    int match = 0;
    int i, j, k;
    int cr, cc;

    if (gretl_is_null_matrix(A) || A->rows % 2 != 0 ||
	gretl_is_null_matrix(B) || B->rows % 2 != 0) {
	*err = E_INVARG;
	return NULL;
    }

    cr = A->rows;
    cc = A->cols;

    if (A->rows == B->rows && A->cols == B->cols) {
	match = 1; /* fine */
    } else if (A->rows == B->rows) {
	if (B->cols == 1) {
	    cc = A->cols;
	    match = 2;
	} else if (A->cols == 1) {
	    cc = B->cols;
	    L = B;
	    R = A;
	    match = 2;
	}
    } else if (A->cols == B->cols) {
	if (B->rows == 2) {
	    match = 3;
	} else if (A->rows == 2) {
	    cr = B->rows;
	    L = B;
	    R = A;
	    match = 3;
	}
    }

    if (match == 0) {
	*err = E_NONCONF;
	return NULL;
    }

    C = gretl_matrix_alloc(cr, cc);
    if (C == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    cr /= 2; /* rows of 16-byte values */

    a = (complex double *) L->val;
    b = (complex double *) R->val;
    c = (complex double *) C->val;

    if (match == 1) {
	int n = cr * cc;

	for (k=0; k<n; k++) {
	    c[k] = a[k] * b[k];
	}
    } else if (match == 2) {
	/* b has just one column */
	k = 0;
	for (j=0; j<cc; j++) {
	    for (i=0; i<cr; i++) {
		c[k++] = a[j*cr+i] * b[i];
	    }
	}
    } else {
	/* b has just one row */
	k = 0;
	for (j=0; j<cc; j++) {
	    for (i=0; i<cr; i++) {
		c[k++] = a[j*cr+i] * b[j];
	    }
	}
    }

    return C;
}

#define cmatrix_get_re(m,i,j) (m->val[(j)*m->rows+(i)*2])
#define cmatrix_get_im(m,i,j) (m->val[(j)*m->rows+(i)*2+1])

int complex_matrix_print (gretl_matrix *A, const char *name, PRN *prn)
{
    double re, im;
    char s[4] = "   ";
    int r, c, i, j;

    if (gretl_is_null_matrix(A) || A->rows % 2 != 0) {
	return E_INVARG;
    }

    r = A->rows / 2;
    c = A->cols;

    if (name != NULL && *name != '\0') {
	pprintf(prn, "%s (%d x %d)", name, r, c);
	pputs(prn, "\n\n");
    }

    for (i=0; i<r; i++) {
	for (j=0; j<c; j++) {
	    re = cmatrix_get_re(A, i, j);
	    im = cmatrix_get_im(A, i, j);
	    s[1] = (im >= 0) ? '+' : '-';
	    pprintf(prn, "%7.4f%s%6.4fi", re, s, fabs(im));
	    if (j < c - 1) {
		pputs(prn, "  ");
	    }
        }
        pputc(prn, '\n');
    }
    pputc(prn, '\n');

    return 0;
}

/* Compose a complex matrix from its real and imaginary
   components */

gretl_matrix *gretl_cmatrix (const gretl_matrix *Re,
			     const gretl_matrix *Im,
			     int *err)
{
    gretl_matrix *C = NULL;
    int i, ic, n;

    if (gretl_is_null_matrix(Re) ||
	gretl_is_null_matrix(Im) ||
	Re->rows != Im->rows ||
	Re->cols != Im->cols) {
	*err = E_INVARG;
	return NULL;
    }

    n = Re->rows * Re->cols;

    C = gretl_matrix_alloc(2 * Re->rows, Re->cols);
    if (C == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    ic = 0;
    for (i=0; i<n; i++) {
	C->val[ic] = Re->val[i];
	C->val[ic+1] = Im->val[i];
	ic += 2;
    }

    return C;
}

/* Extract the real part of @A (if im==0) or the imaginary part
   (if im==1)
*/

gretl_matrix *gretl_cxtract (const gretl_matrix *A, int im,
			     int *err)
{
    gretl_matrix *C = NULL;
    int i, j, r, n;

    if (gretl_is_null_matrix(A)) {
	*err = E_INVARG;
	return NULL;
    }

    r = A->rows / 2;
    n = r * A->cols;

    C = gretl_matrix_alloc(r, A->cols);
    if (C == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    j = im != 0;
    for (i=0; i<n; i++) {
	C->val[i] = A->val[j];
	j += 2;
    }

    return C;
}

/* Return conjugate transpose of complex matrix @A */

gretl_matrix *gretl_ctran (const gretl_matrix *A, int *err)
{
    gretl_matrix *C = NULL;
    cmplx *a, *c, aij;
    int i, j, ra, rc;

    if (gretl_is_null_matrix(A) || A->rows % 2 != 0) {
	*err = E_INVARG;
	return NULL;
    }

    C = gretl_matrix_alloc(2 * A->cols, A->rows / 2);
    if (C == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    a = (cmplx*) A->val;
    c = (cmplx*) C->val;
    ra = A->rows / 2;
    rc = C->rows / 2;

    for (j=0; j<A->cols; j++) {
	for (i=0; i<ra; i++) {
	    aij = a[j*ra+i];
	    aij.i = -aij.i;
	    c[i*rc+j] = aij;
	}
    }

    return C;
}
