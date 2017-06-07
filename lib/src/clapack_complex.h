#ifndef CLAPACK_COMPLEX_H
#define CLAPACK_COMPLEX_H

/* LAPACK subroutines: double-precision complex versions only */

void zheev_ (const char *jobz, const char *uplo, integer *n,
	     cmplx *a, integer *lda, double *w, cmplx *work,
	     integer *lwork, double *rwork, integer *info);

void zgetrf_ (integer *m, integer *n, cmplx *a, integer *lda,
	      integer *ipiv, integer *info);

void zgetri_ (integer *n, cmplx *a, integer *lda, integer *ipiv,
	      cmplx *work, integer *lwork, integer *info);

void zgemm_ (const char *transa, const char *transb,
	     integer *m, integer *n, integer *k,
	     cmplx *alpha, cmplx *a, integer *lda,
	     cmplx *b, integer *ldb, cmplx *beta,
	     cmplx *c, integer *ldc);

#endif /* CLAPACK_COMPLEX_H */
