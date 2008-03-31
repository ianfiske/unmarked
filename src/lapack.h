/* 
 * Scythe Statistical Library Copyright (C) 2000-2002 Andrew D. Martin
 * and Kevin M. Quinn; 2002-present Andrew D. Martin, Kevin M. Quinn,
 * and Daniel Pemstein.  All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify under the terms of the GNU General Public License as
 * published by Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.  See the text files
 * COPYING and LICENSE, distributed with this source code, for further
 * information.
 * --------------------------------------------------------------------
 *  scythe/lapack.h
 *
 */

/*!
 * \file lapack.h
 * \brief Definitions that provide access to LAPACK/BLAS fortran
 *        routines for internal library functions.  
 *
 * This file provides function definitions that help provide
 * LAPACK/BLAS support to Scythe functions.  These definitions are not
 * part of Scythe's public interface and are used exclusively from
 * within the library.
 *
 */


#ifndef SCYTHE_LAPACK_H
#define SCYTHE_LAPACK_H

#ifdef SCYTHE_COMPILE_DIRECT
#endif

namespace scythe {

  namespace lapack {
    inline void
    make_symmetric(double* matrix, int rows)
    {
      for (int i = 1; i < rows; ++i)
        for (int j = 0; j < i; ++j)
          matrix[i * rows + j] = matrix[j * rows + i];
    }

    extern "C" {

      /* Matrix multiplication and gaxpy */
      void dgemm_ (char* transa, char* transb, const int* m,
                   const int* n, const int* k, const double* alpha,
                   const double* a, const int* lda, const double* b,
                   const int* ldb, const double* beta, double* c, 
                   const int* ldc);

      /* Matrix cross product A'A */
      void dsyrk_(const char* uplo, const char* trans, const int* n,
                  const int* k, const double* alpha, const double* a,
                  const int* lda, const double* beta, double* c,
                  const int* ldc);

      /* LU decomposition */
      void dgetrf_ (const int* rows, const int* cols, double* a,
                    const int* lda, int* ipiv, int *info);
      
      /* General inversion (given LU decomposion)*/
      void dgetri_ (const int* n, double* a, const int* lda,
                    const int* ipiv, double* work, const int* lwork,
                    int* info);

      /* Cholesky decomposition */
      void dpotrf_(const char* uplo, const int* n, double* a,
                   const int* lda, int* info);

      /* chol_solve give cholesky */
      void dpotrs_ (const char* uplo, const int* n, const int* nrhs,
                    const double* a, const int* lda, double *b,
                    const int* ldb, int* info);
      
      /* chol_solve from A and b */
      void dposv_ (const char* uplo, const int* n, const int* nrhs,
                   double* a, const int* lda, double* b, const int* ldb,
                   int* info);

      /* Positive Definite Inversion (given LU decomposition) */
      void dpotri_(const char* uplo, const int* n, double* a,
                   const int* lda, int* info);

      /* Eigenvalues/vectors for general (nonsymmetric) square matrices
      void dgeev_(const char* jobvl, const char* jobvr, const int* n,
                  double* a, const int* lda, double* wr, double* wi,
                  double* vl, const int* ldvl, double* vr, const int* ldvr,
                  double* work, const int* lwork, int* info);
                  */

      /* Eigenvalues/vectors for symmetric matrices */
      void dsyevr_ (const char* jobz, const char* range, const char* uplo,
                    const int* n, double* a, const int* lda, double* vl,
                    double* vu, const int* il, const int* iu,
                    const double* abstol, const int* m, double* w,
                    double* z, const int* ldz, int* isuppz, double*
                    work, int* lwork, int* iwork, const int* liwork,
                    int* info);

      /* QR decomposition */
      void dgeqp3_ (const int* m, const int* n, double* a, const int* lda,
                    int* jpvt, double* tau, double* work, const int* lwork,
                    int* info);

      /* QR solve routines */
      void dormqr_ (const char* side, const char* trans, const int* m,
                    const int* n, const int* k, const double* a,
                    const int* lda, const double* tau, double* c,
                    const int* ldc, double* work, const int* lwork,
                    int* info);

      void dtrtrs_ (const char* uplo, const char* trans, const char* diag,
                    const int* n, const int* nrhs, const double* a,
                    const int* lda, double* b, const int* ldb, int* info);

      /* SVD */
      void dgesdd_ (const char* jobz, const int* m, const int* n, double* a,
                    const int* lda, double* s, double* u, const int* ldu,
                    double* vt, const int* ldvt, double* work,
                    const int* lwork, int* iwork, int* info);

    } // end extern
  } // end namespace lapack
} // end namespace scythe

#endif /* SCYTHE_LAPACK_H */
