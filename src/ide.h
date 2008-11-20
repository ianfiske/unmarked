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
 *  scythestat/ide.h
 *
 *
 */
 
/*!  \file ide.h
 *
 * \brief Definitions for inversion and decomposition functions that
 * operate on Scythe's Matrix objects.
 *
 * This file provides a number of common inversion and decomposition
 * routines that operate on Matrix objects.  It also provides related
 * functions for solving linear systems of equations and calculating
 * the determinant of a Matrix.
 *
 * Scythe will use LAPACK/BLAS routines to perform these operations on
 * concrete column-major matrices of double-precision floating point
 * numbers if LAPACK/BLAS is available and you compile your program
 * with the SCYTHE_LAPACK flag enabled.
 *
 * \note As is the case throughout the library, we provide both
 * general and default template definitions of the Matrix-returning
 * functions in this file, explicitly providing documentation for only
 * the general template versions. As is also often the case, Doxygen
 * does not always correctly add the default template definition to
 * the function list below; there is always a default template
 * definition available for every function.
 */

/* TODO: This interface exposes the user to too much implementation.
 * We need a solve function and a solver object.  By default, solve
 * would run lu_solve and the solver factory would return lu_solvers
 * (or perhaps a solver object encapsulating an lu_solver).  Users
 * could choose cholesky when appropriate.  Down the road, qr or svd
 * would become the default and we'd be able to handle non-square
 * matrices.  Instead of doing an lu_decomp or a cholesky and keeping
 * track of the results to repeatedly solve for different b's with A
 * fixed in Ax=b, you'd just call the operator() on your solver object
 * over and over, passing the new b each time.  No decomposition
 * specific solvers (except as toggles to the solver object and
 * solve function).  We'd still provide cholesky and lu_decomp.  We
 * could also think about a similar approach to inversion (one
 * inversion function with an option for method).
 *
 * If virtual dispatch in C++ wasn't such a performance killer (no
 * compiler optimization across virtual calls!!!) there would be an
 * obvious implementation of this interface using simple polymorphism.
 * Unfortunately, we need compile-time typing to maintain performance
 * and makes developing a clean interface that doesn't force users to
 * be template wizards much harder.  Initial experiments with the
 * Barton and Nackman trick were ugly.  The engine approach might work
 * a bit better but has its problems too.  This is not going to get
 * done for the 1.0 release, but it is something we should come back
 * to.
 *
 */

#ifndef SCYTHE_IDE_H
#define SCYTHE_IDE_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#include "error.h"
#include "defs.h"
#ifdef SCYTHE_LAPACK
#include "lapack.h"
#include "stat.h"
#endif
#else
#include "scythestat/matrix.h"
#include "scythestat/error.h"
#include "scythestat/defs.h"
#ifdef SCYTHE_LAPACK
#include "scythestat/lapack.h"
#include "scythestat/stat.h"
#endif
#endif

#include <cmath>
#include <algorithm>

namespace scythe {

  namespace {
    typedef unsigned int uint;
  }
  
  /*! 
   * \brief Cholesky decomposition of a symmetric positive-definite
   * matrix.
   *
   * This function performs Cholesky decomposition.  That is, given a
   * symmetric positive definite Matrix, \f$A\f$, cholesky() returns a
   * lower triangular Matrix \f$L\f$ such that \f$A = LL^T\f$.  This
   * function is faster than lu_decomp() and, therefore, preferable in
   * cases where one's Matrix is symmetric positive definite.
   *
   * \param A The symmetric positive definite Matrix to decompose.
   *
   * \see chol_solve(const Matrix<T,PO1,PS1> &, const Matrix<T,PO2,PS2> &)
   * \see chol_solve(const Matrix<T,PO1,PS1> &, const Matrix<T,PO2,PS2> &, const Matrix<T,PO3,PS3> &)
   * \see lu_decomp(Matrix<T,PO1,PS1>, Matrix<T,PO2,Concrete>&, Matrix<T,PO3,Concrete>&, Matrix<unsigned int, PO4, Concrete>&)
   *
   * \throw scythe_alloc_error (Level 1)
   * \throw scythe_dimension_error (Level 1)
   * \throw scythe_null_error (Level 1)
   * \throw scythe_type_error (Level 2)
   * \throw scythe_alloc_error (Level 1)
   *
   */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  cholesky (const Matrix<T, PO, PS>& A)
  {
    SCYTHE_CHECK_10(! A.isSquare(), scythe_dimension_error,
        "Matrix not square");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "Matrix is NULL");
    // Rounding errors can make this problematic.  Leaving out for now
    //SCYTHE_CHECK_20(! A.isSymmetric(), scythe_type_error,
    //    "Matrix not symmetric");
    
    Matrix<T,RO,Concrete> temp (A.rows(), A.cols(), false);
    T h;
    
    if (PO == Row) { // row-major optimized
      for (uint i = 0; i < A.rows(); ++i) {
        for (uint j = i; j < A.cols(); ++j) {
          h = A(i,j);
          for (uint k = 0; k < i; ++k)
            h -= temp(i, k) * temp(j, k);
          if (i == j) {
            SCYTHE_CHECK_20(h <= (T) 0, scythe_type_error,
                "Matrix not positive definite");

            temp(i,i) = std::sqrt(h);
          } else {
            temp(j,i) = (((T) 1) / temp(i,i)) * h;
            temp(i,j) = (T) 0;
          }
        }
      }
    } else { // col-major optimized
      for (uint j = 0; j < A.cols(); ++j) {
        for (uint i = j; i < A.rows(); ++i) {
          h = A(i, j);
          for (uint k = 0; k < j; ++k)
            h -= temp(j, k) * temp(i, k);
          if (i == j) {
            SCYTHE_CHECK_20(h <= (T) 0, scythe_type_error,
                "Matrix not positive definite");
            temp(j,j) = std::sqrt(h);
          } else {
            temp(i,j) = (((T) 1) / temp(j,j)) * h;
            temp(j,i) = (T) 0;
          }
        }
      }
    }
  
    SCYTHE_VIEW_RETURN(T, RO, RS, temp)
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  cholesky (const Matrix<T,O,S>& A)
  {
    return cholesky<O,Concrete>(A);
  }

  namespace {
    /* This internal routine encapsulates the  
     * algorithm used within chol_solve and lu_solve.  
     */
    template <typename T,
              matrix_order PO1, matrix_style PS1,
              matrix_order PO2, matrix_style PS2,
              matrix_order PO3, matrix_style PS3>
    inline void
    solve(const Matrix<T,PO1,PS1>& L, const Matrix<T,PO2,PS2>& U,
          Matrix<T,PO3,PS3> b, T* x, T* y)
    {
      T sum;

      /* TODO: Consider optimizing for ordering.  Experimentation
       * shows performance gains are probably minor (compared col-major
       * with and without lapack solve routines).
       */
      // solve M*y = b
      for (uint i = 0; i < b.size(); ++i) {
        sum = T (0);
        for (uint j = 0; j < i; ++j) {
          sum += L(i,j) * y[j];
        }
        y[i] = (b[i] - sum) / L(i, i);
      }
          
      // solve M'*x = y
      if (U.isNull()) { // A= LL^T
        for (int i = b.size() - 1; i >= 0; --i) {
          sum = T(0);
          for (uint j = i + 1; j < b.size(); ++j) {
            sum += L(j,i) * x[j];
          }
          x[i] = (y[i] - sum) / L(i, i);
        }
      } else { // A = LU
        for (int i = b.size() - 1; i >= 0; --i) {
          sum = T(0);
          for (uint j = i + 1; j < b.size(); ++j) {
            sum += U(i,j) * x[j];
          }
          x[i] = (y[i] - sum) / U(i, i);
        }
      }
    }
  }

 /*!\brief Solve \f$Ax=b\f$ for x via backward substitution, given a
  * lower triangular matrix resulting from Cholesky decomposition
  *
  * This function solves the system of equations \f$Ax = b\f$ via
  * backward substitution. \a L is the lower triangular matrix generated
  * by Cholesky decomposition such that \f$A = LL'\f$.
  *
  * This function is intended for repeatedly solving systems of
  * equations based on \a A.  That is \a A stays constant while \a
  * b varies.
  *
  * \param A A symmetric positive definite Matrix.
  * \param b A column vector with as many rows as \a A.
  * \param M The lower triangular matrix from the Cholesky decomposition of \a A.
  *
  * \see chol_solve(const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&)
  * \see cholesky(const Matrix<T, PO, PS>&)
  * \see lu_solve (const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&, const Matrix<T,PO3,PS3>&, const Matrix<T,PO4,PS4>&, const Matrix<unsigned int, PO5, PS5>&)
  * \see lu_solve (Matrix<T,PO1,PS1>, const Matrix<T,PO2,PS2>&)
  *
  * \throw scythe_alloc_error (Level 1)
  * \throw scythe_null_error (Level 1)
  * \throw scythe_dimension_error (Level 1)
  * \throw scythe_conformation_error (Level 1)
  *
  */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2,
            matrix_order PO3, matrix_style PS3>
  Matrix<T,RO,RS>
  chol_solve (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& b,
              const Matrix<T,PO3,PS3>& M)
  {
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "A is NULL")
    SCYTHE_CHECK_10(! b.isColVector(), scythe_dimension_error,
        "b must be a column vector");
    SCYTHE_CHECK_10(A.rows() != b.rows(), scythe_conformation_error,
        "A and b do not conform");
    SCYTHE_CHECK_10(A.rows() != M.rows(), scythe_conformation_error,
        "A and M do not conform");
    SCYTHE_CHECK_10(! M.isSquare(), scythe_dimension_error,
        "M must be square");

    T *y = new T[A.rows()];
    T *x = new T[A.rows()];
    
    solve(M, Matrix<>(), b, x, y);

    Matrix<T,RO,RS> result(A.rows(), 1, x);
     
    delete[]x;
    delete[]y;
   
    return result;
  }

  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2,
            matrix_order PO3, matrix_style PS3>
  Matrix<T,PO1,Concrete>
  chol_solve (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& b,
              const Matrix<T,PO3,PS3>& M)
  {
    return chol_solve<PO1,Concrete>(A,b,M);
  }

 /*!\brief Solve \f$Ax=b\f$ for x via backward substitution, 
  * using Cholesky decomposition
  *
  * This function solves the system of equations \f$Ax = b\f$ via
  * backward substitution and Cholesky decomposition. \a A must be a
  * symmetric positive definite matrix for this method to work. This
  * function calls cholesky() to perform the decomposition.
  *
  * \param A A symmetric positive definite matrix.
  * \param b A column vector with as many rows as \a A.
  *
  * \see chol_solve(const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&, const Matrix<T,PO3,PS3>&)
  * \see cholesky(const Matrix<T, PO, PS>&)
  * \see lu_solve (const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&, const Matrix<T,PO3,PS3>&, const Matrix<T,PO4,PS4>&, const Matrix<unsigned int, PO5, PS5>&)
  * \see lu_solve (Matrix<T,PO1,PS1>, const Matrix<T,PO2,PS2>&)
  *
  * \throw scythe_alloc_error (Level 1)
  * \throw scythe_null_error (Level 1)
  * \throw scythe_conformation_error (Level 1)
  * \throw scythe_dimension_error (Level 1)
  * \throw scythe_type_error (Level 2)
  * \throw scythe_alloc_error (Level 1)
  *
  */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,RO,RS>
  chol_solve (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& b)
  {
    /* NOTE: cholesky() call does check for square/posdef of A,
     * and the overloaded chol_solve call handles dimensions
     */
  
    return chol_solve<RO,RS>(A, b, cholesky<RO,Concrete>(A));
  }
 
  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,PO1,Concrete>
  chol_solve (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& b)
  {
    return chol_solve<PO1,Concrete>(A, b);
  }

  
  /*!\brief Calculates the inverse of a symmetric positive definite
   * matrix, given a lower triangular matrix resulting from Cholesky
   * decomposition.
   *
   * This function returns the inverse of a symmetric positive
   * definite matrix. Unlike the one-parameter version, this function
   * requires the caller to perform Cholesky decomposition on the
   * matrix to invert, ahead of time.
   *
   * \param A The symmetric positive definite matrix to invert.
   * \param M The lower triangular matrix from the Cholesky decomposition of \a A.
   *
   * \see invpd(const Matrix<T, PO, PS>&)
   * \see inv(const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&, const Matrix<T,PO3,PS3>&, const Matrix<unsigned int,PO4,PS4>&)
   * \see inv(const Matrix<T, PO, PS>&)
   * \see cholesky(const Matrix<T, PO, PS>&)
   *
   * \throw scythe_alloc_error (Level 1)
   * \throw scythe_null_error (Level 1)
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_dimension_error (Level 1)
   */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,RO,RS>
  invpd (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& M)
  {
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "A is NULL")
    SCYTHE_CHECK_10(! A.isSquare(), scythe_dimension_error,
        "A is not square")
    SCYTHE_CHECK_10(A.rows() != M.cols() || A.cols() != M.rows(), 
        scythe_conformation_error, "A and M do not conform");
      
    // for chol_solve block
    T *y = new T[A.rows()];
    T *x = new T[A.rows()];
    Matrix<T, RO, Concrete> b(A.rows(), 1); // full of zeros
    Matrix<T, RO, Concrete> null;
    
    // For final answer
    Matrix<T, RO, Concrete> Ainv(A.rows(), A.cols(), false);

    for (uint k = 0; k < A.rows(); ++k) {
      b[k] = (T) 1;

      solve(M, null, b, x, y);

      b[k] = (T) 0;
      for (uint l = 0; l < A.rows(); ++l)
        Ainv(l,k) = x[l];
    }

    delete[] y;
    delete[] x;

    SCYTHE_VIEW_RETURN(T, RO, RS, Ainv)
  }

  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,PO1,Concrete>
  invpd (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& M)
  {
    return invpd<PO1,Concrete>(A, M);
  }

  /*!\brief Calculate the inverse of a symmetric positive definite
   * matrix.
  *
  * This function returns the inverse of a symmetric positive definite
  * matrix, using cholesky() to do the necessary decomposition. This
  * method is significantly faster than the generalized inverse
  * function.
  *
  * \param A The symmetric positive definite matrix to invert.
  *
  * \see invpd(const Matrix<T, PO1, PS1>&, const Matrix<T, PO2, PS2>&)
  * \see inv (const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&, const Matrix<T,PO3,PS3>&, const Matrix<unsigned int,PO4,PS4>&)
  * \see inv (const Matrix<T, PO, PS>&)
  *
  * \throw scythe_alloc_error (Level 1)
  * \throw scythe_null_error (Level 1)
  * \throw scythe_conformation_error (Level 1)
  * \throw scythe_dimension_error (Level 1)
  * \throw scythe_type_error (Level 2)
  */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  invpd (const Matrix<T, PO, PS>& A)
  { 
    // Cholesky checks to see if A is square and symmetric
  
    return invpd<RO,RS>(A, cholesky<RO,Concrete>(A));
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  invpd (const Matrix<T,O,S>& A)
  {
    return invpd<O,Concrete>(A);
  }

  /* This code is based on  Algorithm 3.4.1 of Golub and Van Loan 3rd
   * edition, 1996. Major difference is in how the output is
   * structured.  Returns the sign of the row permutation (used by
   * det).  Internal function, doesn't need doxygen.
   */
  namespace {
    template <matrix_order PO1, matrix_style PS1, typename T,
              matrix_order PO2, matrix_order PO3, matrix_order PO4>
    inline T
    lu_decomp_alg(Matrix<T,PO1,PS1>& A, Matrix<T,PO2,Concrete>& L, 
                  Matrix<T,PO3,Concrete>& U, 
                  Matrix<unsigned int, PO4, Concrete>& perm_vec)
    {
      if (A.isRowVector()) {
        L = Matrix<T,PO2,Concrete> (1, 1, true, 1); // all 1s
        U = A;
        perm_vec = Matrix<uint, PO4, Concrete>(1, 1);  // all 0s
        return (T) 0;
      }
      
      L = U = Matrix<T, PO2, Concrete>(A.rows(), A.cols(), false);
      perm_vec = Matrix<uint, PO3, Concrete> (A.rows() - 1, 1, false);

      uint pivot;
      T temp;
      T sign = (T) 1;

      for (uint k = 0; k < A.rows() - 1; ++k) {
        pivot = k;
        // find pivot
        for (uint i = k; i < A.rows(); ++i) {
          if (std::fabs(A(pivot,k)) < std::fabs(A(i,k)))
            pivot = i;
        }
        
        SCYTHE_CHECK_20(A(pivot,k) == (T) 0, scythe_type_error,
            "Matrix is singular");

        // permute
        if (k != pivot) {
          sign *= -1;
          for (uint i = 0; i < A.rows(); ++i) {
            temp = A(pivot,i);
            A(pivot,i) = A(k,i);
            A(k,i) = temp;
          }
        }
        perm_vec[k] = pivot;

        for (uint i = k + 1; i < A.rows(); ++i) {
          A(i,k) = A(i,k) / A(k,k);
          for (uint j = k + 1; j < A.rows(); ++j)
            A(i,j) = A(i,j) - A(i,k) * A(k,j);
        }
      }

      L = A;

      for (uint i = 0; i < A.rows(); ++i) {
        for (uint j = i; j < A.rows(); ++j) {
          U(i,j) = A(i,j);
          L(i,j) = (T) 0;
          L(i,i) = (T) 1;
        }
      }
      return sign;
    }
  }

  /* Calculates the LU Decomposition of a square Matrix */

  /* Note that the L, U, and perm_vec must be concrete. A is passed by
   * value, because it is changed during the decomposition.  If A is a
   * view, it will get mangled, but the decomposition will work fine.
   * Not sure what the copy/view access trade-off is, but passing a
   * view might speed things up if you don't care about messing up
   * your matrix.
   */
    /*! \brief LU decomposition of a square matrix.
     *
     * This function performs LU decomposition. That is, given a
     * non-singular square matrix \a A and three matrix references, \a
     * L, \a U, and \a perm_vec, lu_decomp fills the latter three
     * matrices such that \f$LU = A\f$. This method does not actually
     * calculate the LU decomposition of \a A, but of a row-wise
     * permutation of \a A. This permutation is recorded in perm_vec.
     *
     * \note Note that \a L, \a U, and \a perm_vec must be concrete.
     * \a A is passed by value because the function modifies it during
     * the decomposition.  Users should generally avoid passing Matrix
     * views as the first argument to this function because this
     * results in modification to the Matrix being viewed.
     *
     * \param A Non-singular square matrix to decompose.
     * \param L Lower triangular portion of LU decomposition of A.
     * \param U Upper triangular portion of LU decomposition of A.
     * \param perm_vec Permutation vector recording the row-wise permutation of A actually decomposed by the algorithm.
     *
     * \see cholesky (const Matrix<T, PO, PS>&)
     * \see lu_solve (const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&, const Matrix<T,PO3,PS3>&, const Matrix<T,PO4,PS4>&, const Matrix<unsigned int, PO5, PS5>&)
     * \see lu_solve (Matrix<T,PO1,PS1>, const Matrix<T,PO2,PS2>&)
     *
     * \throw scythe_null_error (Level 1)
     * \throw scythe_dimension_error (Level 1)
     * \throw scythe_type_error (Level 2)
     */
  template <matrix_order PO1, matrix_style PS1, typename T,
            matrix_order PO2, matrix_order PO3, matrix_order PO4>
  void
  lu_decomp(Matrix<T,PO1,PS1> A, Matrix<T,PO2,Concrete>& L, 
            Matrix<T,PO3,Concrete>& U, 
            Matrix<unsigned int, PO4, Concrete>& perm_vec)
  {
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "A is NULL")
    SCYTHE_CHECK_10(! A.isSquare(), scythe_dimension_error,
        "Matrix A not square");

    lu_decomp_alg(A, L, U, perm_vec);
  }

  /* lu_solve overloaded: you need A, b + L, U, perm_vec from
   * lu_decomp.
   *
   */
    /*! \brief Solve \f$Ax=b\f$ for x via forward and backward
     * substitution, given the results of a LU decomposition.
     *
     * This function solves the system of equations \f$Ax = b\f$ via
     * forward and backward substitution and LU decomposition. \a A
     * must be a non-singular square matrix for this method to work.
     * This function requires the actual LU decomposition to be
     * performed ahead of time; by lu_decomp() for example.
     *
     * This function is intended for repeatedly solving systems of
     * equations based on \a A.  That is \a A stays constant while \a
     * b varies.
     *
     * \param A Non-singular square Matrix to decompose, passed by reference.
     * \param b Column vector with as many rows as \a A.
     * \param L Lower triangular portion of LU decomposition of \a A.
     * \param U Upper triangular portion of LU decomposition of \a A.
     * \param perm_vec Permutation vector recording the row-wise permutation of \a A actually decomposed by the algorithm.
     *
     * \see lu_solve (Matrix<T,PO1,PS1>, const Matrix<T,PO2,PS2>&)
     * \see lu_decomp(Matrix<T,PO1,PS1>, Matrix<T,PO2,Concrete>&, Matrix<T,PO3,Concrete>&, Matrix<unsigned int, PO4, Concrete>&)
   * \see chol_solve(const Matrix<T,PO1,PS1> &, const Matrix<T,PO2,PS2> &)
   * \see chol_solve(const Matrix<T,PO1,PS1> &, const Matrix<T,PO2,PS2> &, const Matrix<T,PO3,PS3> &)
     *
     * \throw scythe_null_error (Level 1)
     * \throw scythe_dimension_error (Level 1)
     * \throw scythe_conformation_error (Level 1)
     */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2,
            matrix_order PO3, matrix_style PS3,
            matrix_order PO4, matrix_style PS4,
            matrix_order PO5, matrix_style PS5>
  Matrix<T, RO, RS>
  lu_solve (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& b,
            const Matrix<T,PO3,PS3>& L, const Matrix<T,PO4,PS4>& U,
            const Matrix<unsigned int, PO5, PS5> &perm_vec) 
  {
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "A is NULL")
    SCYTHE_CHECK_10(! b.isColVector(), scythe_dimension_error,
        "b is not a column vector");
    SCYTHE_CHECK_10(! A.isSquare(), scythe_dimension_error,
        "A is not square");
    SCYTHE_CHECK_10(A.rows() != b.rows(), scythe_conformation_error,
        "A and b have different row sizes");
    SCYTHE_CHECK_10(A.rows() != L.rows() || A.rows() != U.rows() ||
                    A.cols() != L.cols() || A.cols() != U.cols(),
                    scythe_conformation_error,
                    "A, L, and U do not conform");
    SCYTHE_CHECK_10(perm_vec.rows() + 1 != A.rows(),
        scythe_conformation_error,
        "perm_vec does not have exactly one less row than A");

    T *y = new T[A.rows()];
    T *x = new T[A.rows()];
    
    Matrix<T,RO,Concrete> bb = row_interchange(b, perm_vec);
    solve(L, U, bb, x, y);

    Matrix<T,RO,RS> result(A.rows(), 1, x);
     
    delete[]x;
    delete[]y;
   
    return result;
  }

  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2,
            matrix_order PO3, matrix_style PS3,
            matrix_order PO4, matrix_style PS4,
            matrix_order PO5, matrix_style PS5>
  Matrix<T, PO1, Concrete>
  lu_solve (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& b,
            const Matrix<T,PO3,PS3>& L, const Matrix<T,PO4,PS4>& U,
            const Matrix<unsigned int, PO5, PS5> &perm_vec) 
  {
    return lu_solve<PO1,Concrete>(A, b, L, U, perm_vec);
  }

    /*! \brief Solve \f$Ax=b\f$ for x via forward and backward
     * substitution, using LU decomposition 
     *
     * This function solves the system of equations \f$Ax = b\f$ via
     * forward and backward substitution and LU decomposition. \a A
     * must be a non-singular square matrix for this method to work.
     *
     * \param A A non-singular square Matrix to decompose.
     * \param b A column vector with as many rows as \a A.
     *
     * \see lu_solve (const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&, const Matrix<T,PO3,PS3>&, const Matrix<T,PO4,PS4>&, const Matrix<unsigned int, PO5, PS5>&) 
     * \see lu_decomp(Matrix<T,PO1,PS1>, Matrix<T,PO2,Concrete>&, Matrix<T,PO3,Concrete>&, Matrix<unsigned int, PO4, Concrete>&)
   * \see chol_solve(const Matrix<T,PO1,PS1> &, const Matrix<T,PO2,PS2> &)
   * \see chol_solve(const Matrix<T,PO1,PS1> &, const Matrix<T,PO2,PS2> &, const Matrix<T,PO3,PS3> &)
   *
     * \throw scythe_null_error (Level 1)
     * \throw scythe_dimension_error (Level 1)
     * \throw scythe_conformation_error (Level 1)
     * \throw scythe_type_error (Level 2)
     */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,RO,RS>
  lu_solve (Matrix<T,PO1,PS1> A, const Matrix<T,PO2,PS2>& b)
  {
    // step 1 compute the LU factorization 
    Matrix<T, RO, Concrete> L, U;
    Matrix<uint, RO, Concrete> perm_vec;
    lu_decomp_alg(A, L, U, perm_vec);

    return lu_solve<RO,RS>(A, b, L, U, perm_vec);
  }

  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,PO1,Concrete>
  lu_solve (Matrix<T,PO1,PS1> A, const Matrix<T,PO2,PS2>& b)
  {
    // Slight code rep here, but very few lines
   
    // step 1 compute the LU factorization 
    Matrix<T, PO1, Concrete> L, U;
    Matrix<uint, PO1, Concrete> perm_vec;
    lu_decomp_alg(A, L, U, perm_vec);
    
    return lu_solve<PO1,Concrete>(A, b, L, U, perm_vec);
  }

 /*!\brief Calculates the inverse of a non-singular square matrix,
  * given an LU decomposition.
  *
  * This function returns the inverse of an arbitrary, non-singular,
  * square matrix \a A when passed a permutation of an LU
  * decomposition, such as that returned by lu_decomp().  A
  * one-parameter version of this function exists that does not
  * require the user to pre-decompose the system.
  *
  * \param A The Matrix to be inverted.
  * \param L A Lower triangular matrix resulting from decomposition.
  * \param U An Upper triangular matrix resulting from decomposition.
  * \param perm_vec The permutation vector recording the row-wise permutation of \a A actually decomposed by the algorithm.
  *
  * \see inv (const Matrix<T, PO, PS>&)
  * \see invpd(const Matrix<T, PO, PS>&)
  * \see invpd(const Matrix<T, PO1, PS1>&, const Matrix<T, PO2, PS2>&)
  * \see lu_decomp(Matrix<T,PO1,PS1>, Matrix<T,PO2,Concrete>&, Matrix<T,PO3,Concrete>&, Matrix<unsigned int, PO4, Concrete>&)
  *
  * \throw scythe_null_error(Level 1)
  * \throw scythe_dimension_error (Level 1)
  * \throw scythe_conformation_error (Level 1)
  */
  template<matrix_order RO, matrix_style RS, typename T,
           matrix_order PO1, matrix_style PS1,
           matrix_order PO2, matrix_style PS2,
           matrix_order PO3, matrix_style PS3,
           matrix_order PO4, matrix_style PS4>
  Matrix<T,RO,RS>
  inv (const Matrix<T,PO1,PS1>& A, 
       const Matrix<T,PO2,PS2>& L, const Matrix<T,PO3,PS3>& U,
       const Matrix<unsigned int,PO4,PS4>& perm_vec)
  {
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "A is NULL")
    SCYTHE_CHECK_10 (! A.isSquare(), scythe_dimension_error,
        "A is not square");
    SCYTHE_CHECK_10(A.rows() != L.rows() || A.rows() != U.rows() ||
                    A.cols() != L.cols() || A.cols() != U.cols(),
                    scythe_conformation_error,
                    "A, L, and U do not conform");
    SCYTHE_CHECK_10(perm_vec.rows() + 1 != A.rows() 
        && !(A.isScalar() && perm_vec.isScalar()),
        scythe_conformation_error,
        "perm_vec does not have exactly one less row than A");

    // For the final result
    Matrix<T,RO,Concrete> Ainv(A.rows(), A.rows(), false);

    // for the solve block
    T *y = new T[A.rows()];
    T *x = new T[A.rows()];
    Matrix<T, RO, Concrete> b(A.rows(), 1); // full of zeros
    Matrix<T,RO,Concrete> bb;
    
    for (uint k = 0; k < A.rows(); ++k) {
      b[k] = (T) 1;
      bb = row_interchange(b, perm_vec);

      solve(L, U, bb, x, y);

      b[k] = (T) 0;
      for (uint l = 0; l < A.rows(); ++l)
        Ainv(l,k) = x[l];
    }

    delete[] y;
    delete[] x;

    SCYTHE_VIEW_RETURN(T, RO, RS, Ainv)
  }

  template<typename T,
           matrix_order PO1, matrix_style PS1,
           matrix_order PO2, matrix_style PS2,
           matrix_order PO3, matrix_style PS3,
           matrix_order PO4, matrix_style PS4>
  Matrix<T,PO1,Concrete>
  inv (const Matrix<T,PO1,PS1>& A, 
       const Matrix<T,PO2,PS2>& L, const Matrix<T,PO3,PS3>& U,
       const Matrix<unsigned int,PO4,PS4>& perm_vec)
  {
    return inv<PO1,Concrete>(A, L, U, perm_vec);
  }

 /*!\brief Invert an arbitrary, non-singular, square matrix.
  *
  * This function returns the inverse of a non-singular square matrix,
  * using lu_decomp() to do the necessary decomposition.  This method
  * is significantly slower than the inverse function for symmetric
  * positive definite matrices, invpd().
  *
  * \param A The Matrix to be inverted.
  *
  * \see inv (const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&, const Matrix<T,PO3,PS3>&, const Matrix<unsigned int,PO4,PS4>&)
  * \see invpd(const Matrix<T, PO, PS>&)
  * \see invpd(const Matrix<T, PO1, PS1>&, const Matrix<T, PO2, PS2>&)
  *
  * \throw scythe_null_error(Level 1)
  * \throw scythe_dimension_error (Level 1)
  * \throw scythe_conformation_error (Level 1)
  * \throw scythe_type_error (Level 2)
  */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  inv (const Matrix<T, PO, PS>& A)
  {
    // Make a copy of A for the decomposition (do it with an explicit
    // copy to a concrete case A is a view)
    Matrix<T,RO,Concrete> AA = A;
    
    // step 1 compute the LU factorization 
    Matrix<T, RO, Concrete> L, U;
    Matrix<uint, RO, Concrete> perm_vec;
    lu_decomp_alg(AA, L, U, perm_vec);

    return inv<RO,RS>(A, L, U, perm_vec);
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  inv (const Matrix<T, O, S>& A)
  {
    return inv<O,Concrete>(A);
  }

  /* Interchanges the rows of A with those in vector p */
  /*!\brief Interchange the rows of a Matrix according to a
  * permutation vector.
  *
  * This function permutes the rows of Matrix \a A according to \a
  * perm_vec.  Each element i of perm_vec contains a row-number, r.
  * For each row, i, in \a A, A[i] is interchanged with A[r].
  *
  * \param A The matrix to permute.
  * \param p The column vector describing the permutations to perform
  * on \a A.
  *
  * \see lu_decomp(Matrix<T,PO1,PS1>, Matrix<T,PO2,Concrete>&, Matrix<T,PO3,Concrete>&, Matrix<unsigned int, PO4, Concrete>&)
  *
  * \throw scythe_dimension_error (Level 1)
  * \throw scythe_conformation_error (Level 1)
  */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,RO,RS>
  row_interchange (Matrix<T,PO1,PS1> A, 
                   const Matrix<unsigned int,PO2,PS2>& p)
  {
    SCYTHE_CHECK_10(! p.isColVector(), scythe_dimension_error,
        "p not a column vector");
    SCYTHE_CHECK_10(p.rows() + 1 != A.rows() && ! p.isScalar(), 
        scythe_conformation_error, "p must have one less row than A");

    for (uint i = 0; i < A.rows() - 1; ++i) {
      Matrix<T,PO1,View> vec1 = A(i, _);
      Matrix<T,PO1,View> vec2 = A(p[i], _);
      std::swap_ranges(vec1.begin_f(), vec1.end_f(), vec2.begin_f());
    }
    
    return A;
  }
  
  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,PO1,Concrete>
  row_interchange (const Matrix<T,PO1,PS1>& A, 
                   const Matrix<unsigned int,PO2,PS2>& p)
  {
    return row_interchange<PO1,Concrete>(A, p);
  }

  /*! \brief Calculate the determinant of a square Matrix.
   *
   * This routine calculates the determinant of a square Matrix, using
   * LU decomposition.
   *
   * \param A The Matrix to calculate the determinant of.
   *
  * \see lu_decomp(Matrix<T,PO1,PS1>, Matrix<T,PO2,Concrete>&, Matrix<T,PO3,Concrete>&, Matrix<unsigned int, PO4, Concrete>&)
  *
   * \throws scythe_dimension_error (Level 1)
   * \throws scythe_null_error (Level 1)
   */
  template <typename T, matrix_order PO, matrix_style PS>
  T
  det (const Matrix<T, PO, PS>& A)
  {
    SCYTHE_CHECK_10(! A.isSquare(), scythe_dimension_error,
        "Matrix is not square")
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "Matrix is NULL")
    
    // Make a copy of A for the decomposition (do it here instead of
    // at parameter pass in case A is a view)
    Matrix<T,PO,Concrete> AA = A;
    
    // step 1 compute the LU factorization 
    Matrix<T, PO, Concrete> L, U;
    Matrix<uint, PO, Concrete> perm_vec;
    T sign = lu_decomp_alg(AA, L, U, perm_vec);

    // step 2 calculate the product of diag(U) and sign
    T det = (T) 1;
    for (uint i = 0; i < AA.rows(); ++i)
      det *= AA(i, i);

    return sign * det;
  }

#ifdef SCYTHE_LAPACK

  template<>
  Matrix<>
  cholesky (const Matrix<>& A)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for cholesky");
    SCYTHE_CHECK_10(! A.isSquare(), scythe_dimension_error,
        "Matrix not square");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "Matrix is NULL");

    // We have to do an explicit copy within the func to match the
    // template declaration of the more general template.
    Matrix<> AA = A;

    // Get a pointer to the internal array and set up some vars
    double* Aarray = AA.getArray();  // internal array pointer
    int rows = (int) AA.rows(); // the dim of the matrix
    int err = 0; // The output error condition

    // Cholesky decomposition step
    lapack::dpotrf_("L", &rows, Aarray, &rows, &err);
    SCYTHE_CHECK_10(err > 0, scythe_type_error,
        "Matrix is not positive definite")
    SCYTHE_CHECK_10(err < 0, scythe_invalid_arg,
        "The " << err << "th value of the matrix had an illegal value")

    // Zero out upper triangle
    for (uint j = 1; j < AA.cols(); ++j)
      for (uint i = 0; i < j; ++i)
        AA(i, j) = 0;

    return AA;
  }

  template<>
  Matrix<>
  chol_solve (const Matrix<>& A, const Matrix<>& b, const Matrix<>& M)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for chol_solve");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "A is NULL")
    SCYTHE_CHECK_10(! b.isColVector(), scythe_dimension_error,
        "b must be a column vector");
    SCYTHE_CHECK_10(A.rows() != b.rows(), scythe_conformation_error,
        "A and b do not conform");
    SCYTHE_CHECK_10(A.rows() != M.rows(), scythe_conformation_error,
        "A and M do not conform");
    SCYTHE_CHECK_10(! M.isSquare(), scythe_dimension_error,
        "M must be square");

    // The algorithm modifies b in place.  We make a copy.
    Matrix<> bb = b;

    // Get array pointers and set up some vars
    const double* Marray = M.getArray();
    double* barray = bb.getArray();
    int rows = (int) bb.rows();
    int cols = (int) bb.cols(); // currently always one, but generalizable
    int err = 0;

    // Solve the system
    lapack::dpotrs_("L", &rows, &cols, Marray, &rows, barray, &rows, &err);
    SCYTHE_CHECK_10(err > 0, scythe_type_error,
        "Matrix is not positive definite")
    SCYTHE_CHECK_10(err < 0, scythe_invalid_arg,
        "The " << err << "th value of the matrix had an illegal value")

    return bb;
  }

  template<>
  Matrix<>
  chol_solve (const Matrix<>& A, const Matrix<>& b)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for chol_solve");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "A is NULL")
    SCYTHE_CHECK_10(! b.isColVector(), scythe_dimension_error,
        "b must be a column vector");
    SCYTHE_CHECK_10(A.rows() != b.rows(), scythe_conformation_error,
        "A and b do not conform");

    // The algorithm modifies both A and b in place, so we make copies
    Matrix<> AA =A;
    Matrix<> bb = b;

    // Get array pointers and set up some vars
    double* Aarray = AA.getArray();
    double* barray = bb.getArray();
    int rows = (int) bb.rows();
    int cols = (int) bb.cols(); // currently always one, but generalizable
    int err = 0;

    // Solve the system
    lapack::dposv_("L", &rows, &cols, Aarray, &rows, barray, &rows, &err);
    SCYTHE_CHECK_10(err > 0, scythe_type_error,
        "Matrix is not positive definite")
    SCYTHE_CHECK_10(err < 0, scythe_invalid_arg,
        "The " << err << "th value of the matrix had an illegal value")

    return bb;
  }

  template <matrix_order PO2, matrix_order PO3, matrix_order PO4>
  inline double
  lu_decomp_alg(Matrix<>& A, Matrix<double,PO2,Concrete>& L, 
                Matrix<double,PO3,Concrete>& U, 
                Matrix<unsigned int, PO4, Concrete>& perm_vec)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for lu_decomp_alg");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error, "A is NULL")
    SCYTHE_CHECK_10 (! A.isSquare(), scythe_dimension_error,
        "A is not square");

    if (A.isRowVector()) {
      L = Matrix<double,PO2,Concrete> (1, 1, true, 1); // all 1s
      U = A;
      perm_vec = Matrix<uint, PO4, Concrete>(1, 1);  // all 0s
      return 0.;
    }
      
    L = U = Matrix<double, PO2, Concrete>(A.rows(), A.cols(), false);
    perm_vec = Matrix<uint, PO3, Concrete> (A.rows(), 1, false);

    // Get a pointer to the internal array and set up some vars
    double* Aarray = A.getArray();  // internal array pointer
    int rows = (int) A.rows(); // the dim of the matrix
    int* ipiv = (int*) perm_vec.getArray(); // Holds the lu decomp pivot array
    int err = 0; // The output error condition

    // Do the decomposition
    lapack::dgetrf_(&rows, &rows, Aarray, &rows, ipiv, &err);

    SCYTHE_CHECK_10(err > 0, scythe_type_error, "Matrix is singular");
    SCYTHE_CHECK_10(err < 0, scythe_lapack_internal_error, 
        "The " << err << "th value of the matrix had an illegal value");

    // Now fill in the L and U matrices.
    L = A;
    for (uint i = 0; i < A.rows(); ++i) {
      for (uint j = i; j < A.rows(); ++j) {
        U(i,j) = A(i,j);
        L(i,j) = 0.;
        L(i,i) = 1.;
      }
    }

    // Change to scythe's rows-1 perm_vec format and c++ indexing
    // XXX Cutting off the last pivot term may be buggy if it isn't
    // always just pointing at itself
    if (perm_vec(perm_vec.size() - 1) != perm_vec.size())
      SCYTHE_THROW(scythe_unexpected_default_error,
          "This is an unexpected error.  Please notify the developers.")
    perm_vec = perm_vec(0, 0, perm_vec.rows() - 2, 0) - 1;

    // Finally, figure out the sign of perm_vec
    if (sum(perm_vec > 0) % 2 == 0)
      return 1;

    return -1;
  }

  /*! \brief The result of a QR decomposition.
   *
   * Objects of this type contain three matrices, \a QR, \a tau, and
   * \a pivot, representing the results of a QR decomposition of a
   * \f$m \times n\f$ matrix.  After decomposition, the upper triangle
   * of \a QR contains the min(\f$m\f$, \f$n\f$) by \f$n\f$ upper
   * trapezoidal matrix \f$R\f$, while \a tau and the elements of \a
   * QR below the diagonal represent the orthogonal matrix \f$Q\f$ as
   * a product of min(\f$m\f$, \f$n\f$) elementary reflectors.  The
   * vector \a pivot is a permutation vector containing information
   * about the pivoting strategy used in the factorization.
   *
   * \a QR is \f$m \times n\f$, tau is a vector of dimension
   * min(\f$m\f$, \f$n\f$), and pivot is a vector of dimension
   * \f$n\f$.
   *
   * \see qr_decomp (const Matrix<>& A)
   */

  struct QRdecomp {
    Matrix<> QR;
    Matrix<> tau;
    Matrix<> pivot;
  };

  /*! \brief QR decomposition of a matrix.
   *
   * This function performs QR decomposition.  That is, given a 
   * \f$m \times n \f$ matrix \a A, qr_decomp computes the QR factorization
   * of \a A with column pivoting, such that \f$A \cdot P = Q \cdot
   * R\f$.  The resulting QRdecomp object contains three matrices, \a
   * QR, \a tau, and \a pivot.  The upper triangle of \a QR contains the
   * min(\f$m\f$, \f$n\f$) by \f$n\f$ upper trapezoidal matrix
   * \f$R\f$, while \a tau and the elements of \a QR below the
   * diagonal represent the orthogonal matrix \f$Q\f$ as a product of
   * min(\f$m\f$, \f$n\f$) elementary reflectors.  The vector \a pivot
   * is a permutation vector containing information about the pivoting
   * strategy used in the factorization.
   *
   * \note This function requires BLAS/LAPACK functionality and is
   * only available on machines that provide these libraries.  Make
   * sure you enable the SCYTHE_LAPACK preprocessor flag if you wish
   * to use this function.  Furthermore, note that this function takes
   * and returns only column-major concrete matrices.  Future versions
   * of Scythe will provide a native C++ implementation of this
   * function with support for general matrix templates.
   *
   * \param A A matrix to decompose.
   *
   * \see QRdecomp
   * \see lu_decomp(Matrix<T,PO1,PS1>, Matrix<T,PO2,Concrete>&, Matrix<T,PO3,Concrete>&, Matrix<unsigned int, PO4, Concrete>&)
   * \see cholesky (const Matrix<T, PO, PS>&)
   * \see qr_solve (const Matrix<>& A, const Matrix<>& b, const QRdecomp& QR)
   * \see qr_solve (const Matrix<>& A, const Matrix<>& b);
   *
   * \throw scythe_null_error (Level 1)
   * \throw scythe_lapack_internal_error (Level 1)
   */
  QRdecomp
  qr_decomp (const Matrix<>& A)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for qr_decomp");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error, "A is NULL");

    // Set up working variables
    Matrix<> QR = A;
    double* QRarray = QR.getArray(); // input/output array pointer
    int rows = (int) QR.rows();
    int cols = (int) QR.cols();
    Matrix<unsigned int> pivot(cols, 1); // pivot vector
    int* parray = (int*) pivot.getArray(); // pivot vector array pointer
    Matrix<> tau = Matrix<>(rows < cols ? rows : cols, 1);
    double* tarray = tau.getArray(); // tau output array pointer
    double tmp, *work; // workspace vars
    int lwork, info;   // workspace size var and error info var

    // Get workspace size
    lwork = -1;
    lapack::dgeqp3_(&rows, &cols, QRarray, &rows, parray, tarray, &tmp,
                    &lwork, &info);

    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dgeqp3");

    lwork = (int) tmp;
    work = new double[lwork];

    // run the routine for real
    lapack::dgeqp3_(&rows, &cols, QRarray, &rows, parray, tarray, work,
                    &lwork, &info);

    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dgeqp3");

    delete[] work;

    pivot -= 1;

    QRdecomp result;
    result.QR = QR;
    result.tau = tau;
    result.pivot = pivot;

    return result;
  }

  /*! \brief Solve \f$Ax=b\f$ given a QR decomposition.
   *
   * This function solves the system of equations \f$Ax = b\f$ using
   * the results of a QR decomposition.  This function requires the
   * actual QR decomposition to be performed ahead of time; by
   * qr_decomp() for example.
   *
   * This function is intended for repeatedly solving systems of
   * equations based on \a A.  That is \a A stays constant while \a b
   * varies.
   *
   * \note This function requires BLAS/LAPACK functionality and is
   * only available on machines that provide these libraries.  Make
   * sure you enable the SCYTHE_LAPACK preprocessor flag if you wish
   * to use this function.  Furthermore, note that this function takes
   * and returns only column-major concrete matrices.  Future versions
   * of Scythe will provide a native C++ implementation of this
   * function with support for general matrix templates.
   *
   * \param A A  Matrix to decompose.
   * \param b A Matrix with as many rows as \a A.
   * \param QR A QRdecomp object containing the result of the QR decomposition of \a A.
   *
   * \see QRdecomp
   * \see qr_solve (const Matrix<>& A, const Matrix<>& b)
   * \see qr_decomp (const Matrix<>& A)
   * \see lu_solve (const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&, const Matrix<T,PO3,PS3>&, const Matrix<T,PO4,PS4>&, const Matrix<unsigned int, PO5, PS5>&)
   * \see lu_solve (Matrix<T,PO1,PS1>, const Matrix<T,PO2,PS2>&)
   * \see chol_solve(const Matrix<T,PO1,PS1> &, const Matrix<T,PO2,PS2> &)
   * \see chol_solve(const Matrix<T,PO1,PS1> &, const Matrix<T,PO2,PS2> &, const Matrix<T,PO3,PS3> &)
   *
   * \throw scythe_null_error (Level 1) 
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_type_error (Level 1)
   * \throw scythe_lapack_internal_error (Level 1)
   */
  inline Matrix<> 
  qr_solve(const Matrix<>& A, const Matrix<>& b, const QRdecomp& QR)
  { 
    SCYTHE_DEBUG_MSG("Using lapack/blas for qr_solve");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error, "A is NULL")
    SCYTHE_CHECK_10(A.rows() != b.rows(), scythe_conformation_error,
          "A and b do not conform");
    SCYTHE_CHECK_10(A.rows() != QR.QR.rows() || A.cols() != QR.QR.cols(),
        scythe_conformation_error, "A and QR do not conform"); 
    int taudim = (int) (A.rows() < A.cols() ? A.rows() : A.cols());
    SCYTHE_CHECK_10(QR.tau.size() != taudim, scythe_conformation_error,
        "A and tau do not conform");
    SCYTHE_CHECK_10(QR.pivot.size() != A.cols(), scythe_conformation_error,
        "pivot vector is not the right length");

    int rows = (int) QR.QR.rows();
    int cols = (int) QR.QR.cols();
    int nrhs = (int) b.cols();
    int lwork, info;
    double *work, tmp;
    double* QRarray = QR.QR.getArray();
    double* tarray = QR.tau.getArray();
    Matrix<> bb = b;
    double* barray = bb.getArray();

    // Get workspace size
    lwork = -1;
    lapack::dormqr_("L", "T", &rows, &nrhs, &taudim, QRarray, &rows,
                    tarray, barray, &rows, &tmp, &lwork, &info);

    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dormqr");

    // And now for real
    lwork = (int) tmp;
    work = new double[lwork];
    lapack::dormqr_("L", "T", &rows, &nrhs, &taudim, QRarray, &rows,
                    tarray, barray, &rows, work, &lwork, &info);

    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dormqr");

    lapack::dtrtrs_("U", "N", "N", &taudim, &nrhs, QRarray, &rows, barray,
                    &rows, &info);

    SCYTHE_CHECK_10(info > 0, scythe_type_error, "Matrix is singular");
    SCYTHE_CHECK_10(info < 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dtrtrs");

    delete[] work;

    Matrix<> result(A.cols(), b.cols(), false);
    for (uint i = 0; i < QR.pivot.size(); ++i)
      result(i, _) = bb(QR.pivot(i), _);
    return result;
  }

  /*! \brief Solve \f$Ax=b\f$ using QR decomposition.
   *
   * This function solves the system of equations \f$Ax = b\f$ using
   * QR decomposition.  This function is intended for repeatedly
   * solving systems of equations based on \a A.  That is \a A stays
   * constant while \a b varies.
   *
   * \note This function used BLAS/LAPACK support functionality and is
   * only available on machines that provide these libraries.  Make
   * sure you enable the SCYTHE_LAPACK preprocessor flag if you wish
   * to use this function.  Furthermore, note that the function takes
   * and returns only column-major concrete matrices.  Future versions
   * of Scythe will provide a native C++ implementation of this
   * function with support for general matrix templates.
   *
   * \param A A  Matrix to decompose.
   * \param b A Matrix with as many rows as \a A.
   *
   * \see QRdecomp
   * \see qr_solve (const Matrix<>& A, const Matrix<>& b, const QRdecomp& QR)
   * \see qr_decomp (const Matrix<>& A)
   * \see lu_solve (const Matrix<T,PO1,PS1>&, const Matrix<T,PO2,PS2>&, const Matrix<T,PO3,PS3>&, const Matrix<T,PO4,PS4>&, const Matrix<unsigned int, PO5, PS5>&)
   * \see lu_solve (Matrix<T,PO1,PS1>, const Matrix<T,PO2,PS2>&)
   * \see chol_solve(const Matrix<T,PO1,PS1> &, const Matrix<T,PO2,PS2> &)
   * \see chol_solve(const Matrix<T,PO1,PS1> &, const Matrix<T,PO2,PS2> &, const Matrix<T,PO3,PS3> &)
   *
   * \throw scythe_null_error (Level 1) 
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_type_error (Level 1)
   * \throw scythe_lapack_internal_error (Level 1)
   */
  inline Matrix<>
  qr_solve (const Matrix<>& A, const Matrix<>& b)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for qr_solve");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error, "A is NULL")
    SCYTHE_CHECK_10(A.rows() != b.rows(), scythe_conformation_error,
        "A and b do not conform");

    /* Do decomposition */
   
    // Set up working variables
    Matrix<> QR = A;
    double* QRarray = QR.getArray(); // input/output array pointer
    int rows = (int) QR.rows();
    int cols = (int) QR.cols();
    Matrix<unsigned int> pivot(cols, 1); // pivot vector
    int* parray = (int*) pivot.getArray(); // pivot vector array pointer
    Matrix<> tau = Matrix<>(rows < cols ? rows : cols, 1);
    double* tarray = tau.getArray(); // tau output array pointer
    double tmp, *work; // workspace vars
    int lwork, info;   // workspace size var and error info var

    // Get workspace size
    lwork = -1;
    lapack::dgeqp3_(&rows, &cols, QRarray, &rows, parray, tarray, &tmp,
                    &lwork, &info);

    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dgeqp3");

    lwork = (int) tmp;
    work = new double[lwork];

    // run the routine for real
    lapack::dgeqp3_(&rows, &cols, QRarray, &rows, parray, tarray, work,
                    &lwork, &info);

    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dgeqp3");

    delete[] work;

    pivot -= 1;

    /* Now solve the system */
    
    // working vars
    int nrhs = (int) b.cols();
    Matrix<> bb = b;
    double* barray = bb.getArray();
    int taudim = (int) tau.size();

    // Get workspace size
    lwork = -1;
    lapack::dormqr_("L", "T", &rows, &nrhs, &taudim, QRarray, &rows,
                    tarray, barray, &rows, &tmp, &lwork, &info);

    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dormqr");

    // And now for real
    lwork = (int) tmp;
    work = new double[lwork];
    lapack::dormqr_("L", "T", &rows, &nrhs, &taudim, QRarray, &rows,
                    tarray, barray, &rows, work, &lwork, &info);

    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dormqr");

    lapack::dtrtrs_("U", "N", "N", &taudim, &nrhs, QRarray, &rows, barray,
                    &rows, &info);

    SCYTHE_CHECK_10(info > 0, scythe_type_error, "Matrix is singular");
    SCYTHE_CHECK_10(info < 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dtrtrs");

    delete[] work;

    Matrix<> result(A.cols(), b.cols(), false);
    for (uint i = 0; i < pivot.size(); ++i)
      result(i, _) = bb(pivot(i), _);

    return result;
  }

  template<>
  Matrix<>
  invpd (const Matrix<>& A)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for invpd");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "A is NULL")
    SCYTHE_CHECK_10 (! A.isSquare(), scythe_dimension_error,
        "A is not square");

    // We have to do an explicit copy within the func to match the
    // template declaration of the more general template.
    Matrix<> AA = A;

    // Get a pointer to the internal array and set up some vars
    double* Aarray = AA.getArray();  // internal array pointer
    int rows = (int) AA.rows(); // the dim of the matrix
    int err = 0; // The output error condition

    // Cholesky decomposition step
    lapack::dpotrf_("L", &rows, Aarray, &rows, &err);
    SCYTHE_CHECK_10(err > 0, scythe_type_error,
        "Matrix is not positive definite")
    SCYTHE_CHECK_10(err < 0, scythe_invalid_arg,
        "The " << err << "th value of the matrix had an illegal value")

    // Inversion step
    lapack::dpotri_("L", &rows, Aarray, &rows, &err);
    SCYTHE_CHECK_10(err > 0, scythe_type_error,
        "The (" << err << ", " << err << ") element of the matrix is zero"
        << " and the inverse could not be computed")
    SCYTHE_CHECK_10(err < 0, scythe_invalid_arg,
        "The " << err << "th value of the matrix had an illegal value")
    lapack::make_symmetric(Aarray, rows);

    return AA;
  }

  template<>
  Matrix<>
  invpd (const Matrix<>& A, const Matrix<>& M)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for invpd");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "A is NULL")
    SCYTHE_CHECK_10 (! A.isSquare(), scythe_dimension_error,
        "A is not square");
    SCYTHE_CHECK_10(A.rows() != M.cols() || A.cols() != M.rows(), 
        scythe_conformation_error, "A and M do not conform");

    // We have to do an explicit copy within the func to match the
    // template declaration of the more general template.
    Matrix<> MM = M;

    // Get pointer and set up some vars
    double* Marray = MM.getArray();
    int rows = (int) MM.rows();
    int err = 0;

    // Inversion step
    lapack::dpotri_("L", &rows, Marray, &rows, &err);
    SCYTHE_CHECK_10(err > 0, scythe_type_error,
        "The (" << err << ", " << err << ") element of the matrix is zero"
        << " and the inverse could not be computed")
    SCYTHE_CHECK_10(err < 0, scythe_invalid_arg,
        "The " << err << "th value of the matrix had an illegal value")
    lapack::make_symmetric(Marray, rows);

    return MM;
  }

  template <>
  Matrix<>
  inv(const Matrix<>& A)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for inv");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "A is NULL")
    SCYTHE_CHECK_10 (! A.isSquare(), scythe_dimension_error,
        "A is not square");

    // We have to do an explicit copy within the func to match the
    // template declaration of the more general template.
    Matrix<> AA = A;

    // Get a pointer to the internal array and set up some vars
    double* Aarray = AA.getArray();  // internal array pointer
    int rows = (int) AA.rows(); // the dim of the matrix
    int* ipiv = new int[rows];  // Holds the lu decomp pivot array
    int err = 0; // The output error condition

    // LU decomposition step
    lapack::dgetrf_(&rows, &rows, Aarray, &rows, ipiv, &err);

    SCYTHE_CHECK_10(err > 0, scythe_type_error, "Matrix is singular");
    SCYTHE_CHECK_10(err < 0, scythe_invalid_arg, 
        "The " << err << "th value of the matrix had an illegal value");

    // Inversion step; first do a workspace query, then the actual
    // inversion
    double work_query = 0;
    int work_size = -1;
    lapack::dgetri_(&rows, Aarray, &rows, ipiv, &work_query, 
                    &work_size, &err);
    double* workspace = new double[(work_size = (int) work_query)];
    lapack::dgetri_(&rows, Aarray, &rows, ipiv, workspace, &work_size,
                    &err);
    delete[] ipiv;
    delete[] workspace;

    SCYTHE_CHECK_10(err > 0, scythe_type_error, "Matrix is singular");
    SCYTHE_CHECK_10(err < 0, scythe_invalid_arg, 
        "Internal error in LAPACK routine dgetri");

    return AA;
  }

  /*!\brief The result of a singular value decomposition.
   *
   * Objects of this type hold the results of a singular value
   * decomposition (SVD) of an \f$m \times n\f$ matrix \f$A\f$, as
   * returned by svd().  The SVD takes the form: \f$A = U
   * \cdot \Sigma \cdot V'\f$.  SVD objects contain \a d, which
   * holds the singular values of \f$A\f$ (the diagonal of
   * \f$\Sigma\f$) in descending order.  Furthermore, depending on the
   * options passed to svd(), they may hold some or all of the
   * left singular vectors of \f$A\f$ in \a U and some or all of the
   * right singular vectors of \f$A\f$ in \a Vt.
   *
   * \see svd(const Matrix<>& A, int nu, int nv);
   */
   
  struct SVD {
    Matrix<> d;  // singular values
    Matrix<> U;  // left singular vectors
    Matrix<> Vt; // transpose of right singular vectors
  };

 /*!\brief Calculates the singular value decomposition of a matrix,
  * optionally computing the left and right singular vectors.
  *
  * This function returns the singular value decomposition (SVD) of a
  * \f$m \times n\f$ matrix \a A, optionally computing the left and right
  * singular vectors.  It returns the singular values and vectors in
  * a SVD object.
  *
  * \note This function requires BLAS/LAPACK functionality and is
  * only available on machines that provide these libraries.  Make
  * sure you enable the SCYTHE_LAPACK preprocessor flag if you wish
  * to use this function.  Furthermore, note that this function takes
  * and returns only column-major concrete matrices.  Future versions
  * of Scythe will provide a native C++ implementation of this
  * function with support for general matrix templates.
  *
  * \param A The matrix to decompose.
  * \param nu The number of left singular vectors to compute and return.  Values less than zero are equivalent to min(\f$m\f$, \f$n\f$).
  * \param nv The number of right singular vectors to compute and return.  Values less than zero are equivalent to min(\f$m\f$, \f$n\f$).
  *
  * \throw scythe_null_error (Level 1)
  * \throw scythe_convergence_error (Level 1)
  * \throw scythe_lapack_internal_error (Level 1)
  *
  * \see SVD
  * \see eigen(const Matrix<>& A, bool vectors)
  */
 
 inline SVD
 svd (const Matrix<>& A, int nu = -1, int nv = -1)
 {
   SCYTHE_DEBUG_MSG("Using lapack/blas for eigen");
   SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
       "Matrix is NULL");

   char* jobz;
   int m = (int) A.rows();
   int n = (int) A.cols();
   int mn = (int) std::min(A.rows(), A.cols());
   Matrix<> U;
   Matrix<> V;
   if (nu < 0) nu = mn;
   if (nv < 0) nv = mn;
   if (nu <= mn && nv<= mn) {
     jobz = "S";
     U = Matrix<>(m, mn, false);
     V = Matrix<>(mn, n, false);
   } else if (nu == 0 && nv == 0) {
     jobz = "N";
   } else {
     jobz = "A";
     U = Matrix<>(m, m, false);
     V = Matrix<>(n, n, false);
   }
   double* Uarray = U.getArray();
   double* Varray = V.getArray();

   int ldu = (int) U.rows();
   int ldvt = (int) V.rows();
   Matrix<> X = A;
   double* Xarray = X.getArray();
   Matrix<> d(mn, 1, false);
   double* darray = d.getArray();

   double tmp, *work;
   int lwork, info;
   int *iwork = new int[8 * mn];

   // get optimal workspace
   lwork = -1;
   lapack::dgesdd_(jobz, &m, &n, Xarray, &m, darray, Uarray, &ldu,
                   Varray, &ldvt, &tmp, &lwork, iwork, &info);
   SCYTHE_CHECK_10(info < 0, scythe_lapack_internal_error,
       "Internal error in LAPACK routine dgessd");
   SCYTHE_CHECK_10(info > 0, scythe_convergence_error, "Did not converge");

   lwork = (int) tmp;
   work = new double[lwork];

   // Now for real
   lapack::dgesdd_(jobz, &m, &n, Xarray, &m, darray, Uarray, &ldu,
                   Varray, &ldvt, work, &lwork, iwork, &info);
   SCYTHE_CHECK_10(info < 0, scythe_lapack_internal_error,
       "Internal error in LAPACK routine dgessd");
   SCYTHE_CHECK_10(info > 0, scythe_convergence_error, "Did not converge");
   delete[] work;
  
   if (nu < mn && nu > 0)
     U = U(0, 0, U.rows() - 1, (unsigned int) std::min(m, nu) - 1);
   if (nv < mn && nv > 0)
     V = V(0, 0, (unsigned int) std::min(n, nv) - 1, V.cols() - 1);
   SVD result;
   result.d = d;
   result.U = U;
   result.Vt = V;

   return result;
 }

 /*!\brief The result of an eigenvalue/vector decomposition.
  *
  * Objects of this type hold the results of the eigen() function.
  * That is the eigenvalues and, optionally, the eigenvectors of a
  * symmetric matrix of order \f$n\f$.  The eigenvalues are stored in
  * ascending order in the member column vector \a values.  The
  * vectors are stored in the \f$n \times n\f$ matrix \a vectors.
  *
  * \see eigen(const Matrix<>& A, bool vectors)
  */
 
 struct Eigen {
   Matrix<> values;
   Matrix<> vectors;
 };

 /*!\brief Calculates the eigenvalues and eigenvectors of a symmetric
  * matrix.
  *
  * This function returns the eigenvalues and, optionally,
  * eigenvectors of a symmetric matrix \a A of order \f$n\f$.  It
  * returns an Eigen object containing the vector of values, in
  * ascending order, and, optionally, a matrix holding the vectors.
  *
  * \note This function requires BLAS/LAPACK functionality and is
  * only available on machines that provide these libraries.  Make
  * sure you enable the SCYTHE_LAPACK preprocessor flag if you wish
  * to use this function.  Furthermore, note that this function takes
  * and returns only column-major concrete matrices.  Future versions
  * of Scythe will provide a native C++ implementation of this
  * function with support for general matrix templates.
  *
  * \param A The Matrix to be decomposed.
  * \param vectors This boolean value indicates whether or not to
  * return eigenvectors in addition to eigenvalues.  It is set to true
  * by default.
  *
  * \throw scythe_null_error (Level 1)
  * \throw scythe_dimension_error (Level 1)
  * \throw scythe_lapack_internal_error (Level 1)
  *
  * \see Eigen
  * \see svd(const Matrix<>& A, int nu, int nv);
  */
  inline Eigen
  eigen (const Matrix<>& A, bool vectors=true)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for eigen");
    SCYTHE_CHECK_10(! A.isSquare(), scythe_dimension_error,
        "Matrix not square");
    SCYTHE_CHECK_10(A.isNull(), scythe_null_error,
        "Matrix is NULL");
    // Should be symmetric but rounding errors make checking for this
    // difficult.

    // Make a copy of A
    Matrix<> AA = A;

    // Get a point to the internal array and set up some vars
    double* Aarray = AA.getArray(); // internal array points
    int order = (int) AA.rows();    // input matrix is order x order
    double dignored = 0;            // we don't use this option
    int iignored = 0;               // or this one
    double abstol = 0.0;            // tolerance (default)
    int m;                          // output value
    Matrix<> result;                // result matrix
    char getvecs[1];                // are we getting eigenvectors?
    if (vectors) {
      getvecs[0] = 'V';
      result = Matrix<>(order, order + 1, false);
    } else {
      result = Matrix<>(order, 1, false);
      getvecs[0] = 'N';
    }
    double* eigenvalues = result.getArray(); // pointer to result array
    int* isuppz = new int[2 * order];        // indices of nonzero eigvecs
    double tmp;   // inital temporary value for getting work-space info
    int lwork, liwork, *iwork, itmp; // stuff for workspace
    double *work; // and more stuff for workspace
    int info = 0;  // error code holder

    // get optimal size for work arrays
    lwork = -1;
    liwork = -1;
    lapack::dsyevr_(getvecs, "A", "L", &order, Aarray, &order, &dignored,
        &dignored, &iignored, &iignored, &abstol, &m, eigenvalues, 
        eigenvalues + order, &order, isuppz, &tmp, &lwork, &itmp,
        &liwork, &info);
    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dsyevr");
    lwork = (int) tmp;
    liwork = itmp;
    work = new double[lwork];
    iwork = new int[liwork];

    // do the actual operation
    lapack::dsyevr_(getvecs, "A", "L", &order, Aarray, &order, &dignored,
        &dignored, &iignored, &iignored, &abstol, &m, eigenvalues, 
        eigenvalues + order, &order, isuppz, work, &lwork, iwork,
        &liwork, &info);
    SCYTHE_CHECK_10(info != 0, scythe_lapack_internal_error,
        "Internal error in LAPACK routine dsyevr");

    delete[] isuppz;
    delete[] work;
    delete[] iwork;

    Eigen resobj;
    if (vectors) {
      resobj.values = result(_, 0);
      resobj.vectors = result(0, 1, result.rows() -1, result.cols() - 1);
    } else {
      resobj.values = result;
    }

    return resobj;
  }

#endif

} // end namespace scythe

#endif /* SCYTHE_IDE_H */
