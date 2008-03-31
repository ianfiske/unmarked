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
 *  scythestat/la.h
 *
 */

/*!
 * \file la.h
 * \brief Definitions and implementations for functions that perform
 * common linear algebra manipulations on Scythe Matrix objects.
 *
 * This file provides a number of common linear algebraic functions
 * for use with the Matrix class.  These functions include common
 * operations such as transposition, a number of utility functions for
 * creating useful matrices like the identity matrix, and efficient
 * implementations for common operations like the cross-product.
 *
 * \note As is the case throughout the library, we provide both
 * general and default template definitions of the Matrix-returning
 * functions in this file, explicitly providing documentation for only
 * the general template versions.
 */

#ifndef SCYTHE_LA_H
#define SCYTHE_LA_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#include "algorithm.h"
#include "error.h"
#ifdef SCYTHE_LAPACK
#include "lapack.h"
#endif
#else
#include "scythestat/matrix.h"
#include "scythestat/algorithm.h"
#include "scythestat/error.h"
#ifdef SCYTHE_LAPACK
#include "scythestat/lapack.h"
#endif
#endif

#include <numeric>
#include <algorithm>
#include <set>

namespace scythe {

  namespace {
    typedef unsigned int uint;
  }

  /* Matrix transposition */

  /*!\brief Transpose a Matrix.
   *
   * This function transposes \a M, returning a Matrix \a R where each
   * element of \a M, \f$M_ij\f$ is placed in position \f$R_ji\f$.
   * Naturally, the returned Matrix has M.cols() rows and M.rows()
   * columns.
   *
   * \param M The Matrix to transpose.
   *
   * \throw scythe_alloc_error (Level 1)
   *
   */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  t (const Matrix<T,PO,PS>& M)
  {
    uint rows = M.rows();
    uint cols = M.cols();
    Matrix<T,RO,Concrete> ret(cols, rows, false);
    if (PO == Col)
      copy<Col,Row>(M, ret);
    else
      copy<Row,Col>(M, ret);

    SCYTHE_VIEW_RETURN(T, RO, RS, ret)
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  t (const Matrix<T,O,S>& M)
  {
    return t<O,Concrete>(M);
  }

  /* Ones matrix generation */
  
  /*! 
	 * \brief Create a matrix of ones.
	 *
	 * This function creates a matrix of ones, with the given dimensions
	 * \a rows and \a cols.
	 *
	 * \param rows The number of rows in the resulting Matrix.
	 * \param cols The number of columns in the resulting Matrix.
	 *
	 * \see eye (unsigned int k)
   *
   * \throw scythe_alloc_error (Level 1)
	 */
  template <typename T, matrix_order O, matrix_style S>
  Matrix<T,O,S>
  ones (unsigned int rows, unsigned int cols)
  {
    return Matrix<T,O,S> (rows, cols, true, (T) 1);
  }

  template <typename T, matrix_order O>
  Matrix<T, O, Concrete>
  ones (unsigned int rows, unsigned int cols)
  {
    return ones<T, O, Concrete>(rows, cols);
  }

  template <typename T>
  Matrix<T, Col, Concrete>
  ones (unsigned int rows, unsigned int cols)
  {
    return ones<T, Col, Concrete>(rows, cols);
  }

  inline Matrix<double, Col, Concrete>
  ones (unsigned int rows, unsigned int cols)
  {
    return ones<double, Col, Concrete>(rows, cols);
  }

  /* Identity  Matrix generation */

  // This functor contains the working parts of the eye algorithm.
  namespace {
    template <class T> struct eye_alg {
        T operator() (uint i, uint j) {
          if (i == j)
            return (T) 1.0;
          return (T) 0.0;
        }
    };
  }
  
  /*!\brief Create a \a k by \a k identity Matrix.
   *
   * This function creates a \a k by \a k Matrix with 1s along the
   * diagonal and 0s on the off-diagonal.  This template is overloaded
   * multiple times to provide default type, matrix_order, and
   * matrix_style.  The default call to eye returns a Concrete Matrix
   * containing double precision floating point numbers, in
   * column-major order.  The user can write explicit template calls
   * to generate matrices with other orders and/or styles.
   *
   * \param k The dimension of the identity Matrix.
   *
   * \see diag(const Matrix<T,O,S>& M)
   * \see ones(unsigned int rows, unsigned int cols)
   *
   * \throw scythe_alloc_error (Level 1)
   * 
   */
  template <typename T, matrix_order O, matrix_style S>
  Matrix<T,O,S>
  eye (unsigned int k)
  {
    Matrix<T,O,Concrete> ret(k, k, false);
    for_each_ij_set(ret, eye_alg<T>());
    SCYTHE_VIEW_RETURN(T, O, S, ret)
  }

  template <typename T, matrix_order O>
  Matrix<T, O, Concrete>
  eye (uint k)
  {
    return eye<T, O, Concrete>(k);
  }

  template <typename T>
  Matrix<T, Col, Concrete>
  eye (uint k)
  {
    return eye<T, Col, Concrete>(k);
  }

  inline Matrix<double, Col, Concrete>
  eye (uint k)
  {
    return eye<double, Col, Concrete>(k);
  }

  /* Create a k x 1 vector-additive sequence matrix */

  // The seqa algorithm
  namespace {
    template <typename T> struct seqa_alg {
      T cur_; T inc_;
      seqa_alg(T start, T inc) : cur_ (start), inc_ (inc)  {}
      T operator() () { T ret = cur_; cur_ += inc_; return ret; }
    };
  }
  
 /*! 
	* \brief Create a \a rows x 1 vector-additive sequence Matrix.
	*
	* This function creates a \a rows x 1 Matrix \f$v\f$, where 
  * \f$v_i = \mbox{start} + i \cdot \mbox{incr}\f$.
  *
  * This function is defined by a series of templates.  This template
  * is the most general, requiring the user to explicitly instantiate
  * the template in terms of element type, matrix_order and
  * matrix_style.  Further versions allow for explicit instantiation
  * based just on type and matrix_order (with matrix_style defaulting
  * to Concrete) and just on type (with matrix_style defaulting to
  * Col).  Finally, the default version of th function generates
  * column-major concrete Matrix of doubles.
	*
	* \param start Desired start value.
	* \param incr Amount to add in each step of the sequence.
	* \param rows Total number of rows in the Matrix.
  *
   * \throw scythe_alloc_error (Level 1)
	*/

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T,O,S>
  seqa (T start, T incr, uint rows)
  {
    Matrix<T,O,Concrete> ret(rows, 1, false);
    generate(ret.begin_f(), ret.end_f(), seqa_alg<T>(start, incr));
    SCYTHE_VIEW_RETURN(T, O, S, ret)
  }

  template <typename T, matrix_order O>
  Matrix<T, O, Concrete>
  seqa (T start, T incr, uint rows)
  {
    return seqa<T, O, Concrete>(start, incr, rows);
  }

  template <typename T>
  Matrix<T, Col, Concrete>
  seqa (T start, T incr, uint rows)
  {
    return seqa<T, Col, Concrete>(start, incr, rows);
  }

  inline Matrix<double, Col, Concrete>
  seqa (double start, double incr, uint rows)
  {
    return seqa<double, Col, Concrete>(start, incr, rows);
  }

  /* Uses the STL sort to sort a Matrix in ascending row-major order */
  
  /*! 
	 * \brief Sort a Matrix.
	 *
   * This function returns a copy of \a M, sorted in ascending order.
   * The sorting order is determined by the template parameter
   * SORT_ORDER or, by default, to matrix_order of \a M.
   *
	 * \param M The Matrix to sort.
   *
   * \see sortc
   * 
   * \throw scythe_alloc_error (Level 1)
	 */
  template <matrix_order SORT_ORDER,
            matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  sort (const Matrix<T, PO, PS>& M)
  {
    Matrix<T,RO,Concrete> ret = M;

    std::sort(ret.template begin<SORT_ORDER>(), 
              ret.template end<SORT_ORDER>());

    SCYTHE_VIEW_RETURN(T, RO, RS, ret)
  }

  template <matrix_order SORT_ORDER,
            typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  sort (const Matrix<T,O,S>& M)
  {
    return sort<SORT_ORDER, O, Concrete>(M);
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  sort (const Matrix<T,O,S>& M)
  {
    return sort<O, O, Concrete>(M);
  }

  /*!\brief Sort the columns of a Matrix.
   *
   * This function returns a copy of \a M, with each column sorted in
   * ascending order.
	 *
	 * \param M The Matrix to sort.
   * 
   * \see sort
   *
   * \throw scythe_alloc_error (Level 1)
	 */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  sortc (const Matrix<T, PO, PS>& M)
  {
    Matrix<T,RO,Concrete> ret = M;

    // TODO need to figure out a way to do fully optimized
    // vector iteration
    for (uint col = 0; col < ret.cols(); ++col) {
      Matrix<T,PO,View> column = ret(_, col);
      std::sort(column.begin(), column.end());
    }

    SCYTHE_VIEW_RETURN(T, RO, RS, ret)
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  sortc(const Matrix<T,O,S>& M)
  {
    return sortc<O,Concrete>(M);
  }

  /* Column bind two matrices */
  
  /*! 
	 * \brief Column bind two matrices.
	 *
	 * This function column binds two matrices, \a A and \a B.
	 *
	 * \param A The left-hand Matrix.
	 * \param B The right-hand Matrix.
	 *
	 * \see rbind(const Matrix<T,PO1,PS1>& A, 
   *            const Matrix<T,PO2,PS2>& B) 
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
	 */

  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,RO,RS>
  cbind (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& B)
  {
    SCYTHE_CHECK_10(A.rows() != B.rows(), scythe_conformation_error,
        "Matrices have different numbers of rows");

    Matrix<T,RO,Concrete> ret(A.rows(), A.cols() + B.cols(), false);
    std::copy(B.template begin_f<Col>(), B.template end_f<Col>(),
              std::copy(A.template begin_f<Col>(), 
                        A.template end_f<Col>(), 
                        ret.template begin_f<Col>()));
    SCYTHE_VIEW_RETURN(T, RO, RS, ret)
  }


  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,PO1,Concrete>
  cbind (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& B)
  {
    return cbind<PO1,Concrete>(A, B);
  }

  /* Row bind two matrices */
  
  /*! 
	 * \brief Row bind two matrices.
	 *
	 * This function row binds two matrices, \a A and \a B.
	 *
	 * \param A The upper Matrix.
	 * \param B The lower Matrix.
	 *
	 * \see cbind(const Matrix<T,PO1,PS1>& A, 
   *            const Matrix<T,PO2,PS2>& B) 
   *
   * \throw scythe_alloc_error (Level 1)
   * \throw scythe_conformation_error (Level 1)
	 */

  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,RO,RS>
  rbind (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& B)
  {
    SCYTHE_CHECK_10(A.cols() != B.cols(), scythe_conformation_error,
        "Matrices have different numbers of columns");

    Matrix<T,RO,Concrete> ret(A.rows() + B.rows(), A.cols(), false);
    std::copy(B.template begin_f<Row>(), B.template end_f<Row>(),
              std::copy(A.template begin_f<Row>(), 
                        A.template end_f<Row>(), 
                        ret.template begin_f<Row>()));
    SCYTHE_VIEW_RETURN(T, RO, RS, ret)
  }

  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,PO1,Concrete>
  rbind (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& B)
  {
    return rbind<PO1,Concrete>(A, B);
  }

  /* Calculates the order of each element in a Matrix */

  // Functor encapsulating the meat of the algorithm
  namespace {
    template <class T,matrix_order O,matrix_style S> struct order_alg {
      Matrix<T,O> M_;
      order_alg (const Matrix<T,O,S>& M) : M_ (M) {}
      uint operator() (T x) {
        Matrix<bool,O> diff = (M_ < x);
        return std::accumulate(diff.begin_f(), diff.end_f(), (uint) 0);
      }
    };
  }
    
  /*! 
	 * \brief Calculate the rank-order of each element in a Matrix.
	 *
	 * This function calculates the rank-order of each element in a
   * Matrix, returning a Matrix in which the \e i'th element
   * indicates the order position of the \e i'th element of \a M.
   * The returned Matrix contains unsigned integers.
	 * 
	 * \param M A column vector.
   *
   * \throw scythe_alloc_error (Level 1)
	 */

  /* NOTE This function used to only work on column vectors.  I see no
   * reason to maintain this restriction.
   */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<unsigned int, RO, RS>
  order (const Matrix<T, PO, PS>& M)
  {
    Matrix<uint, RO, Concrete> ranks(M.rows(), M.cols(), false);
    std::transform(M.begin_f(), M.end_f(), ranks.template begin_f<PO>(),
                   order_alg<T, PO, PS>(M));
    SCYTHE_VIEW_RETURN(uint, RO, RS, ranks)
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<unsigned int, O, Concrete>
  order (const Matrix<T,O,S>& M)
  {
    return order<O,Concrete>(M);
  }
  
  /* Selects all the rows of Matrix A for which binary column vector e
   * has an element equal to 1
   */
   
  /*! 
	 * \brief Locate rows for which a binary column vector equals 1
	  
   * This function identifies all the rows of a Matrix \a M for which
   * the binary column vector \a e has an element equal to 1,
   * returning a Matrix 
	 
	 * \param M The Matrix of interest.
	 * \param e A boolean column vector.
	 *
	 * \see unique(const Matrix<T>& M)
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_dimension_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
	 */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,RO,RS>
  selif (const Matrix<T,PO1,PS1>& M, const Matrix<bool,PO2,PS2>& e)
  {
    SCYTHE_CHECK_10(M.rows() != e.rows(), scythe_conformation_error,
     "Data matrix and selection vector have different number of rows");
    SCYTHE_CHECK_10(! e.isColVector(), scythe_dimension_error,
        "Selection matrix is not a column vector");

    uint N = std::accumulate(e.begin_f(), e.end_f(), (uint) 0);
    Matrix<T,RO,Concrete> res(N, M.cols(), false);
    int cnt = 0;
    for (uint i = 0; i < e.size(); ++i) {
      if (e[i]) {
        Matrix<T,RO,View> Mvec = M(i, _);
        // TODO again, need optimized vector iteration
        std::copy(Mvec.begin_f(), Mvec.end_f(), 
            res(cnt++, _).begin_f());
      }
    }

    SCYTHE_VIEW_RETURN(T, RO, RS, res)
  }
 
  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2>
  Matrix<T,PO1,Concrete>
  selif (const Matrix<T,PO1,PS1>& M, const Matrix<bool,PO2,PS2>& e)
  {
    return selif<PO1,Concrete>(M, e);
  }

  /* Find unique elements in a matrix and return a sorted row vector */
  /*! 
	 * \brief Find unique elements in a Matrix.
	 *
   * This function identifies all of the unique elements in a Matrix,
   * and returns them in a sorted row vector.
	 *
	 * \param M The Matrix to search.
	 *
	 * \see selif(const Matrix<T>& M, const Matrix<bool>& e)
   *
   * \throw scythe_alloc_error (Level 1)
	 */

  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  unique (const Matrix<T, PO, PS>& M)
  {
    std::set<T> u(M.begin_f(), M.end_f());
    Matrix<T,RO,Concrete> res(1, u.size(), false);
    std::copy(u.begin(), u.end(), res.begin_f());

    SCYTHE_VIEW_RETURN(T, RO, RS, res)
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  unique (const Matrix<T,O,S>& M)
  {
    return unique<O,Concrete>(M);
  }

  /* NOTE I killed reshape.  It seems redundant with resize. DBP */


  /* Make vector out of unique elements of a symmetric Matrix.  
   */
   
  /*! 
	 * \brief Vectorize a symmetric Matrix.
	 * 
   * This function returns a column vector containing only those
   * elements necessary to reconstruct the symmetric Matrix, \a M.  In
   * practice, this means extracting one triangle of \a M and
   * returning it as a vector.
   *
   * Note that the symmetry check in this function (active at error
   * level 3) is quite costly.
   *
	 * \param M A symmetric Matrix.
	 *
   * \throw scythe_dimension_error (Level 3)
   * \throw scythe_alloc_error (Level 1)
   *
   * \see xpnd(const Matrix<T,PO,PS>& v)
	 */

  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  vech (const Matrix<T, PO, PS>& M)
  {
    SCYTHE_CHECK_20(! M.isSymmetric(), scythe_dimension_error,
        "Matrix not symmetric");

    Matrix<T,RO,Concrete> 
      res((uint) (0.5 * (M.size() - M.rows())) + M.rows(), 1, false);
    typename Matrix<T,RO,Concrete>::forward_iterator it = res.begin_f();

    /* We want to traverse M in storage order if possible so we take
     * the upper triangle of row-order matrices and the lower triangle
     * of column-order matrices.
     */
    if (M.storeorder() == Col) {
      for (uint i = 0; i < M.rows(); ++i) {
        Matrix<T,PO,View> strip = M(i, i, M.rows() - 1, i);
        it = std::copy(strip.begin_f(), strip.end_f(), it);
      }
    } else {
      for (uint j = 0; j < M.cols(); ++j) {
        Matrix<T,PO,View> strip = M(j, j, j, M.cols() - 1);
        it = std::copy(strip.begin_f(), strip.end_f(), it);
      }
    }

    SCYTHE_VIEW_RETURN(T, RO, RS, res)
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  vech (const Matrix<T,O,S>& M)
  {
    return vech<O,Concrete>(M);
  }

  /*! Expand a vector into a symmetric Matrix.
   *
   * This function takes the vector \a v and returns a symmetric
   * Matrix containing the elements of \a v within each triangle.
   *
   * \param \a v The vector expand.
   *
   * \see vech(const Matrix<T,PO,PS>& M)
   *
   * \throw scythe_dimensions_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  xpnd (const Matrix<T, PO, PS>& v)
  {
    double size_d = -.5 + .5 * ::sqrt(1 + 8 * v.size());
    SCYTHE_CHECK_10(std::fmod(size_d, 1.) != 0., 
        scythe_dimension_error, 
        "Input vector can't generate square matrix");

    uint size = (uint) size_d;
    Matrix<T,RO,Concrete> res(size, size, false);

    /* It doesn't matter if we travel in order here.
     * TODO Might want to use iterators.
     */
    uint cnt = 0;
    for (uint i = 0; i < size; ++i)
      for (uint j = i; j < size; ++j)
        res(i, j) = res(j, i) = v[cnt++];

    SCYTHE_VIEW_RETURN(T, RO, RS, res)
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  xpnd (const Matrix<T,O,S>& v)
  {
    return xpnd<O,Concrete>(v);
  }

  /* Get the diagonal of a Matrix. */
  
  /*! 
	 * \brief Return the diagonal of a Matrix.
	 *
	 * This function returns the diagonal of a Matrix in a row vector.
	 *
	 * \param M The Matrix one wishes to extract the diagonal of.
	 *
	 * \see crossprod (const Matrix<T,PO,PS> &M)
   *
   * \throw scythe_alloc_error (Level 1)
	 */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  diag (const Matrix<T, PO, PS>& M)
  {
    Matrix<T,RO,Concrete> res(std::min(M.rows(), M.cols()), 1, false);
    
    /* We want to use iterators to maximize speed for both concretes
     * and views, but we always want to tranvers M in order to avoid
     * slowing down concretes.
     */
    uint incr = 1;
    if (PO == Col)
      incr += M.rows();
    else
      incr += M.cols();

    typename Matrix<T,PO,PS>::const_iterator pit;
    typename Matrix<T,RO,Concrete>::forward_iterator rit 
      = res.begin_f();
    for (pit = M.begin(); pit < M.end(); pit += incr)
      *rit++ = *pit;

    SCYTHE_VIEW_RETURN(T, RO, RS, res)
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  diag (const Matrix<T,O,S>& M)
  {
    return diag<O,Concrete>(M);
  }

  /* Fast calculation of A*B+C. */

  namespace {
    // Algorithm when one matrix is 1x1
    template <matrix_order RO, typename T,
              matrix_order PO1, matrix_style PS1,
              matrix_order PO2, matrix_style PS2>
    void
    gaxpy_alg(Matrix<T,RO,Concrete>& res, const Matrix<T,PO1,PS1>& X,
              const Matrix<T,PO2,PS2>& B, T constant)
    {
      res = Matrix<T,RO,Concrete>(X.rows(), X.cols(), false);
      if (maj_col<RO,PO1,PO2>())
        std::transform(X.template begin_f<Col>(), 
                       X.template end_f<Col>(),
                       B.template begin_f<Col>(),
                       res.template begin_f<Col>(),
                       ax_plus_b<T>(constant));
      else
        std::transform(X.template begin_f<Row>(), 
                       X.template end_f<Row>(),
                       B.template begin_f<Row>(),
                       res.template begin_f<Row>(),
                       ax_plus_b<T>(constant));
    }
  }

  /*! Fast caclulation of \f$AB + C\f$.
   *
   * This function calculates \f$AB + C\f$ efficiently, traversing the
   * matrices in storage order where possible, and avoiding the use of
   * extra temporary matrix objects.
   *
   * Matrices conform when \a A, \a B, and \a C are chosen with
   * dimensions
   * \f$((m \times n), (1 \times 1), (m \times n))\f$, 
   * \f$((1 \times 1), (n \times k), (n \times k))\f$, or
   * \f$((m \times n), (n \times k), (m \times k))\f$.
   *
   * Scythe will use LAPACK/BLAS routines to compute \f$AB+C\f$
   * with column-major matrices of double-precision floating point
   * numbers if LAPACK/BLAS is available and you compile your program
   * with the SCYTHE_LAPACK flag enabled.
   *
   * \param A A \f$1 \times 1\f$ or \f$m \times n\f$ Matrix.
   * \param B A \f$1 \times 1\f$ or \f$n \times k\f$ Matrix.
   * \param C A \f$m \times n\f$ or \f$n \times k\f$ or
   *          \f$m \times k\f$ Matrix.
   *
   * \throw scythe_conformation_error (Level 0)
   * \throw scythe_alloc_error (Level 1)
   */

  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2,
            matrix_order PO3, matrix_style PS3>
  Matrix<T,RO,RS>
  gaxpy (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& B,
         const Matrix<T,PO3,PS3>& C)
  {
    
    Matrix<T, RO, Concrete> res;

    if (A.isScalar() && B.rows() == C.rows() && B.cols() == C.cols()) {
      // Case 1: 1x1 * nXk + nXk
      gaxpy_alg(res, B, C, A[0]);
    } else if (B.isScalar() && A.rows() == C.rows() &&
               A.cols() == C.cols()) {
      // Case 2: m x n  *  1 x 1  +  m x n
      gaxpy_alg(res, A, C, B[0]);
    } else if (A.cols() == B.rows() && A.rows() == C.rows() &&
               B.cols() == C.cols()) {
      // Case 3: m x n  *  n x k  +  m x k

      res = Matrix<T,RO,Concrete> (A.rows(), B.cols(), false);

      /* These are identical to matrix mult, one optimized for
       * row-major and one for col-major.
       */
      
      T tmp;
      if (RO == Col) { // col-major optimized
       for (uint j = 0; j < B.cols(); ++j) {
         for (uint i = 0; i < A.rows(); ++i)
          res(i, j) = C(i, j);
         for (uint l = 0; l < A.cols(); ++l) {
           tmp = B(l, j);
           for (uint i = 0; i < A.rows(); ++i)
             res(i, j) += tmp * A(i, l);
         }
       }
      } else { // row-major optimized
       for (uint i = 0; i < A.rows(); ++i) {
         for (uint j = 0; j < B.cols(); ++j)
           res(i, j) = C(i, j);
         for (uint l = 0; l < B.rows(); ++l) {
           tmp = A(i, l);
           for (uint j = 0; j < B.cols(); ++j)
             res(i, j) += tmp * B(l,j);
         }
       }
      }

    } else {
      SCYTHE_THROW(scythe_conformation_error,
          "Expects (m x n  *  1 x 1  +  m x n)"
          << "or (1 x 1  *  n x k  +  n x k)"
          << "or (m x n  *  n x k  +  m x k)");
    }

    SCYTHE_VIEW_RETURN(T, RO, RS, res)
  }

  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2,
            matrix_order PO3, matrix_style PS3>
  Matrix<T,PO1,Concrete>
  gaxpy (const Matrix<T,PO1,PS1>& A, const Matrix<T,PO2,PS2>& B,
         const Matrix<T,PO3,PS3>& C)
  {
    return gaxpy<PO1,Concrete>(A,B,C);
  }

  /*! Fast caclulation of \f$A'A\f$.
   *
   * This function calculates \f$A'A\f$ efficiently, traversing the
   * matrices in storage order where possible, and avoiding the use of
   * the temporary matrix objects.
   *
   * Scythe will use LAPACK/BLAS routines to compute the cross-product
   * of column-major matrices of double-precision floating point
   * numbers if LAPACK/BLAS is available and you compile your program
   * with the SCYTHE_LAPACK flag enabled.
   *
   * \param A The Matrix to return the cross product of.
   *
   * \see diag (const Matrix<T, PO, PS>& M)
   */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS>
  Matrix<T, RO, RS>
  crossprod (const Matrix<T, PO, PS>& A)
  {
    /* When rows > 1, we provide differing implementations of the
     * algorithm depending on A's ordering to maximize strided access.
     *
     * The non-vector version of the algorithm fills in a triangle and
     * then copies it over.
     */
    Matrix<T, RO, Concrete> res;
    T tmp;
    if (A.rows() == 1) {
      res = Matrix<T,RO,Concrete>(A.cols(), A.cols(), true);
      for (uint k = 0; k < A.rows(); ++k) {
        for (uint i = 0; i < A.cols(); ++i) {
          tmp = A(k, i);
          for (uint j = i; j < A.cols(); ++j) {
            res(j, i) =
              res(i, j) += tmp * A(k, j);
          }
        }
      }
    } else {
      if (PO == Row) { // row-major optimized
        /* TODO: This is a little slower than the col-major.  Improve.
         */
        res = Matrix<T,RO,Concrete>(A.cols(), A.cols(), true);
        for (uint k = 0; k < A.rows(); ++k) {
          for (uint i = 0; i < A.cols(); ++i) {
            tmp = A(k, i);
            for (uint j = i; j < A.cols(); ++j) {
                res(i, j) += tmp * A(k, j);
            }
          }
        }
        for (uint i = 0; i < A.cols(); ++i)
          for (uint j = i + 1; j < A.cols(); ++j)
            res(j, i) = res(i, j);
      } else { // col-major optimized
        res = Matrix<T,RO,Concrete>(A.cols(), A.cols(), false);
        for (uint j = 0; j < A.cols(); ++j) {
          for (uint i = j; i < A.cols(); ++i) {
            tmp = (T) 0;
            for (uint k = 0; k < A.rows(); ++k)
              tmp += A(k, i) * A(k, j);
            res(i, j) = tmp;
          }
        }
        for (uint i = 0; i < A.cols(); ++i)
          for (uint j = i + 1; j < A.cols(); ++j)
            res(i, j) = res(j, i);
      }
    }

    SCYTHE_VIEW_RETURN(T, RO, RS, res)
  }

  template <typename T, matrix_order O, matrix_style S>
  Matrix<T, O, Concrete>
  crossprod (const Matrix<T,O,S>& M)
  {
    return crossprod<O,Concrete>(M);
  }

#ifdef SCYTHE_LAPACK
  /* Template specializations of for col-major, concrete
   * matrices of doubles that are only available when a lapack library
   * is available.
   */

  template<>
  Matrix<>
  gaxpy<Col,Concrete,double,Col,Concrete,Col,Concrete,Col,Concrete>
  (const Matrix<>& A, const Matrix<>& B, const Matrix<>& C)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for gaxpy");
    Matrix<> res;

    if (A.isScalar() && B.rows() == C.rows() && B.cols() == C.cols()) {
      // Case 1: 1x1 * nXk + nXk
      gaxpy_alg(res, B, C, A[0]);
    } else if (B.isScalar() && A.rows() == C.rows() &&
               A.cols() == C.cols()) {
      // Case 2: m x n  *  1 x 1  +  m x n
      gaxpy_alg(res, A, C, B[0]);
    } else if (A.cols() == B.rows() && A.rows() == C.rows() &&
               B.cols() == C.cols()) {
      res = C; // NOTE: this copy may eat up speed gains, but can't be
               //       avoided.
      
      // Case 3: m x n  *  n x k  +  m x k
      double* Apnt = A.getArray();
      double* Bpnt = B.getArray();
      double* respnt = res.getArray();
      const double one(1.0);
      int rows = (int) res.rows();
      int cols = (int) res.cols();
      int innerDim = A.cols();

      lapack::dgemm_("N", "N", &rows, &cols, &innerDim, &one, Apnt,
                     &rows, Bpnt, &innerDim, &one, respnt, &rows);

    }
      return res;
  }

  template<>
  Matrix<>
  crossprod(const Matrix<>& A)
  {
    SCYTHE_DEBUG_MSG("Using lapack/blas for crossprod");
    // Set up some constants
    const double zero = 0.0;
    const double one = 1.0;

    // Set up return value and arrays
    Matrix<> res(A.cols(), A.cols(), false);
    double* Apnt = A.getArray();
    double* respnt = res.getArray();
    int rows = (int) A.rows();
    int cols = (int) A.cols();

    lapack::dsyrk_("L", "T", &cols, &rows, &one, Apnt, &rows, &zero, respnt,
                   &cols);
    lapack::make_symmetric(respnt, cols);

    return res;
  }

#endif


} // end namespace scythe


#endif /* SCYTHE_LA_H */
