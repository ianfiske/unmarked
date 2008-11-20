/* 
 * Scythe Statistical Library
 * Copyright (C) 2000-2002 Andrew D. Martin and Kevin M. Quinn;
 * 2002-present Andrew D. Martin, Kevin M. Quinn, and Daniel
 * Pemstein.  All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * under the terms of the GNU General Public License as published by
 * Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.  See the text files COPYING
 * and LICENSE, distributed with this source code, for further
 * information.
 * --------------------------------------------------------------------
 * scythestat/algorithm.h
 */

/*!  \file algorithm.h 
 *
 * \brief Generic algorithms for Scythe objects.
 *
 * This file provides implementations of a few algorithms that operate
 * on Scythe objects and also contains the definitions of a handful of
 * useful function objects.  These functions and functors are primarily
 * intended for use within the library.  We add algorithms to this
 * header as need arises and do not currently attempt to provide a
 * comprehensive set of generic algorithms for working with Scythe
 * matrices.
 *
 */

#ifndef SCYTHE_ALGORITHM_H
#define SCYTHE_ALGORITHM_H

#include <cmath>
#include <functional>
#include <algorithm>

#ifdef SCYTHE_COMPILE_DIRECT
#include "defs.h"
#include "matrix.h"
#include "matrix_random_access_iterator.h"
#else
#include "scythestat/defs.h"
#include "scythestat/matrix.h"
#include "scythestat/matrix_random_access_iterator.h"
#endif

namespace scythe {
  namespace {
    typedef unsigned int uint;
  }

  /* Matrix forward declaration */
  template <typename T_type, matrix_order ORDER, matrix_style STYLE>
  class Matrix;

  /*! \brief A Functor encapsulating exponentiation.
   *
   * This function object wraps exponentiation operations for use in
   * generic algorithms.
   */
  template <typename T>
  struct exponentiate : std::binary_function<T, T, T>
  {
    T operator() (T base, T exp) const
    {
      return std::pow(base, exp);
    }
  };

  /*! \brief A Functor encapsulating \f$ax+b\f$.
   *
   * This function object wraps the operation \f$ax+b\f$ for use in
   * generic algorithms, where a is some constant.
   */
  template <typename T>
  struct ax_plus_b : std::binary_function<T,T,T>
  {
    T a_;
    ax_plus_b (T a) : a_ (a) {}
    T operator() (T x, T b) const
    {
      return (a_ * x + b);
    }
  };

  /*! \brief Iterate through a Matrix in order.
   *
   * This function iterates through a Matrix, \a M, in order,
   * setting each element in the Matrix to the result of an invocation
   * of the function object, \a func.  The () operator of \a func
   * should take two unsigned integer parameters (i - the row offset
   * into \a M; j - the column offset into \a M) and return a result
   * of type T.
   *
   * \param M The Matrix to iterate over.
   * \param func The functor to execute on each iteration.
   *
   */
   
  template <typename T, matrix_order O, matrix_style S, class FUNCTOR>
  void 
  for_each_ij_set (Matrix<T,O,S>& M, FUNCTOR func)
  {
    if (O == Col) {
      for (uint j = 0; j < M.cols(); ++j)
        for (uint i = 0; i < M.rows(); ++i)
          M(i, j) = func(i, j);
    } else {
      for (uint i = 0; i < M.cols(); ++i)
        for (uint j = 0; j < M.rows(); ++j)
          M(i, j) = func(i, j);
    }
  }

  /*! \brief Copy the contents of one Matrix into another.
   *
   * This function copies the contents of one Matrix into
   * another, traversing each Matrix in the order specified by the
   * template terms ORDER1 and ORDER2.  This function requires an
   * explicit template call that specifies ORDER1 and ORDER2.
   *
   * \param source The Matrix to copy.
   * \param dest   The Matrix to copy into.
   */

  template <matrix_order ORDER1, matrix_order ORDER2,
            typename T, typename S, matrix_order SO, matrix_style SS,
            matrix_order DO, matrix_style DS>
  void 
  copy(const Matrix<T,SO,SS>& source, Matrix<S,DO,DS>& dest)
  {
    std::copy(source.template begin_f<ORDER1>(), 
              source.template end_f<ORDER1>(),
              dest.template begin_f<ORDER2>());
  }

  /*! \brief Copy the contents of one Matrix into another.
   *
   * This function copies the contents of one Matrix into
   * another, traversing each Matrix in the order specified by the
   * template terms ORDER1 and ORDER2.  If \a source is larger than \a
   * dest, the function only copies as many elements from \a source as
   * will fit in \a dest.  On the other hand, if \a source is smaller
   * than \a dest, the function will start over at the beginning of
   * \a source, recycling the contents of \a source as many times as
   * necessary to fill \a dest.  This function requires an explicit
   * template call that specifies ORDER1 and ORDER2.
   *
   * \param source The Matrix to copy.
   * \param dest   The Matrix to copy into.
   */
  template <matrix_order ORDER1, matrix_order ORDER2,
            typename T, matrix_order SO, matrix_style SS,
            matrix_order DO, matrix_style DS>
  void 
  copy_recycle (const Matrix<T,SO,SS>& source, Matrix<T,DO,DS>& dest)
  {
    if (source.size() == dest.size()) {
      copy<ORDER1,ORDER2> (source, dest);
    } else if (source.size() > dest.size()) {
      const_matrix_random_access_iterator<T,ORDER1,SO,SS> s_iter 
        = source.template begin<ORDER1>();
      std::copy(s_iter, s_iter + dest.size(),
                dest.template begin_f<ORDER2>());
    } else {
      const_matrix_random_access_iterator<T,ORDER1,SO,SS> s_begin
        = source.template begin<ORDER1> ();
      matrix_random_access_iterator<T,ORDER2,DO,DS> d_iter 
        = dest.template begin<ORDER2>();
      matrix_random_access_iterator<T,ORDER2,DO,DS> d_end
        = dest.template end<ORDER2>();
      while (d_iter != d_end) {
        unsigned int span = std::min(source.size(), 
            (unsigned int) (d_end - d_iter));
        d_iter = std::copy(s_begin, s_begin + span, d_iter);
      }
    }
  }

  /*! \brief Determine the sign of a number.
   *
   * This function compares \a x to (T) 0, returning (T) 1 if \a x is
   * greater than zero, (T) -1 if \a x is less than zero, and (T) 0
   * otherwise.
   *
   * \param x The value to check.
   */
  template <class T>
  inline T sgn (const T & x)
  {
    if (x > (T) 0)
      return (T) 1;
    else if (x < (T) 0)
      return (T) -1;
    else
      return (T) 0;
  }

}  // end namespace scythe

#endif /* SCYTHE_ALGORITHM_H */
