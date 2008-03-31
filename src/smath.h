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
 *  scythestat/smath.h
 *
 */

/*!
 * \file smath.h
 * \brief Definitions for functions that perform common mathematical
 * operations on every element of a Matrix.
 * 
 * \note As is the case throughout the library, we provide both
 * general and default template definitions of the Matrix-returning
 * functions in this file, explicitly providing documentation for only
 * the general template versions. As is also often the case, Doxygen
 * does not always correctly add the default template definition to
 * the function list below; there is always a default template
 * definition available for every function.
 *
 */

#ifndef SCYTHE_MATH_H
#define SCYTHE_MATH_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#include "algorithm.h"
#include "error.h"
#else
#include "scythestat/matrix.h"
#include "scythestat/algorithm.h"
#include "scythestat/error.h"
#endif

#include <cmath>
#include <numeric>
#include <set>

namespace scythe {

  namespace {
    typedef unsigned int uint;
  }

/* Almost every function in this file follows one of the two patterns
 * described by these macros.  The first macro handles single-argument
 * functions.  The second handles two-matrix-argument functions (or
 * scalar-matrix, matrix-scalar.  The second macro also permits
 * cross-type operations (these are limited only by the capabilities
 * of the underlying functions).
 */
#define SCYTHE_MATH_OP(NAME, OP)                                      \
  template <matrix_order RO, matrix_style RS, typename T,             \
            matrix_order PO, matrix_style PS>                         \
  Matrix<T,RO,RS>                                                     \
  NAME (const Matrix<T,PO,PS>& A)                                     \
  {                                                                   \
    Matrix<T,RO,RS> res(A.rows(), A.cols(), false);                   \
    std::transform(A.begin_f(), A.end_f(), res.begin_f(), OP);        \
    return res;                                                       \
  }                                                                   \
                                                                      \
  template <typename T, matrix_order O, matrix_style S>               \
  Matrix<T,O,Concrete>                                                \
  NAME (const Matrix<T,O,S>& A)                                       \
  {                                                                   \
    return NAME<O,Concrete>(A);                                       \
  }

#define SCYTHE_MATH_OP_2ARG(NAME, OP)                                 \
  template <matrix_order RO, matrix_style RS, typename T,             \
            matrix_order PO1, matrix_style PS1,                       \
            matrix_order PO2, matrix_style PS2, typename S>           \
  Matrix<T,RO,RS>                                                     \
  NAME (const Matrix<T,PO1,PS1>& A, const Matrix<S,PO2,PS2>& B)       \
  {                                                                   \
    SCYTHE_CHECK_10 (A.size() != 1 && B.size() != 1 &&                \
        A.size() != B.size(), scythe_conformation_error,              \
        "Matrices with dimensions (" << A.rows()                      \
        << ", " << A.cols()                                           \
        << ") and (" << B.rows() << ", " << B.cols()                  \
        << ") are not conformable");                                  \
                                                                      \
    Matrix<T,RO,RS> res;                                              \
                                                                      \
    if (A.size() == 1) {                                              \
      res.resize2Match(B);                                            \
      std::transform(B.template begin_f<RO>(), B.template end_f<RO>(),\
          res.begin_f(), std::bind1st(OP, A(0)));                     \
    } else if (B.size() == 1) {                                       \
      res.resize2Match(A);                                            \
      std::transform(A.template begin_f<RO>(), A.template end_f<RO>(),\
                     res.begin_f(), std::bind2nd(OP, B(0)));          \
    } else {                                                          \
      res.resize2Match(A);                                            \
      std::transform(A.template begin_f<RO>(), A.template end_f<RO>(),\
                     B.template begin_f<RO>(), res.begin_f(), OP);    \
    }                                                                 \
                                                                      \
    return res;                                                       \
  }                                                                   \
                                                                      \
  template <typename T, matrix_order PO1, matrix_style PS1,           \
                        matrix_order PO2, matrix_style PS2,           \
                        typename S>                                   \
  Matrix<T,PO1,Concrete>                                              \
  NAME (const Matrix<T,PO1,PS1>& A, const Matrix<S,PO2,PS2>& B)       \
  {                                                                   \
    return NAME<PO1,Concrete>(A, B);                                  \
  }                                                                   \
                                                                      \
  template<matrix_order RO, matrix_style RS, typename T,              \
           matrix_order PO, matrix_style PS, typename S>              \
  Matrix<T,RO,RS>                                                     \
  NAME (const Matrix<T,PO,PS>& A, S b)                                \
  {                                                                   \
    return NAME<RO,RS>(A, Matrix<S,RO,Concrete>(b));                  \
  }                                                                   \
                                                                      \
  template <typename T, typename S, matrix_order PO, matrix_style PS> \
  Matrix<T,PO,Concrete>                                               \
  NAME (const Matrix<T,PO,PS>& A, S b)                                \
  {                                                                   \
    return NAME<PO,Concrete>(A, Matrix<S,PO,Concrete>(b));            \
  }                                                                   \
                                                                      \
  template<matrix_order RO, matrix_style RS, typename T,              \
           matrix_order PO, matrix_style PS, typename S>              \
  Matrix<T,RO,RS>                                                     \
  NAME (T a, const Matrix<S,PO,PS>& B)                                \
  {                                                                   \
    return NAME<RO,RS>(Matrix<S, RO,Concrete>(a), B);                 \
  }                                                                   \
                                                                      \
  template <typename T, typename S, matrix_order PO, matrix_style PS> \
  Matrix<T,PO,Concrete>                                               \
  NAME (T a, const Matrix<S,PO,PS>& B)                                \
  {                                                                   \
    return NAME<PO,Concrete>(Matrix<S,PO,Concrete>(a), B);            \
  }


  /* calc the inverse cosine of each element of a Matrix */
  
 /*! 
	* \brief Calculate the inverse cosine of each element of a Matrix
	*
	* This function calculates the inverse cosine of each element in a Matrix
	*
	* \param A The matrix whose inverse cosines are of interest.
	*
	* \see tan()
	* \see tanh()
	* \see sin()
	* \see sinh()
	* \see cos()
	* \see cosh()
	* \see acosh()
	* \see asin()
	* \see asinh()
	* \see atan()
	* \see atanh()
	* \see atan2()
	*/

  SCYTHE_MATH_OP(acos, ::acos)
  
  /* calc the inverse hyperbolic cosine of each element of a Matrix */
   /*! 
	* \brief Calculate the inverse hyperbolic cosine of each element of a Matrix
	*
	* This function calculates the inverse hyperbolic cosine of each element
	* in a Matrix
	*
	* \param A The matrix whose inverse hyperbolic cosines are of interest.
	*
	* \see tan()
	* \see tanh()
	* \see sin()
	* \see sinh()
	* \see cos()
	* \see cosh()
	* \see acos()
	* \see asin()
	* \see asinh()
	* \see atan()
	* \see atanh()
	* \see atan2()
	*/

  SCYTHE_MATH_OP(acosh, ::acosh)

  /* calc the inverse sine of each element of a Matrix */
  
   /*! 
	* \brief Calculate the inverse sine of each element of a Matrix
	*
	* This function calculates the inverse sine of each element
	* in a Matrix
	*
	* \param A The matrix whose inverse sines are of interest.
	*
	* \see tan()
	* \see tanh()
	* \see sin()
	* \see sinh()
	* \see cos()
	* \see cosh()
	* \see acos()
	* \see acosh()
	* \see asinh()
	* \see atan()
	* \see atanh()
	* \see atan2()
	*/

  SCYTHE_MATH_OP(asin, ::asin)
  
  /* calc the inverse hyperbolic sine of each element of a Matrix */
  
  /*! 
	* \brief Calculate the inverse hyperbolic sine of each element of a Matrix
	*
	* This function calculates the inverse hyperbolic sine of each element
	* in a Matrix
	*
	* \param A The matrix whose inverse hyperbolic sines are of interest.
	*
	* \see tan()
	* \see tanh()
	* \see sin()
	* \see sinh()
	* \see cos()
	* \see cosh()
	* \see acos()
	* \see acosh()
	* \see asin()
	* \see atan()
	* \see atanh()
	* \see atan2()
	*/
	
  SCYTHE_MATH_OP(asinh, ::asinh)
  
  /* calc the inverse tangent of each element of a Matrix */
  
   /*! 
	* \brief Calculate the inverse tangent of each element of a Matrix
	*
	* This function calculates the inverse tangent of each element
	* in a Matrix
	*
	* \param A The matrix whose inverse tangents are of interest.
	*
	* \see tan()
	* \see tanh()
	* \see sin()
	* \see sinh()
	* \see cos()
	* \see cosh()
	* \see acos()
	* \see acosh()
	* \see asin()
	* \see asin()
	* \see atanh()
	* \see atan2()
	*/
	
  SCYTHE_MATH_OP(atan, ::atan)
  
  /* calc the inverse hyperbolic tangent of each element of a Matrix */
   /*! 
	* \brief Calculate the inverse hyperbolic tangent of each element of a Matrix
	*
	* This function calculates the inverse hyperbolic tangent of each element
	* in a Matrix
	*
	* \param A The matrix whose inverse hyperbolic tangents are of interest.
	*
	* \see tan()
	* \see tanh()
	* \see sin()
	* \see sinh()
	* \see cos()
	* \see cosh()
	* \see acos()
	* \see acosh()
	* \see asin()
	* \see asinh()
	* \see atan()
	* \see atan2()
	*/
	
  SCYTHE_MATH_OP(atanh, ::atanh)
  
  /* calc the angle whose tangent is y/x  */
  
   /*! 
	* \brief Calculate the angle whose tangent is y/x
	*
	* This function calculates the angle whose tangent is y/x, given two 
	* matrices A and B (where y is the ith element of A, and x is the jth element
	* of matrix B).
	*
	* \param A The matrix of y values 
	* \param B The matrix of x values
	*
	* \see tan()
	* \see tanh()
	* \see sin()
	* \see sinh()
	* \see cos()
	* \see cosh()
	* \see acos()
	* \see acosh()
	* \see asin()
	* \see asinh()
	* \see atan()
	* \see atanh()
	*/

  SCYTHE_MATH_OP_2ARG(atan2, std::ptr_fun(::atan2))

  /* calc the cube root of each element of a Matrix */
   /*! 
	* \brief Calculate the cube root of each element of a Matrix
	*
	* This function calculates the cube root of each element
	* in a Matrix
	*
	* \param A The matrix whose cube roots are of interest.
	*
	* \see sqrt()
	*/

  SCYTHE_MATH_OP(cbrt, ::cbrt)
  
  /* calc the ceil of each element of a Matrix */
  /*! 
	* \brief Calculate the ceiling of each element of a Matrix
	*
	* This function calculates the ceiling of each element
	* in a Matrix
	*
	* \param A The matrix whose ceilings are of interest.
	*
	* \see floor()
	*/

  SCYTHE_MATH_OP(ceil, ::ceil)
  
  /* create a matrix containing the absval of the first input and the
   * sign of the second
   */
    /*! 
	* \brief Create a matrix containing the absolute value of the first input
	* and the sign of the second input
	*
	* This function creates a matrix containing the absolute value of the first
	* input, a matrix called A, and the sign of the second input, matrix B.
	*
	* \param A The matrix whose absolute values will comprise the resultant matrix.
	* \param B The matrix whose signs will comprise the resultant matrix
	*/

  SCYTHE_MATH_OP_2ARG(copysign, std::ptr_fun(::copysign))
  
  /* calc the cosine of each element of a Matrix */
    
 /*! 
	* \brief Calculate the cosine of each element of a Matrix
	*
	* This function calculates the cosine of each element in a Matrix
	*
	* \param A The matrix whose cosines are of interest.
	*
	* \see tan()
	* \see tanh()
	* \see sin()
	* \see sinh()
	* \see cosh()
	* \see acos()
	* \see acosh()
	* \see asin()
	* \see asinh()
	* \see atan()
	* \see atanh()
	* \see atan2()
	*/

  SCYTHE_MATH_OP(cos, ::cos)
  
  /* calc the hyperbolic cosine of each element of a Matrix */
   /*! 
	* \brief Calculate the hyperbolic cosine of each element of a Matrix
	*
	* This function calculates the hyperbolic cosine of each element in a Matrix
	*
	* \param A The matrix whose hyperbolic cosines are of interest.
	*
	* \see tan()
	* \see tanh()
	* \see sin()
	* \see sinh()
	* \see cos()
	* \see acos()
	* \see acosh()
	* \see asin()
	* \see asinh()
	* \see atan()
	* \see atanh()
	* \see atan2()
	*/

  SCYTHE_MATH_OP(cosh, ::cosh)
  
  /* calc the error function of each element of a Matrix */
   /*! 
	* \brief Calculate the error function of each element of a Matrix
	*
	* This function calculates the error function of each element in a Matrix
	*
	* \param A The matrix whose error functions are of interest.
	*
	* \see erfc()
	*/

  SCYTHE_MATH_OP(erf, ::erf)
  
  /* calc the complementary error function of each element of a Matrix */
   /*! 
	* \brief Calculate the complementary error function of each element of a Matrix
	*
	* This function calculates the complemenatry error function of each 
	* element in a Matrix
	*
	* \param A The matrix whose complementary error functions are of interest.
	*
	* \see erf()
	*/

  SCYTHE_MATH_OP(erfc, ::erfc)
  
  /* calc the vaue e^x of each element of a Matrix */
   /*! 
	* \brief Calculate the value e^x for each element of a Matrix
	*
	* This function calculates the value e^x for each element of a matrix, where
	* x is the ith element of the matrix A
	*
	* \param A The matrix whose elements are to be exponentiated.
	*
	* \see expm1()
	*/

  SCYTHE_MATH_OP(exp, ::exp)
  
  /* calc the exponent - 1 of each element of a Matrix */
  /*! 
	* \brief Calculate the value e^(x-1) for each element of a Matrix
	*
	* This function calculates the value e^(x-1) for each element of a matrix, where
	* x is the ith element of the matrix A
	*
	* \param A The matrix whose elements are to be exponentiated.
	*
	* \see exp()
	*/

  SCYTHE_MATH_OP(expm1, ::expm1)
  
  /* calc the absval of each element of a Matrix */
   /*! 
	* \brief Calculate the absolute value of each element of a Matrix
	*
	* This function calculates the absolute value of each element in a Matrix
	*
	* \param A The matrix whose absolute values are to be taken.
	*/

  SCYTHE_MATH_OP(fabs, ::fabs)

  /* calc the floor of each element of a Matrix */
  /*! 
	* \brief Calculate the floor of each element of a Matrix
	*
	* This function calculates the floor of each element
	* in a Matrix
	*
	* \param A The matrix whose floors are of interest.
	*
	* \see ceil()
	*/

  SCYTHE_MATH_OP(floor, ::floor)
  
  /* calc the remainder of the division of each matrix element */
   /*! 
	* \brief Calculate the remainder of the division of each matrix element
	*
	* This function calculates the remainder when the elements of Matrix A are
	* divided by the elements of Matrix B.  
	*
	* \param A The matrix to serve as dividend
	* \param B the matrix to serve as divisor
	*/

  SCYTHE_MATH_OP_2ARG(fmod, std::ptr_fun(::fmod))

  /* calc the fractional val of input and return exponents in int
   * matrix reference
   */
   
   /*! 
	*/
  template <matrix_order RO, matrix_style RS, typename T,
	    matrix_order PO1, matrix_style PS1,
	    matrix_order PO2, matrix_style PS2>
  Matrix<T,RO,RS>
  frexp (const Matrix<T,PO1,PS1>& A, Matrix<int,PO2,PS2>& ex)
  {
    SCYTHE_CHECK_10(A.size() != ex.size(), scythe_conformation_error,
        "The input matrix sizes do not match");
    Matrix<T,PO1,Concrete> res(A.rows(), A.cols());
    
    typename Matrix<T,PO1,PS1>::const_forward_iterator it;
    typename Matrix<T,PO1,Concrete>::forward_iterator rit 
      = res.begin_f();
    typename Matrix<int,PO2,PS2>::const_forward_iterator it2
      = ex.begin_f();
    for (it = A.begin_f(); it != A.end_f(); ++it) {
      *rit = ::frexp(*it, &(*it2));
      ++it2; ++rit;
    }

    return res;
  }
  
  template <typename T, matrix_order PO1, matrix_style PS1,
	    matrix_order PO2, matrix_style PS2>
  Matrix<T,PO1,Concrete>
  frexp (Matrix<T,PO1,PS1>& A, Matrix<int,PO2,PS2>& ex)
  {
    return frexp<PO1,Concrete>(A,ex);
  }

  /* calc the euclidean distance between the two inputs */
  /*! 
	* \brief Calculate the euclidean distance between two inputs
	*
	* This function calculates the euclidean distance between the elements of Matrix
	* A and the elements of Matrix B.
	*
	* \param A Input matrix
	* \param B Input matrix
	*/

  SCYTHE_MATH_OP_2ARG(hypot, std::ptr_fun(::hypot))

  /*  return (int) logb */
  SCYTHE_MATH_OP(ilogb, ::ilogb)
  
  /* compute the bessel func of the first kind of the order 0 */
   /*! 
	* \brief Compute the Bessel function of the first kind of the order 0
	*
	* This function computes the Bessel function of the first kind of order 0
	* for each element in the input matrix, A.
	* 
	* \param A Matrix for which the Bessel function is of interest
	*
	* \see j1()
	* \see jn()
	* \see y0()
	* \see y1()
	* \see yn()
	*/
  
  SCYTHE_MATH_OP(j0, ::j0)
  
  /* compute the bessel func of the first kind of the order 1 */
  /*! 
	* \brief Compute the Bessel function of the first kind of the order 1
	*
	* This function computes the Bessel function of the first kind of order 1
	* for each element in the input matrix, A.
	* 
	* \param A Matrix for which the Bessel function is of interest
	*
	* \see j0()
	* \see jn()
	* \see y0()
	* \see y1()
	* \see yn()
	*/

  SCYTHE_MATH_OP(j1, ::j1)
  
  /* compute the bessel func of the first kind of the order n 
   * TODO: This definition causes the compiler to issue some warnings.
   * Fix
   */
   /*!
	* \brief Compute the Bessel function of the first kind of the order n
	*
	* This function computes the Bessel function of the first kind of order n
	* for each element in the input matrix, A.
	* 
	* \param n Order of the Bessel function
	* \param A Matrix for which the Bessel function is of interest
	*
	* \see j0()
	* \see j1()
	* \see y0()
	* \see y1()
	* \see yn()
	*/

  SCYTHE_MATH_OP_2ARG(jn, std::ptr_fun(::jn))

  /* calc x * 2 ^ex */
   /*!
	* \brief Compute x * 2^ex
	*
	* This function computes the value of x * 2^ex, where x is the ith element of
	* the input matrix A, and ex is the desired value of the exponent.
	* 
	* \param A Matrix whose elements are to be multiplied
	* \param ex Matrix of powers to which 2 will be raised.
	*/
  SCYTHE_MATH_OP_2ARG(ldexp, std::ptr_fun(::ldexp))
  
  /*  compute the natural log of the absval of gamma function */
  
   /*!
	* \brief Compute the natural log of the absolute value of the gamma function
	*
	* This function computes the absolute value of the Gamma Function, evaluated at
	* each element of the input matrix A.
	* 
	* \param A Matrix whose elements will serve as inputs for the Gamma Function
	*
	* \see log()
	*/

  SCYTHE_MATH_OP(lgamma, ::lgamma)
  
  /* calc the natural log of each element of a Matrix */
   /*!
	* \brief Compute the natural log of each element of a Matrix
	*
	* This function computes the natural log of each element in a matrix, A.
	* 
	* \param A Matrix whose natural logs are of interest
	*
	* \see log10()
	* \see log1p()
	* \see logb()
	*/

  SCYTHE_MATH_OP(log, ::log)
  
  /* calc the base-10 log of each element of a Matrix */
   /*!
	* \brief Compute the log base 10 of each element of a Matrix
	*
	* This function computes the log base 10 of each element in a matrix, A.
	* 
	* \param A Matrix whose logs are of interest
	*
	* \see log()
	* \see log1p()
	* \see logb()
	*/

  SCYTHE_MATH_OP(log10, ::log10)
  
  /* calc the natural log of 1 + each element of a Matrix */
  /*!
	* \brief Compute the natural log of 1 + each element of a Matrix
	*
	* This function computes the natural log of 1 + each element of a Matrix.
	* 
	* \param A Matrix whose logs are of interest
	*
	* \see log()
	* \see log10()
	* \see logb()
	*/
  
  SCYTHE_MATH_OP(log1p, ::log1p)
  
  /* calc the logb of each element of a Matrix */
  /*!
	* \brief Compute the logb each element of a Matrix
	*
	* This function computes the log base b of each element of a Matrix.
	* 
	* \param A Matrix whose logs are of interest
	*
	* \see log()
	* \see log10()
	* \see log1p()
	*/

  SCYTHE_MATH_OP(logb, ::logb)
  
  /* x = frac + i, return matrix of frac and place i in 2nd matrix
   */
  template <matrix_order RO, matrix_style RS, typename T,
	    matrix_order PO1, matrix_style PS1,
	    matrix_order PO2, matrix_style PS2>
  Matrix<T,RO,RS>
  modf (const Matrix<T,PO1,PS1>& A, Matrix<double,PO2,PS2>& ipart)
  {
    SCYTHE_CHECK_10(A.size() != ipart.size(), scythe_conformation_error,
        "The input matrix sizes do not match");
    Matrix<T,PO1,Concrete> res(A.rows(), A.cols());
    
    typename Matrix<T,PO1,PS1>::const_forward_iterator it;
    typename Matrix<T,PO1,Concrete>::forward_iterator rit 
      = res.begin_f();
    typename Matrix<double,PO2,PS2>::const_forward_iterator it2
      = ipart.begin_f();
    for (it = A.begin_f(); it != A.end_f(); ++it) {
      *rit = ::modf(*it, &(*it2));
      ++it2; ++rit;
    }

    return res;
  }
  
  template <typename T, matrix_order PO1, matrix_style PS1,
	    matrix_order PO2, matrix_style PS2>
  Matrix<T,PO1,Concrete>
  modf (Matrix<T,PO1,PS1>& A, Matrix<double,PO2,PS2>& ipart)
  {
    return modf<PO1,Concrete>(A,ipart);
  }

  /* calc x^ex of each element of a Matrix */
  
   /*!
	* \brief Compute x^ex for each element of a matrix
	*
	* This function computes x^ex, where x is the ith element of the matrix A, 
	* and ex is the desired exponent.
	* 
	* \param A Matrix to be exponentiated
	* \param ex Desired exponent
	*/
  SCYTHE_MATH_OP_2ARG(pow, std::ptr_fun(::pow))

  /* calc rem == x - n * y */
  SCYTHE_MATH_OP_2ARG(remainder, std::ptr_fun(::remainder))

  /* return x rounded to nearest int */
  
  /*!
	* \brief Return x rounded to the nearest integer
	*
	* This function returns x, where x is the ith element of the Matrix A, 
	* rounded to the nearest integer.
	* 
	* \param A Matrix whose elements are to be rounded
	*/

  SCYTHE_MATH_OP(rint, ::rint)

  /* returns x * FLT_RADIX^ex */
  SCYTHE_MATH_OP_2ARG(scalbn, std::ptr_fun(::scalbn))

  /*  calc the sine of x */
  
  /*! 
	* \brief Calculate the sine of each element of a Matrix
	*
	* This function calculates the sine of each element in a Matrix
	*
	* \param A The matrix whose sines are of interest.
	*
	* \see tan()
	* \see tanh()
	* \see sinh()
	* \see cos()
	* \see cosh()
	* \see acos()
	* \see acosh()
	* \see asin()
	* \see asinh()
	* \see atan()
	* \see atanh()
	* \see atan2()
	*/

  SCYTHE_MATH_OP(sin, ::sin)

  /* calc the hyperbolic sine of x */
   /*! 
	* \brief Calculate the hyperbolic sine of each element of a Matrix
	*
	* This function calculates the hyperbolic sine of each element in a Matrix
	*
	* \param A The matrix whose hyperbolic sines are of interest.
	*
	* \see tan()
	* \see tanh()
	* \see sin()
	* \see cos()
	* \see cosh()
	* \see acos()
	* \see acosh()
	* \see asin()
	* \see asinh()
	* \see atan()
	* \see atanh()
	* \see atan2()
	*/

  SCYTHE_MATH_OP(sinh, ::sinh)
  
  /* calc the sqrt of x */
  /*! 
	* \brief Calculate the square root of each element in a matrix
	*
	* This function calculates the square root of each element in a Matrix
	*
	* \param A The matrix whose roots are of interest.
	*
	* \see cbrt()

	*/
	
  SCYTHE_MATH_OP(sqrt, ::sqrt)

  /* calc the tangent of x */
  
  /*! 
	* \brief Calculate the tangent of each element of a Matrix
	*
	* This function calculates the tangent of each element in a Matrix
	*
	* \param A The matrix whose tangents are of interest.
	*
	* \see sinh()
	* \see tanh()
	* \see sin()
	* \see cos()
	* \see cosh()
	* \see acos()
	* \see acosh()
	* \see asin()
	* \see asinh()
	* \see atan()
	* \see atanh()
	* \see atan2()
	*/

  SCYTHE_MATH_OP(tan, ::tan)

  /* calc the hyperbolic tangent of x */
  /*! 
	* \brief Calculate the hyperbolic tangent of each element of a Matrix
	*
	* This function calculates the hyperbolic tangent of each element in a Matrix
	*
	* \param A The matrix whose hyperbolic tangents are of interest.
	*
	* \see sinh()
	* \see tan()
	* \see sin()
	* \see cos()
	* \see cosh()
	* \see acos()
	* \see acosh()
	* \see asin()
	* \see asinh()
	* \see atan()
	* \see atanh()
	* \see atan2()
	*/

  SCYTHE_MATH_OP(tanh, ::tanh)

  /* bessel function of the second kind of order 0*/
   /*! 
	* \brief Compute the Bessel function of the second kind of order 0
	*
	* This function computes the Bessel function of the second kind of order 0
	* for each element in the input matrix, A.
	* 
	* \param A Matrix for which the Bessel function is of interest
	*
	* \see j0()
	* \see j1()
	* \see jn()
	* \see y1()
	* \see yn()
	*/

  SCYTHE_MATH_OP(y0, ::y0)

  /* bessel function of the second kind of order 1*/
   /*! 
	* \brief Compute the Bessel function of the second kind of order 1
	*
	* This function computes the Bessel function of the second kind of order 1
	* for each element in the input matrix, A.
	* 
	* \param A Matrix for which the Bessel function is of interest
	*
	* \see j0()
	* \see j1()
	* \see jn()
	* \see y0()
	* \see yn()
	*/

  SCYTHE_MATH_OP(y1, ::y1)

  /* bessel function of the second kind of order n
   * TODO: This definition causes the compiler to issue some warnings.
   * Fix
   */
  /*!
	* \brief Compute the Bessel function of the second kind of order n
	*
	* This function computes the Bessel function of the second kind of order n
	* for each element in the input matrix, A.
	* 
	* \param n Order of the Bessel function
	* \param A Matrix for which the Bessel function is of interest
	*
	* \see j0()
	* \see j1()
	* \see jn()
	* \see y0()
	* \see y1()
	*/

  SCYTHE_MATH_OP_2ARG(yn, std::ptr_fun(::yn))
  
} // end namespace scythe

#endif /* SCYTHE_MATH_H */
