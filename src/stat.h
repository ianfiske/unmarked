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
 *  scythestat/stat.h
 *
 */

/*!
 * \file stat.h
 * \brief Definitions for functions that perform common
 * statistical operations on Scythe Matrix objects.
 *
 * \note As is the case throughout the library, we provide both
 * general and default template definitions of the Matrix-returning
 * functions in this file, explicitly providing documentation for only
 * the general template versions.
 */

#ifndef SCYTHE_STAT_H
#define SCYTHE_STAT_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#include "algorithm.h"
#include "error.h"
#else
#include "scythestat/matrix.h"
#include "scythestat/algorithm.h"
#include "scythestat/error.h"
#endif

#include <numeric>
#include <set>


namespace scythe {

  namespace {
    typedef unsigned int uint;
  }

/* A macro for defining column versions of a function.  That is,
 * when expanded, this macro produces general and default template
 * functions that compute function NAME on each column in a matrix and
 * return a row vector with the results.  We use this to generate
 * column versions of every function in this header file.
 */
#define SCYTHE_STATMETH_COL(NAME)                                     \
  template <matrix_order RO, matrix_style RS, typename T,             \
            matrix_order PO, matrix_style PS>                         \
  Matrix<T,RO,RS>                                                     \
  NAME ## c (const Matrix<T,PO,PS>& A)                                \
  {                                                                   \
    Matrix<T,RO,RS> res (1, A.cols(), false);                         \
                                                                      \
    for (uint j = 0; j < A.cols(); ++j)                               \
      res[j] = NAME(A(_, j));                                         \
                                                                      \
    return res;                                                       \
  }                                                                   \
                                                                      \
  template <typename T, matrix_order O, matrix_style S>               \
  Matrix<T,O,Concrete>                                                \
  NAME ## c (const Matrix<T,O,S>& A)                                  \
  {                                                                   \
    return NAME ## c<O,Concrete>(A);                                  \
  }         


  /* Calculate the sum of a Matrix */
  
  /*! 
	* \brief Calculate the sum of a Matrix
	*
	* This function calculates the sum of a matrix by adding each element
	* in turn.
	*
	* \param A The matrix to be summed.
	*
	* \see prod(const Matrix<T,PO,PS> &A)
	* \see sumc(const Matrix<T,PO,PS> &A)
	* \see prodc(const Matrix<T,PO,PS> &A)
	*/


  template <typename T, matrix_order PO, matrix_style PS>
  T
  sum (const Matrix<T,PO,PS> &A)
  {
    return (std::accumulate(A.begin_f(), A.end_f(), (T) 0));
  }

  /* Calculate the sum of each column in a Matrix */
   
  /*! 
	* \brief Calculate the sum of each column in a Matrix
	*
	* This function calculates the sum of each column in a matrix by 
	* consecutively adding elements in a single column, looping through all 
	* columns, and returning the results.
	*
	* \param A The matrix to be summed.
	*
	* \see prod(const Matrix<T,PO,PS> &A)
	* \see sum(const Matrix<T,PO,PS> &A)
	* \see prodc(const Matrix<T,PO,PS> &A)
	*/

  SCYTHE_STATMETH_COL(sum);
  
  /* Calculate the product of a Matrix */
  
   /*! 
	* \brief Calculate the product of a Matrix
	*
	* This function calculates the product of a matrix by beginning with the 
	* first element of a matrix, and consecutively multiplying each entry.
	*
	* \param A The matrix to be multiplied.
	*
	* \see sumc(const Matrix<T,PO,PS> &A)
	* \see sum(const Matrix<T,PO,PS> &A)
	* \see prodc(const Matrix<T,PO,PS> &A)
	*/
	
  template <typename T, matrix_order PO, matrix_style PS>
  T
  prod (const Matrix<T,PO,PS> &A)
  {
    return std::accumulate(A.begin_f(), A.end_f(), (T) 1, 
                           std::multiplies<T> ());
  }

  /* Calculate the product of each column of a matrix */
  
   /*! 
	* \brief Calculate the product of each column of a Matrix
	*
	* This function calculates the product of each column of a matrix by 
	* multiplying all elements of a single column, looping through all columns,
	* and returning the results.
	*
	* \param A The matrix to be multiplied.
	*
	* \see sumc(const Matrix<T,PO,PS> &A)
	* \see sum(const Matrix<T,PO,PS> &A)
	* \see prod(const Matrix<T,PO,PS> &A)
	*/
	
  SCYTHE_STATMETH_COL(prod);
  
  /* Calculate the mean of a Matrix */
    
   /*! 
	* \brief Calculate the mean of a Matrix
	*
	* This function calculates the mean of a matrix by summing all elements of 
	* the matrix, and dividing by the total number of elements in the matrix.
	*
	* \param A The matrix to be averaged.
	*
	* \see sum(const Matrix<T,PO,PS> &A)
	* \see meanc(const Matrix<T,PO,PS> &A)
	* \see median(const Matrix<T,PO,PS> &A)
	* \see mode(const Matrix<T,PO,PS> &A)
	* \see variance(const Matrix<T,PO,PS> &A)
	*/

  template <typename T, matrix_order PO, matrix_style PS>
  T
  mean (const Matrix<T,PO,PS> &A)
  {
    return (std::accumulate(A.begin_f(), A.end_f(), (T) 0) / A.size());
  }

  /* Calculate the mean of each column of a Matrix */
  
   /*! 
	* \brief Calculate the mean of each column of a Matrix
	*
	* This function calculates the mean of each column of a matrix by summing 
	* all elements of a column in the matrix, divding by the total number of 
	* elements in the column, and looping over every column in the matrix.
	*
	* \param A The matrix to be averaged.
	*
	* \see sumc(const Matrix<T,PO,PS> &A)
	* \see mean(const Matrix<T,PO,PS> &A)
	* \see medianc(const Matrix<T,PO,PS> &A)
	* \see modec(const Matrix<T,PO,PS> &A)
	* \see variancec(const Matrix<T,PO,PS> &A)
	*/

  SCYTHE_STATMETH_COL(mean);
  
  /* Calculate the median of a matrix.  Uses a sort but I'll implement
   * the randomized alg when I figure out how to generalize it to
   * even-length lists
   */
   
   /*! 
	* \brief Calculate the median of a Matrix
	*
	* This function calculates the median of a matrix by first sorting the elements
	* of the matrix, and then finding the middle element.
	*
	* \param A The matrix whose median is of interest.
	*
	* \see medianc(const Matrix<T,PO,PS> &A)
	* \see mean(const Matrix<T,PO,PS> &A)
	* \see mode(const Matrix<T,PO,PS> &A)
	*/

  template <typename T, matrix_order PO, matrix_style PS>
  T
  median (const Matrix<T,PO,PS> &A)
  {
    Matrix<T, PO, PS> temp(A);
    uint n = temp.size();

    sort(temp.begin(), temp.end());
    if (n % 2 == 0)
      return ((temp[n / 2] + temp[n / 2 - 1]) / 2);
    else
      return temp[(uint) ::floor(n / 2)];
  }

  /* Calculate the median of each column of a matrix */
  
   /*! 
	* \brief Calculate the median of each column a Matrix
	*
	* This function calculates the median of each column of a matrix by first 
	* sorting the elements and locating the middle in a single column, and then
	* looping over all columns.
	*
	* \param A The matrix whose medians are of interest.
	*
	* \see median(const Matrix<T,PO,PS> &A)
	* \see meanc(const Matrix<T,PO,PS> &A)
	* \see modec(const Matrix<T,PO,PS> &A)
	*/

  SCYTHE_STATMETH_COL(median);

  /* Calculate the mode of a matrix */
  
   /*! 
	* \brief Calculate the mode of a Matrix
	*
	* This function calculates the mode of a matrix by determining which value of
	* the matrix occurs with the highest frequency.
	*
	* \param A The matrix whose mode is of interest.
	*
	* \see modec(const Matrix<T,PO,PS> &A)
	* \see mean(const Matrix<T,PO,PS> &A)
	* \see median(const Matrix<T,PO,PS> &A)
	*/
	
  template <typename T, matrix_order PO, matrix_style PS>
  T
  mode (const Matrix<T,PO,PS> &A)
  {
    Matrix<T, PO, PS> temp(A);
    
    sort(temp.begin(), temp.end());

    T last = temp[0];
    uint cnt = 1;
    T cur_max = temp[0];
    uint max_cnt = 1;
    
    for (uint i = 1; i < temp.size(); ++i) {
      if (last == temp[i]) {
        ++cnt;
      } else {
        last = temp[i];
        cnt = 1;
      }
      if (cnt > max_cnt) {
        max_cnt = cnt;
        cur_max = temp[i];
      }
    }

    return cur_max;
  }

   /*! 
	* \brief Calculate the mode of the columns of a Matrix
	*
	* This function calculates the mode of the columns of a matrix by 
	* determining which value in a single column of the matrix occurs 
	* most frequently, and then looping over all columns.
	*
	* \param A The matrix whose modes are of interest.
	*
	* \see mode(const Matrix<T,PO,PS> &A)
	* \see meanc(const Matrix<T,PO,PS> &A)
	* \see medianc(const Matrix<T,PO,PS> &A)
	*/

  SCYTHE_STATMETH_COL(mode);

  /* Calculate the variance of a Matrix */

  /* A functor that encapsulates a single variance calculation step.
   * Also used by skew and kurtosis. */
  namespace {
    template <typename T, typename T2>
    struct var_step : std::binary_function<T, T, T>
    {
      T constant_;
      T2 divisor_;
      T exponent_;
      var_step (T c, T2 d, T e) : constant_ (c), divisor_ (d),
                                    exponent_ (e) {}
      T operator() (T last, T x) const
      {
        return (last + std::pow(constant_ - x, exponent_) / divisor_);
      }
    };
  }
  
   /*! 
	* \brief Calculate the variance of a Matrix
	*
	* This function calculates the variance of a matrix.
	*
	* \param A The matrix whose variance is of interest.
	*
  * \see var(cons Matrix<T,PO,PS> &A, T mu)
	* \see varc(const Matrix<T,PO,PS> &A)
	* \see sd(const Matrix<T,PO,PS> &A)
	* \see mean(const Matrix<T,PO,PS> &A)
	*/
  template <typename T, matrix_order PO, matrix_style PS>
  T
  var (const Matrix<T,PO,PS> &A)
  {
    return var(A, mean(A));
  }

  /* Calculate the variances of each column of a Matrix. */
  
  /*! 
	* \brief Calculate the variance of each column of a Matrix
	*
	* This function calculates the variance of each column of a matrix.
	*
	* \param A The matrix whose variances are of interest.
	*
	* \see var(const Matrix<T,PO,PS> &A)
  * \see var(cons Matrix<T,PO,PS> &A, T mu)
	* \see sdc(const Matrix<T,PO,PS> &A)
	* \see meanc(const Matrix<T,PO,PS> &A)
	*/
	
  SCYTHE_STATMETH_COL(var);

 /*! 
	* \brief Calculate the variance of a Matrix
	*
	* This function calculates the variance of a matrix when the mean is
  * already known.
	*
	* \param A The matrix whose variance is of interest.
  * \param mu The mean of the values in the matrix.
	*
  * \see var(cons Matrix<T,PO,PS> &A)
	* \see varc(const Matrix<T,PO,PS> &A)
	* \see sd(const Matrix<T,PO,PS> &A)
	* \see mean(const Matrix<T,PO,PS> &A)
	*/
  template <typename T, matrix_order PO, matrix_style PS>
  T
  var (const Matrix<T,PO,PS> &A, T mu)
  {
    return std::accumulate(A.begin_f(), A.end_f(), (T) 0, 
                         var_step<T, uint> (mu, A.size() - 1, 2));
  }
  
  /* Calculate the standard deviation of a Matrix (not std cause of namespace std:: */
  
  /*! 
	* \brief Calculate the standard deviation of a Matrix
	*
	* This function calculates the standard deviation of a matrix by 
	* taking the square root of the matrix's variance. 	
	*
	* \param A The matrix whose standard deviation is of interest.
	*
  * \see sd(const Matrix<T,PO,PS) &A, T mu)
	* \see sdc(const Matrix<T,PO,PS> &A)
	* \see variance(const Matrix<T,PO,PS> &A)
	*/
	
  template <typename T, matrix_order PO, matrix_style PS>
  T
  sd (const Matrix<T,PO,PS> &A)
  {
    return std::sqrt(var(A));
  }
  
  /* Calculate the standard deviation of each column of a Matrix */
   /*! 
	* \brief Calculate the standard deviation of each column of a Matrix
	*
	* This function calculates the standard deviation of each column of a matrix by 
	* taking the square root of each column's variance. 	
	*
	* \param A The matrix whose standard deviations are of interest.
	*
	* \see sd(const Matrix<T,PO,PS> &A)
	* \see variancec(const Matrix<T,PO,PS> &A)
	*/
	
  SCYTHE_STATMETH_COL(sd);

  /*! 
	* \brief Calculate the standard deviation of a Matrix
	*
	* This function calculates the standard deviation of a matrix
  * when the matrix's mean is already known.
	*
	* \param A The matrix whose standard deviation is of interest.
  * \param mu The matrix mean.
	*
  * \see sd(const Matrix<T,PO,PS) &A)
	* \see sdc(const Matrix<T,PO,PS> &A)
	* \see variance(const Matrix<T,PO,PS> &A)
	*/
	
  template <typename T, matrix_order PO, matrix_style PS>
  T
  sd (const Matrix<T,PO,PS> &A, T mu)
  {
    return std::sqrt(var(A, mu));
  }

  /* Calculate the skew of a Matrix */

  /*! 
	 * \brief Calculate the skew of a Matrix
	 *
	 * This function calculates the skew of a matrix.
	 *
	 * \param A The matrix whose skew is of interest.
	 *
	 * \see skewc(const Matrix<T,PO,PS> &A)
	 * \see kurtosis(const Matrix<T,PO,PS> &A)
	 */

  template <typename T, matrix_order PO, matrix_style PS>
  T
  skew (const Matrix<T,PO,PS> &A)
  {
    T mu = mean(A);
    T sde = sd(A, mu);
    return std::accumulate(A.begin_f(), A.end_f(), (T) 0, 
              var_step<T, T> (mu, A.size() * std::pow(sde, 3), 3));
  }

  /* Calculate the skew of each column of a Matrix. */
  
   /*! 
	* \brief Calculate the skew of each column of a Matrix
	*
	* This function calculates the skew of each column of a matrix.
	*
	* \param A The matrix whose skews are of interest.
	*
	* \see skew(const Matrix<T,PO,PS> &A)
	* \see kurtosisc(const Matrix<T,PO,PS> &A)
	*/

  SCYTHE_STATMETH_COL(skew);
  
  /* Calculate the kurtosis of a Matrix */
    
   /*! 
	* \brief Calculate the kurtosis of a Matrix
	*
	* This function calculates the kurtosis of a matrix.
	*
	* \param A The matrix whose kurtosis is of interest.
	*
	* \see skew(const Matrix<T,PO,PS> &A)
	* \see kurtosisc(const Matrix<T,PO,PS> &A)
	*/
  template <typename T, matrix_order PO, matrix_style PS>
  T
  kurtosis (const Matrix<T,PO,PS> &A)
  {
    T mu = mean(A);
    T sde = sd(A, mu);
    return (std::accumulate(A.begin_f(), A.end_f(), (T) 0, 
              var_step<T, T> (mu, A.size() * std::pow(sde, 4), 4))
            - 3);
  }
  
  /* Calculate the kurtosis of each column of a Matrix. */
  
   /*! 
	* \brief Calculate the kurtosis of each column of a Matrix
	*
	* This function calculates the kurtosis of each column of a matrix.
	*
	* \param A The matrix whose kurtoses are of interest.
	*
	* \see skewc(const Matrix<T,PO,PS> &A)
	* \see kurtosis(const Matrix<T,PO,PS> &A)
	*/

  SCYTHE_STATMETH_COL(kurtosis);

  /* Calculates the maximum element in a Matrix */
  /*! 
	* \brief Calculate the maximum element in a Matrix
	*
	* This function identifies the maximum element in a matrix.
	*
	* \param A The matrix whose maximum element is of interest.
	*
	* \see min(const Matrix<T,PO,PS> &A)
	* \see maxc (const Matrix<T,PO,PS> &A)

	*/
	
  template <typename T, matrix_order PO, matrix_style PS>
  T
  max (const Matrix<T,PO,PS> &A)
  {
    return *(max_element(A.begin_f(), A.end_f()));
  }

   /*! 
	* \brief Calculate the maximum of each column of a Matrix
	*
	* This function identifies the maximum of each column in a matrix.
	*
	* \param A The matrix whose maximae are of interest.
	*
	* \see max(const Matrix<T,PO,PS> &A)
	* \see minc(const Matrix<T,PO,PS> &A)
	*/
  
  SCYTHE_STATMETH_COL(max);

  /* Calculates the minimum element in a Matrix */
  
	/*! 
	* \brief Calculate the maximum element in a Matrix
	*
	* This function identifies the maximum element in a matrix.
	*
	* \param A The matrix whose maximum element is of interest.
	*
	* \see max(const Matrix<T,PO,PS> &A)
	* \see minc(const Matrix<T,PO,PS> &A)
	*/
  template <typename T, matrix_order PO, matrix_style PS>
  T
  min (const Matrix<T,PO,PS> &A)
  {
    return *(min_element(A.begin_f(), A.end_f()));
  }
  
   /*! 
	* \brief Calculate the minimum of each column of a Matrix
	*
	* This function identifies the minimum of each column in a matrix.
	*
	* \param A The matrix whose minimae are of interest.
	*
	* \see min(const Matrix<T,PO,PS> &A)
	* \see maxc(const Matrix<T,PO,PS> &A)
	*/
  
  SCYTHE_STATMETH_COL(min);

  /* Find the index of the max element */  
	/*! 
	* \brief Calculate the index of the maximum element in a Matrix
	*
	* This function identifies the index of the maximum element in a matrix.
	*
	* \param A The matrix whose maximum element indices are of interest.
	*
	* \see minind(const Matrix<T,PO,PS> &A)
	* \see max(const Matrix<T,PO,PS> &A)
	* \see maxindc(const Matrix<T,PO,PS> &A)
	*/

  template <typename T, matrix_order PO, matrix_style PS>
  unsigned int
  maxind (const Matrix<T,PO,PS> &A)
  {
    return (max_element(A.begin_f(), A.end_f())).get_index();
  }
  
   /*! 
	* \brief Calculate the index of the maximum for each column of a Matrix
	*
	* This function identifies the index of the maximum for each column of a Matrix.
	*
	* \param A The matrix whose maximum indices are of interest.
	*
	* \see maxc(const Matrix<T,PO,PS> &A)
	* \see minindc(const Matrix<T,PO,PS> &A)
	*/

  SCYTHE_STATMETH_COL(maxind);
  
  /* Find the index of the min element */
  
  /*! 
	* \brief Calculate the index of the minimum element in a Matrix
	*
	* This function identifies the index of the minimum element in a matrix.
	*
	* \param A The matrix whose minimum element indices are of interest.
	*
	* \see maxind(const Matrix<T,PO,PS> &A)
	* \see min(const Matrix<T,PO,PS> &A)
	* \see minindc(const Matrix <T> &A)
	*/
  template <typename T, matrix_order PO, matrix_style PS>
  unsigned int
  minind (const Matrix<T,PO,PS> &A)
  {
    return (min_element(A.begin_f(), A.end_f())).get_index();
  }

    /*! 
	* \brief Calculate the index of the minimum for each column of a Matrix
	*
	* This function identifies the index of the minimum for each column of a Matrix.
	*
	* \param A The matrix whose minimum indices are of interest.
	*
	* \see minc(const Matrix<T,PO,PS> &A)
	* \see maxindc(const Matrix<T,PO,PS> &A)
	*/
  SCYTHE_STATMETH_COL(minind);

} // end namespace scythe


#endif /* SCYTHE_STAT_H */
