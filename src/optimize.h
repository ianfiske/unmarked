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
 *  scythestat/optimize.h
 *
 */

/*!
 * \file optimize.h
 * \brief Definitions of functions for doing numerical optimization
 * and related operations.
 *
 * This file contains a number of functions that are useful for
 * numerical optimization and maximum likelihood estimation.  In
 * addition, it contains some basic facilities for evaluating definite
 * integrals.
 *
 * As is the case across Scythe, we provide both general and default
 * template definitions for the functions in this file that return
 * Matrix objects.  The general definitions allow the user to
 * customize the matrix_order and matrix_style of the returned Matrix,
 * while the default versions return concrete matrices of the same
 * matrix_order as the first (or only) Matrix argument to the
 * function.  In cases where we supply these two types of definitions,
 * we explicitly document only the general version, although the
 * default definition will typically appear in the function list
 * below.
 *
 * \note
 * Doxygen has some difficulty dealing with overloaded templates.
 * Under certain circumstances it does not correctly process the
 * definitions of default templates.  In these cases, the definition
 * for the default template will not even appear in the function list.
 * We provide default templates for all of the Matrix-returning
 * functions in this file.
 *
 */

#ifndef SCYTHE_OPTIMIZE_H
#define SCYTHE_OPTIMIZE_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#include "algorithm.h"
#include "error.h"
#include "rng.h"
#include "distributions.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "stat.h"
#else
#include "scythestat/matrix.h"
#include "scythestat/algorithm.h"
#include "scythestat/error.h"
#include "scythestat/rng.h"
#include "scythestat/distributions.h"
#include "scythestat/la.h"
#include "scythestat/ide.h"
#include "scythestat/smath.h"
#include "scythestat/stat.h"
#endif

/* We want to use an anonymous namespace to make the following consts
 * and functions local to this file, but mingw doesn't play nice with
 * anonymous namespaces so we do things differently when using the
 * cross-compiler.
 */
#ifdef __MINGW32__
#define SCYTHE_MINGW32_STATIC static
#else
#define SCYTHE_MINGW32_STATIC
#endif

namespace scythe {
#ifndef __MINGW32__
  namespace {
#endif

    /* Functions (private to this file) that do very little... */
    template <typename T, matrix_order O, matrix_style S>
    SCYTHE_MINGW32_STATIC T donothing (const Matrix<T,O,S>& x)
    {
      return (T) 0.0;
    }

    template <typename T>
    SCYTHE_MINGW32_STATIC T donothing (T& x)
    {
      return (T) 0.0;
    }
#ifndef __MINGW32__
  }
#endif


  /* Return the machine epsilon 
   * Notes: Algorithm taken from Sedgewick, Robert. 1992. Algorithms
   * in C++. Addison Wesley. pg. 561
   */
   /*! \brief Compute the machine epsilon.
    *
    * The epsilon function returns the machine epsilon: the smallest
    * number that, when summed with 1, produces a value greater than
    * one.
    */
  template <typename T>
  T
  epsilon()
  {
    T eps, del, neweps;
    del    = (T) 0.5;
    eps    = (T) 0.0;
    neweps = (T) 1.0;
  
    while ( del > 0 ) {
      if ( 1 + neweps > 1 ) {  /* Then the value might be too large */
        eps = neweps;    /* ...save the current value... */
        neweps -= del;    /* ...and decrement a bit */
      } else {      /* Then the value is too small */
        neweps += del;    /* ...so increment it */
      }
      del *= 0.5;      /* Reduce the adjustment by half */
    }

    return eps;
  }
  
   /*! \brief Calculate the definite integral of a function from a to b.
    *
    * This function calculates the definite integral of a univariate
    * function on the interval \f$[a,b]\f$.
    *
    * \param fun The function (or functor) whose definite integral is
    * to be calculated.  This function should both take and return a
    * single argument of type T.
    * \param a The starting value of the interval.
    * \param b The ending value of the interval.
    * \param N The number of subintervals to calculate.  Increasing
    * this number will improve the accuracy of the estimate but will
    * also increase run-time.
    *
    * \throw scythe_invalid_arg (Level 1)
    *
    * \see adaptsimp(FUNCTOR fun, T a, T b, unsigned int& N, double tol = 1e-5)
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    */
  template <typename T, typename FUNCTOR>
  T  intsimp (FUNCTOR fun, T a, T b, unsigned int N)
  {
    SCYTHE_CHECK_10(a > b, scythe_invalid_arg,
        "Lower limit larger than upper");
    
    T I = (T) 0;
    T w = (b - a) / N;
    for (unsigned int i = 1; i <= N; ++i)
      I += w * (fun(a +(i - 1) *w) + 4 * fun(a - w / 2 + i * w) +
          fun(a + i * w)) / 6;
   
    return I;
  }
  
   /*! \brief Calculate the definite integral of a function from a to b.
    *
    * This function calculates the definite integral of a univariate
    * function on the interval \f$[a,b]\f$.
    *
    * \param fun The function (or functor) whose definite integral is 
    * to be calculated.  This function should both take and return a
    * single argument of type T.
    * \param a The starting value of the interval.
    * \param b The ending value of the interval.
    * \param N The number of subintervals to calculate.  Increasing
    * this number will improve the accuracy of the estimate but will
    * also increase run-time.
    * \param tol The accuracy required.  Both accuracy and run-time
    * decrease as this number increases.
    *
    * \throw scythe_invalid_arg (Level 1)
    *
    * \see intsimp(FUNCTOR fun, T a, T b, unsigned int& N)
    *
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    */
  template <typename T, typename FUNCTOR>
  T adaptsimp(FUNCTOR fun, T a, T b, unsigned int N, double tol = 1e-5)
  {
    SCYTHE_CHECK_10(a > b, scythe_invalid_arg,
        "Lower limit larger than upper");

    T I = intsimp(fun, a, b, N);
    if (std::fabs(I - intsimp(fun, a, b, N / 2)) > tol)
      return adaptsimp(fun, a, (a + b) / 2, N, tol)
        + adaptsimp(fun, (a + b) / 2, b, N, tol);

    return I;
  }

   /*! \brief Calculate gradient of a function using a forward
    * difference formula.
    *
    * This function numerically calculates the gradient of a
    * vector-valued function at \a theta using a forward difference
    * formula.
    *
    * \param fun The function to calculate the gradient of.  This
    * function should both take and return a single Matrix (vector) of 
    * type T.
    * \param theta The column vector of values at which to calculate 
    * the gradient of the function.
    *
    * \see gradfdifls(FUNCTOR fun, T alpha, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p)
    * \see jacfdif(FUNCTOR fun, const Matrix<T,PO,PS>& theta)
    * \see hesscdif(FUNCTOR fun, const Matrix<T,PO,PS>& theta)
    *
    * \throw scythe_dimension_error (Level 1)
    *
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS, typename FUNCTOR>
  Matrix<T, RO, RS>
  gradfdif (FUNCTOR fun, const Matrix<T,PO,PS>& theta)
      
  {
    SCYTHE_CHECK_10(! theta.isColVector(), scythe_dimension_error,
        "Theta not column vector");

    unsigned int k = theta.size();
    T h = std::sqrt(epsilon<T>());
    h = std::sqrt(h);

    Matrix<T,RO,RS> grad(k, 1);
    Matrix<T,RO> e;
    Matrix<T,RO> temp;
    for (unsigned int i = 0; i < k; ++i) {
      e = Matrix<T,RO>(k, 1);
      e[i] = h;
      temp = theta + e;
      donothing(temp); // XXX I don't understand this
      e = temp - theta;
      grad[i] = (fun(theta + e) - fun(theta)) / e[i];
    }

    return grad;
  }

  // Default template version
  template <typename T, matrix_order O, matrix_style S, 
            typename FUNCTOR>
  Matrix<T, O, Concrete>
  gradfdif (FUNCTOR fun, const Matrix<T,O,S>& theta)
  {
    return gradfdif<O,Concrete>(fun, theta);
  }

   /*! \brief Calculate the first derivative of the function using
    * a forward difference formula.
    *
    * This function numerically calculates the first derivative of a
    * function with respect to \a alpha at \f$theta + alpha \cdot p\f$
    * using a forward difference formula.  This function is primarily
    * useful for linesearches.
    *
    * \param fun The function to calculate the first derivative of.
    * This function should take a single Matrix<T> argument and return
    * a value of type T.
    * \param alpha Double the step length.
    * \param theta A Matrix (vector) of parameter values at which to
    * calculate the gradient.
    * \param p A direction vector.
    *
    * \see gradfdif(FUNCTOR fun, const Matrix<T,PO,PS>& theta)
    * \see jacfdif(FUNCTOR fun, const Matrix<T,PO,PS>& theta)
    * \see hesscdif(FUNCTOR fun, const Matrix<T,PO,PS>& theta)
    *
    * \throw scythe_dimension_error (Level 1)
    *
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    */
  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2, typename FUNCTOR>
  T
  gradfdifls (FUNCTOR fun, T alpha, const Matrix<T,PO1,PS1>& theta, 
              const Matrix<T,PO2,PS2>& p)

  {
    SCYTHE_CHECK_10(! theta.isColVector(), scythe_dimension_error,
        "Theta not column vector");
    SCYTHE_CHECK_10(! p.isColVector(), scythe_dimension_error,
        "p not column vector");

    unsigned int k = theta.size();
    T h = std::sqrt(epsilon<T>()); 
    h = std::sqrt(h);
    //T h = std::sqrt(2.2e-16);

    T deriv;

    for (unsigned int i = 0; i < k; ++i) {
      T temp = alpha + h;
      donothing(temp);
      T e = temp - alpha;
      deriv = (fun(theta + (alpha + e) * p) - fun(theta + alpha * p)) 
              / e;
    }
    
    return deriv;
  }

   /*! \brief Calculate the Jacobian of a function using a forward
    * difference formula.
    *
    * This function numerically calculates the Jacobian of a
    * vector-valued function using a forward difference formula.
    *
    * \param fun The function to calculate the Jacobian of.  This
    * function should both take and return a Matrix (vector) of type
    * T.
    * \param theta The column vector of parameter values at which to 
    * take the Jacobian of \a fun.
    *
    * \see gradfdif(FUNCTOR fun, const Matrix<T,PO,PS>& theta)
    * \see gradfdifls(FUNCTOR fun, T alpha, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p)
    * \see hesscdif(FUNCTOR fun, const Matrix<T,PO,PS>& theta)
    *
    * \throw scythe_dimension_error (Level 1)
    *
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS, typename FUNCTOR>
  Matrix<T,RO,RS>
  jacfdif (FUNCTOR fun, const Matrix<T,PO,PS>& theta)
  {
    SCYTHE_CHECK_10(! theta.isColVector(), scythe_dimension_error,
        "Theta not column vector");

    Matrix<T,RO> fval = fun(theta);
    unsigned int k = theta.rows();
    unsigned int n = fval.rows();

    T h = std::sqrt(epsilon<T>()); //2.2e-16
    h = std::sqrt(h);

    Matrix<T,RO,RS> J(n,k);
    Matrix<T,RO> e;
    Matrix<T,RO> temp;
    Matrix<T,RO> fthetae;
    Matrix<T,RO> ftheta;
    
    for (int i = 0; i < k; ++i) {
      e = Matrix<T,RO>(k,1);
      e[i] = h;
      temp = theta + e;
      donothing(temp); /// XXX ??
      e = temp - theta;
      fthetae = fun(theta + e);
      ftheta = fun(theta);
      for (unsigned int j = 0; j < n; ++j) {
        J(j,i) = (fthetae[j] - ftheta[j]) / e[i];
      }
    }
   
    return J;
  }

  // default template
  template <typename T, matrix_order PO, matrix_style PS,
            typename FUNCTOR>
  Matrix<T,PO,PS>
  jacfdif (FUNCTOR fun, const Matrix<T,PO,PS>& theta)
  {
    return jacfdif<PO,Concrete>(fun, theta);
  }


   /*! \brief Calculate the Hessian of a function using a central
    * difference formula.
    *
    * This function numerically calculates the Hessian of a
    * vector-valued function using a central difference formula.
    *
    * \param fun The function to calculate the Hessian of.  This
    * function should take a Matrix (vector) of type T and return a
    * single value of type T.
    * \param theta The column vector of parameter values at which to 
    * calculate the Hessian.
    *
    * \see gradfdif(FUNCTOR fun, const Matrix<T,PO,PS>& theta)
    * \see gradfdifls(FUNCTOR fun, T alpha, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p)
    * \see jacfdif(FUNCTOR fun, const Matrix<T,PO,PS>& theta)
    *
    * \throw scythe_dimension_error
    *
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS, typename FUNCTOR>
  Matrix<T, RO, RS>
  hesscdif (FUNCTOR fun, const Matrix<T,PO,PS>& theta)
  {
    SCYTHE_CHECK_10(! theta.isColVector(), scythe_dimension_error,
        "Theta not column vector");
    
    T fval = fun(theta);

    //std::cout << std::endl;
    //std::cout << "hesscdif theta = " << theta << "\n";
    //std::cout << "hesscdif fun(theta) = " << fval << std::endl;

    unsigned int k = theta.rows();

    // stepsize CAREFUL -- THIS IS MACHINE SPECIFIC !!!!
    T h2 = std::sqrt(epsilon<T>());
    //T h2 = (T) 1e-10;
    T h = std::sqrt(h2); 

    Matrix<T, RO, RS> H(k,k);

    //std::cout << "h2 = " << h2 << "  h = " << h << std::endl;

    Matrix<T,RO> ei;
    Matrix<T,RO> ej;
    Matrix<T,RO> temp;

    for (unsigned int i = 0; i < k; ++i) {
      ei = Matrix<T,RO>(k, 1);
      ei[i] = h;
      temp = theta + ei;
      donothing(temp); // XXX Again, I'm baffled
      ei = temp - theta;
      for (unsigned int j = 0; j < k; ++j){
        ej = Matrix<T,RO>(k,1);
        ej[j] = h;
        temp = theta + ej;
        donothing(temp); // XXX and again
        ej = temp - theta;
        
        if (i == j) {
          H(i,i) = ( -fun(theta + 2.0 * ei) + 16.0 *
              fun(theta + ei) - 30.0 * fval + 16.0 *
              fun(theta - ei) -
              fun(theta - 2.0 * ei)) / (12.0 * h2);
        } else {
          H(i,j) = ( fun(theta + ei + ej) - fun(theta + ei - ej)
              - fun(theta - ei + ej) + fun(theta - ei - ej))
            / (4.0 * h2);
        }
      }
    }
       
    //std::cout << "end of hesscdif, H = " << H << "\n";
    return H;
  }

  // default template
  template <typename T, matrix_order PO, matrix_style PS,
            typename FUNCTOR>
  Matrix<T,PO,PS>
  hesscdif (FUNCTOR fun, const Matrix<T,PO,PS>& theta)
  {
    return hesscdif<PO,Concrete>(fun, theta);
  }

   /*! \brief Find the step length that minimizes an implied 1-dimensional function.
    *
    * This function performs a line search to find the step length
    * that approximately minimizes an implied one dimensional
    * function.
    *
    * \param fun The function to minimize.  This function should take
    * one Matrix (vector) argument of type T and return a single value
    * of type T.
    * \param theta A column vector of parameter values that anchor the
    * 1-dimensional function.
    * \param p A direction vector that creates the 1-dimensional
    * function.
    *
    * \see linesearch2(FUNCTOR fun, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p, rng<RNGTYPE>& runif)
    * \see zoom(FUNCTOR fun, T alpha_lo, T alpha_hi, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p)
    * \see BFGS(FUNCTOR fun, const Matrix<T,PO,PS>& theta, rng<RNGTYPE>& runif, unsigned int maxit, T tolerance, bool trace = false)
    *
    * \throw scythe_dimension_error (Level 1)
    *
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    */
  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2, typename FUNCTOR>
  T linesearch1 (FUNCTOR fun, const Matrix<T,PO1,PS1>& theta,
                 const Matrix<T,PO2,PS2>& p)
  {
    SCYTHE_CHECK_10(! theta.isColVector(), scythe_dimension_error,
        "Theta not column vector");
    SCYTHE_CHECK_10(! p.isColVector(), scythe_dimension_error,
        "p not column vector");

    T alpha_bar = (T) 1.0;
    T rho = (T) 0.9;
    T c   = (T) 0.5;
    T alpha = alpha_bar;
    Matrix<T,PO1> fgrad = gradfdif(fun, theta);

    while (fun(theta + alpha * p) > 
           (fun(theta) + c * alpha * t(fgrad) * p)[0]) {
      alpha = rho * alpha;
    }

    return alpha;
  }

   /*! \brief Find the step length that minimizes an implied 1-dimensional function.
    *
    * This function performs a line search to find the step length
    * that approximately minimizes an implied one dimensional
    * function.
    *
    * \param fun The function to minimize.  This function should take
    * one Matrix (vector) argument of type T and return a single value
    * of type T.
    * \param theta A column vector of parameter values that anchor the
    * 1-dimensional function.
    * \param p A direction vector that creates the 1-dimensional
    * function.
    * \param runif A random uniform number generator function object
    * (an object that returns a random uniform variate on (0,1) when
    * its () operator is invoked).
    *
    * \see linesearch1(FUNCTOR fun, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p)
    * \see zoom(FUNCTOR fun, T alpha_lo, T alpha_hi, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p)
    * \see BFGS(FUNCTOR fun, const Matrix<T,PO,PS>& theta, rng<RNGTYPE>& runif, unsigned int maxit, T tolerance, bool trace = false)
    *
    * \throw scythe_dimension_error (Level 1)
    *
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    */
  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2, typename FUNCTOR,
            typename RNGTYPE>
  T linesearch2 (FUNCTOR fun, const Matrix<T,PO1,PS1>& theta,
                 const Matrix<T,PO2,PS2>& p, rng<RNGTYPE>& runif)
  {
    SCYTHE_CHECK_10(! theta.isColVector(), scythe_dimension_error,
        "Theta not column vector");
    SCYTHE_CHECK_10(! p.isColVector(), scythe_dimension_error,
        "p not column vector");

    T alpha_last = (T) 0.0;
    T alpha_cur = (T) 1.0;
    T alpha_max = (T) 10.0;
    T c1 = (T) 1e-4;
    T c2 = (T) 0.5;
    unsigned int max_iter = 50;
    T fgradalpha0 = gradfdifls(fun, (T) 0, theta, p);

    for (unsigned int i = 0; i < max_iter; ++i) {
      T phi_cur = fun(theta + alpha_cur * p);
      T phi_last = fun(theta + alpha_last * p);
     
      if ((phi_cur > (fun(theta) + c1 * alpha_cur * fgradalpha0))
          || ((phi_cur >= phi_last) && (i > 0))) {
        T alphastar = zoom(fun, alpha_last, alpha_cur, theta, p);
        return alphastar;
      }

      T fgradalpha_cur = gradfdifls(fun, alpha_cur, theta, p);
      if (std::fabs(fgradalpha_cur) <= -1 * c2 * fgradalpha0)
        return alpha_cur;

      if ( fgradalpha_cur >= (T) 0.0) {
        T alphastar = zoom(fun, alpha_cur, alpha_last, theta, p);
        return alphastar;
      }
      
      alpha_last = alpha_cur;
      // runif stuff below is probably not correc KQ 12/08/2006
      // I think it should work now DBP 01/02/2007
      alpha_cur = runif() * (alpha_max - alpha_cur) + alpha_cur;
    }

    return 0.001;
  }

   /*! \brief Find minimum of a function once bracketed.
    *
    * This function finds the minimum of a function, once bracketed.
    *
    * \param fun The function to minimize.  This function should take
    * one Matrix (vector) argument of type T and return a single value
    * of type T.
    * \param alpha_lo The lower bracket.
    * \param alpha_hi The upper bracket.
    * \param theta A column vector of parameter values that anchor the
    * 1-dimensional function.
    * \param p A direction vector that creates the 1-dimensional
    *
    * \see linesearch1(FUNCTOR fun, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p)
    * \see linesearch2(FUNCTOR fun, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p, rng<RNGTYPE>& runif)
    * \see BFGS(FUNCTOR fun, const Matrix<T,PO,PS>& theta, rng<RNGTYPE>& runif, unsigned int maxit, T tolerance, bool trace = false)
    *
    * \throw scythe_dimension_error (Level 1)
    *
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    * 
    */
  template <typename T, matrix_order PO1, matrix_style PS1,
            matrix_order PO2, matrix_style PS2, typename FUNCTOR>
  T zoom (FUNCTOR fun, T alpha_lo, T alpha_hi,
          const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p)
  {
    SCYTHE_CHECK_10(! theta.isColVector(), scythe_dimension_error,
        "Theta not column vector");
    SCYTHE_CHECK_10(! p.isColVector(), scythe_dimension_error,
        "p not column vector");

    T alpha_j = (alpha_lo + alpha_hi) / 2.0;
    T phi_0 = fun(theta);
    T c1 = (T) 1e-4;
    T c2 = (T) 0.5;
    T fgrad0 = gradfdifls(fun, (T) 0, theta, p);

    unsigned int count = 0;
    unsigned int maxit = 20;
    while(count < maxit) {
      T phi_j = fun(theta + alpha_j * p);
      T phi_lo = fun(theta + alpha_lo * p);
     
      if ((phi_j > (phi_0 + c1 * alpha_j * fgrad0))
          || (phi_j >= phi_lo)){
        alpha_hi = alpha_j;
      } else {
        T fgradj = gradfdifls(fun, alpha_j, theta, p);
        if (std::fabs(fgradj) <= -1 * c2 * fgrad0){ 
          return alpha_j;
        }
        if ( fgradj * (alpha_hi - alpha_lo) >= 0){
          alpha_hi = alpha_lo;
        }
        alpha_lo = alpha_j;
      }
      ++count;
    }
   
    return alpha_j;
  }


   /*! \brief Find function minimum using the BFGS algorithm.
    *
    * Numerically find the minimum of a function using the BFGS
    * algorithm.
    *
    * \param fun The function to minimize.  This function should take
    * one Matrix (vector) argument of type T and return a single value
    * of type T.
    * \param theta A column vector of parameter values that anchor the
    * 1-dimensional function.
    * \param runif A random uniform number generator function object
    * (an object that returns a random uniform variate on (0,1) when
    * its () operator is invoked).
    * \param maxit The maximum number of iterations.
    * \param tolerance The convergence tolerance.
    * \param trace Boolean value determining whether BFGS should print 
    *              to stdout (defaults to false).
    *
    * \see linesearch1(FUNCTOR fun, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p)
    * \see linesearch2(FUNCTOR fun, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p, rng<RNGTYPE>& runif)
    * \see zoom(FUNCTOR fun, T alpha_lo, T alpha_hi, const Matrix<T,PO1,PS1>& theta, const Matrix<T,PO2,PS2>& p)
    *
    * \throw scythe_dimension_error (Level 1)
    * \throw scythe_convergence_error (Level 0)
    *
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    */
  // there were 2 versions of linesearch1-- the latter was what we
  // had been calling linesearch2
  template <matrix_order RO, matrix_style RS, typename T, 
            matrix_order PO, matrix_style PS,
            typename FUNCTOR, typename RNGTYPE>
  Matrix<T,RO,RS>
  BFGS (FUNCTOR fun, const Matrix<T,PO,PS>& theta, rng<RNGTYPE>& runif, 
        unsigned int maxit, T tolerance, bool trace = false)
  {
    SCYTHE_CHECK_10(! theta.isColVector(), scythe_dimension_error,
        "Theta not column vector");

    unsigned int n = theta.size();

    // H is initial inverse hessian
    Matrix<T,RO> H = inv(hesscdif(fun, theta));

    // gradient at starting values
    Matrix<T,RO> fgrad = gradfdif(fun, theta);
    Matrix<T,RO> thetamin = theta;
    Matrix<T,RO> fgrad_new = fgrad;
    Matrix<T,RO> I = eye<T,RO>(n); 
    Matrix<T,RO> s;
    Matrix<T,RO> y;

    unsigned int count = 0;
    while( (t(fgrad_new)*fgrad_new)[0] > tolerance) {
      Matrix<T> p = -1.0 * H * fgrad;
      //std::cout << "initial H * fgrad = " << H * fgrad << "\n";
      //std::cout << "initial p = " << p << "\n";

      T alpha = linesearch2(fun, thetamin, p, runif);
      //T alpha = linesearch1(fun, thetamin, p);

      //std::cout << "after linesearch p = " << p << "\n";


      Matrix<T> thetamin_new = thetamin + alpha * p;
      fgrad_new = gradfdif(fun, thetamin_new);
      s = thetamin_new - thetamin;
      y = fgrad_new - fgrad;
      T rho = 1.0 / (t(y) * s)[0];
      H = (I - rho * s * t(y)) * H *(I - rho * y * t(s))
        + rho * s * t(s);

      thetamin = thetamin_new;
      fgrad = fgrad_new;
      ++count;

      if (trace){
	std::cout << "BFGS iteration = " << count << std::endl;
	std::cout << "thetamin = " << (t(thetamin)) ;
	std::cout << "gradient = " << (t(fgrad)) ;
	std::cout << "t(gradient) * gradient = " << 
	  (t(fgrad) * fgrad) ;
	std::cout << "function value = " << fun(thetamin) << 
	  std::endl << std::endl;
      }
      //std::cout << "Hessian = " << hesscdif(fun, theta) << "\n";
      //std::cout << "H = " << H << "\n";
      //std::cout << "alpha = " << alpha << std::endl;
      //std::cout << "p = " << p << "\n";
      //std::cout << "-1 * H * fgrad = " << -1.0 * H * fgrad << "\n";

      SCYTHE_CHECK(count > maxit, scythe_convergence_error,
          "Failed to converge.  Try better starting values");
    }
   
    return thetamin;
  }

  // Default template
  template <typename T, matrix_order O, matrix_style S,
            typename FUNCTOR, typename RNGTYPE>
  Matrix<T,O,Concrete>
  BFGS (FUNCTOR fun, const Matrix<T,O,S>& theta, rng<RNGTYPE>& runif,
        unsigned int maxit, T tolerance, bool trace = false)
  {
    return BFGS<O,Concrete> (fun, theta, runif, maxit, tolerance, trace);
  }

  
  /* Solves a system of n nonlinear equations in n unknowns of the form
   * fun(thetastar) = 0 for thetastar given the function, starting 
   * value theta, max number of iterations, and tolerance.
   * Uses Broyden's method.
   */
   /*! \brief Solve a system of nonlinear equations.
    *
    * Solves a system of n nonlinear equations in n unknowns of the form
    * \f$fun(\theta^*) = 0\f$ for \f$\theta^*\f$.
    *
    * \param fun The function to solve.  The function should both take
    * and return a Matrix of type T.
    * \param theta A column vector of parameter values at which to
    * start the solve procedure.
    * \param maxit The maximum number of iterations.
    * \param tolerance The convergence tolerance.
    *
    * \throw scythe_dimension_error (Level 1)
    * \throw scythe_convergence_error (Level 1)
    *
    * \note
    * Users will typically wish to implement \a fun in terms of a
    * functor.  Using a functor provides a generic way in which to
    * evaluate functions with more than one parameter.  Furthermore,
    * although one can pass a function pointer to this routine,
    * the compiler cannot inline and fully optimize code
    * referenced by function pointers.
    */
  template <matrix_order RO, matrix_style RS, typename T,
            matrix_order PO, matrix_style PS, typename FUNCTOR>
  Matrix<T,RO,RS>
  nls_broyden(FUNCTOR fun, const Matrix<T,PO,PS>& theta,
              unsigned int maxit = 5000, T tolerance = 1e-6)
  {
    SCYTHE_CHECK_10(! theta.isColVector(), scythe_dimension_error,
        "Theta not column vector");


    Matrix<T,RO> thetastar = theta;
    Matrix<T,RO> B = jacfdif(fun, thetastar);

    Matrix<T,RO> fthetastar;
    Matrix<T,RO> p;
    Matrix<T,RO> thetastar_new;
    Matrix<T,RO> fthetastar_new;
    Matrix<T,RO> s;
    Matrix<T,RO> y;

    for (unsigned int i = 0; i < maxit; ++i) {
      fthetastar = fun(thetastar);
      p = lu_solve(B, -1 * fthetastar);
      T alpha = (T) 1.0;
      thetastar_new = thetastar + alpha*p;
      fthetastar_new = fun(thetastar_new);
      s = thetastar_new - thetastar;
      y = fthetastar_new - fthetastar;
      B = B + ((y - B * s) * t(s)) / (t(s) * s);
      thetastar = thetastar_new;
      if (max(fabs(fthetastar_new)) < tolerance)
        return thetastar;
    }
 
    SCYTHE_THROW_10(scythe_convergence_error,  "Failed to converge.  Try better starting values or increase maxit");

    return thetastar;
  }

  // default template
  template <typename T, matrix_order O, matrix_style S,
            typename FUNCTOR>
  Matrix<T,O,Concrete>
  nls_broyden (FUNCTOR fun, const Matrix<T,O,S>& theta,
               unsigned int maxit = 5000, T tolerance = 1e-6)
  {
    return nls_broyden<O,Concrete>(fun, theta, maxit, tolerance);
  }

} // namespace scythe

#endif /* SCYTHE_OPTIMIZE_H */
