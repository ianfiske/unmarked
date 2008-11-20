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
 *  scythestat/rng.h
 *
 * The code for many of the RNGs defined in this file and implemented
 * in rng.cc is based on that in the R project, version 1.6.0-1.7.1.
 * This code is available under the terms of the GNU GPL.  Original
 * copyright:
 * 
 * Copyright (C) 1998      Ross Ihaka
 * Copyright (C) 2000-2002 The R Development Core Team
 * Copyright (C) 2003      The R Foundation
 */

/*!
 * \file rng.h
 *
 * \brief The definition of the random number generator base class.
 *
 */

/* Doxygen doesn't deal well with the macros that we use to make
 * matrix versions of rngs easy to define.
 */

#ifndef SCYTHE_RNG_H
#define SCYTHE_RNG_H

#include <iostream>
#include <cmath>

#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#include "error.h"
#include "algorithm.h"
#include "distributions.h"
#include "ide.h"
#include "la.h"
#else
#include "scythestat/matrix.h"
#include "scythestat/error.h"
#include "scythestat/algorithm.h"
#include "scythestat/distributions.h"
#include "scythestat/ide.h"
#include "scythestat/la.h"
#endif

namespace scythe {

/* Shorthand for the matrix versions of the various distributions'
 * random number generators.
 */
  
#define SCYTHE_RNGMETH_MATRIX(NAME, RTYPE, ARGNAMES, ...)             \
  template <matrix_order O, matrix_style S>                           \
  Matrix<RTYPE, O, S>                                                 \
  NAME (unsigned int rows, unsigned int cols, __VA_ARGS__)            \
  {                                                                   \
    Matrix<RTYPE, O, Concrete> ret(rows, cols, false);                \
    typename Matrix<RTYPE,O,Concrete>::forward_iterator it;           \
    typename Matrix<RTYPE,O,Concrete>::forward_iterator last          \
      = ret.end_f();                                                  \
    for (it = ret.begin_f(); it != last; ++it)                        \
      *it = NAME (ARGNAMES);                                          \
    SCYTHE_VIEW_RETURN(RTYPE, O, S, ret)                              \
  }                                                                   \
                                                                      \
  Matrix<RTYPE, Col, Concrete>                                        \
  NAME (unsigned int rows, unsigned int cols, __VA_ARGS__)            \
  {                                                                   \
    return NAME <Col,Concrete> (rows, cols, ARGNAMES);                \
  }

   /*! \brief Random number generator.
    *
    * This class provides objects capable of generating random numbers
    * from a variety of probability distributions.  This
    * abstract class forms the foundation of random number generation in
    * Scythe.  Specific random number generators should extend this class
    * and implement the virtual void function runif(); this function
    * should take no arguments and return uniformly distributed random
    * numbers on the interval (0, 1).  The rng class provides no
    * interface for seed-setting or initialization, allowing for maximal
    * flexibility in underlying implementation.  This class does provide
    * implementations of functions that return random numbers from a wide
    * variety of commonly (and not-so-commonly) used distributions, by
    * manipulating the uniform variates returned by runif().  See
    * rng/mersenne.h and rng/lecuyer.h for the rng implementations
    * offered by Scythe.
    *
    * Each univariate distribution is represented by three overloaded
    * versions of the same method.  The first is a simple method
    * returning a single value.  The remaining method versions return
    * Matrix values and are equivalent to calling the single-valued
    * method multiple times to fill a Matrix object.  They each take
    * two arguments describing the number of rows and columns in the
    * returned Matrix object and as many subsequent arguments as is
    * necessary to describe the distribution.  As is the case
    * throughout the library, the Matrix-returning versions of the
    * method include both a general and default template.  We
    * explicitly document only the single-valued versions of the
    * univariate methods.  For matrix-valued distributions we provide
    * only a single method per distribution.
    *
    * \note Doxygen incorrectly parses the macros we use to
    * automatically generate the Matrix returning versions of the
    * various univariate methods in this class.  Whenever you see the
    * macro variable __VA_ARGS__ in the public member function list
    * below, simply substitute in the arguments in the explicitly
    * documented single-valued version of the method.
    *
    */
  template <class RNGTYPE>
  class rng
  {
    public:

      /* This declaration allows users to treat rng objects like
       * functors that generate random uniform numbers.  This can be
       * quite convenient.
       */
       /*! \brief Generate uniformly distributed random variates.
        *
        * This operator acts as an alias for runif() and generates
        * pseudo-random variates from the uniform distribution on the
        * interval (0, 1).  We include this operator to allow rng
        * objects to behave as function objects.
        */
      double operator() ()
      {
        return runif();
      }

      /* Returns random uniform numbers on (0, 1).  This function must
       * be implemented by extending classes */
       /*! \brief Generate uniformly distributed random variates.
        *
        * This method generates pseudo-random variates from the
        * uniform distribution on the interval (0, 1).
        *
        * This function is pure virtual and is implemented by
        * extending concrete classes, like scythe::mersenne and
        * scythe::lecuyer.
        */
      double runif ()
      {
        return as_derived().runif();
      }


      /* No point in declaring these virtual because we have to
       * override them anyway because C++ isn't too bright.  Also, it
       * is illegal to make template methods virtual
       */
      template <matrix_order O, matrix_style S>
      Matrix<double,O,S> runif(unsigned int rows, 
                               unsigned int cols) 
      {
        Matrix<double, O, S> ret(rows, cols, false);
        typename Matrix<double,O,S>::forward_iterator it;
        typename Matrix<double,O,S>::forward_iterator last=ret.end_f();
        for (it = ret.begin_f(); it != last; ++it)
          *it = runif();

        return ret;
      }

      Matrix<double,Col,Concrete> runif(unsigned int rows, 
                                        unsigned int cols)
      {
        return runif<Col,Concrete>(rows, cols);
      }

      /*! \brief Generate a beta distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * beta distribution described by the shape parameters \a a and
			 * \a b.
       *
       * \param alpha The first positive beta shape parameter.
       * \param beta the second positive beta shape parameter.
			 * 
			 * \see pbeta(double x, double a, double b)
			 * \see dbeta(double x, double a, double b)
			 * \see betafn(double a, double b)
			 * \see lnbetafn(double a, double b)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rbeta (double alpha, double beta)
      {
        double report;
        double xalpha, xbeta;
        
        // Check for allowable parameters
        SCYTHE_CHECK_10(alpha <= 0, scythe_invalid_arg, "alpha <= 0");
        SCYTHE_CHECK_10(beta <= 0, scythe_invalid_arg, "beta <= 0");
        
        xalpha = rchisq (2 * alpha);
        xbeta = rchisq (2 * beta);
        report = xalpha / (xalpha + xbeta);
        
        return (report);
      }

      SCYTHE_RNGMETH_MATRIX(rbeta, double, SCYTHE_ARGSET(alpha, beta),
          double alpha, double beta);

      /*! \brief Generate a non-central hypergeometric disributed
			 * random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * non-centrial hypergeometric distribution described by the
       * number of positive outcomes \a m1, the two group size
       * parameters \a n1 and \a n2, and the odds ratio \a psi.
       *
       * \param m1 The number of positive outcomes in both groups.
       * \param n1 The size of group one.
       * \param n2 The size of group two.
       * \param psi The odds ratio
       * \param delta The precision.
       *
       * \throw scythe_convergence_error (Level 0)
       */
      double 
      rnchypgeom(double m1, double n1, double n2, double psi,
                 double delta)
      {
        // Calculate mode of mass function
        double a = psi - 1;
        double b = -1 * ((n1+m1+2)*psi + n2 - m1);
        double c = psi * (n1+1) * (m1+1);
        double q = -0.5 * ( b + sgn(b) * 
            std::sqrt(std::pow(b,2) - 4*a*c));
        double root1 = c/q;
        double root2 = q/a;
        double el = std::max(0.0, m1-n2);
        double u = std::min(n1,m1);
        double mode = std::floor(root1);
        int exactcheck = 0;
        if (u<mode || mode<el) {
          mode = std::floor(root2);
          exactcheck = 1;
        }
     

        int size = static_cast<int>(u+1);

        double *fvec = new double[size];
        fvec[static_cast<int>(mode)] = 1.0;
        double s;
        // compute the mass function at y
        if (delta <= 0 || exactcheck==1){  //exact evaluation 
          // sum from mode to u
          double f = 1.0;
          s = 1.0;
          for (double i=(mode+1); i<=u; ++i){
            double r = ((n1-i+1)*(m1-i+1))/(i*(n2-m1+i)) * psi;
            f = f*r;
            s += f;
            fvec[static_cast<int>(i)] = f;
          }
         
          // sum from mode to el
          f = 1.0;
          for (double i=(mode-1); i>=el; --i){
            double r = ((n1-i)*(m1-i))/((i+1)*(n2-m1+i+1)) * psi;
            f = f/r;
            s += f;
            fvec[static_cast<int>(i)] = f;
          }
        } else { // approximation
          double epsilon = delta/10.0;
          // sum from mode to ustar
          double f = 1.0;
          s = 1.0;
          double i = mode+1;
          double r;
          do {
            if (i>u) break;
            r = ((n1-i+1)*(m1-i+1))/(i*(n2-m1+i)) * psi;
            f = f*r;
            s += f;
            fvec[static_cast<int>(i)] = f;
            ++i;
          } while(f>=epsilon || r>=5.0/6.0);
         
          // sum from mode to elstar
          f = 1.0;
          i = mode-1;
          do {
            if (i<el) break;
            r = ((n1-i)*(m1-i))/((i+1)*(n2-m1+i+1)) * psi;
            f = f/r;
            s += f;
            fvec[static_cast<int>(i)] = f;
            --i;
          } while(f>=epsilon || r <=6.0/5.0);         
        }

        double udraw = runif();
        double psum = fvec[static_cast<int>(mode)]/s;
        if (udraw<=psum)
          return mode;
        double lower = mode-1;
        double upper = mode+1;

        do{
          double fl;
          double fu;
          if (lower >= el)
            fl = fvec[static_cast<int>(lower)];
          else 
            fl = 0.0;

          if (upper <= u)
            fu = fvec[static_cast<int>(upper)];
          else
            fu = 0.0;

          if (fl > fu) {
            psum += fl/s;
            if (udraw<=psum)
              return lower;
            --lower;
          } else {
            psum += fu/s;
            if (udraw<=psum)
              return upper;
            ++upper;
          }
        } while(udraw>psum);
       
        delete [] fvec;
        SCYTHE_THROW(scythe_convergence_error,
          "Algorithm did not converge");
      }

      SCYTHE_RNGMETH_MATRIX(rnchypgeom, double,
          SCYTHE_ARGSET(m1, n1, n2, psi, delta), double m1, double n1,
          double n2, double psi, double delta);

      /*! \brief Generate a Bernoulli distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * Bernoulli distribution with probability of success \a p.
       *
       * \param p The probability of success on a trial.
			 * 
       * \throw scythe_invalid_arg (Level 1)
       */
      unsigned int
      rbern (double p)
      {
        unsigned int report;
        double unif;
          
        // Check for allowable paramters
        SCYTHE_CHECK_10(p < 0 || p > 1, scythe_invalid_arg,
            "p parameter not in[0,1]");
        
        unif = runif ();
        if (unif < p)
          report = 1;
        else
          report = 0;
        
        return (report);
      }

      SCYTHE_RNGMETH_MATRIX(rbern, unsigned int, p, double p);
      
      /*! \brief Generate a binomial distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * binomial distribution with \a n trials and \p probability of
			 * success on each trial.
       *
       * \param n The number of trials.
       * \param p The probability of success on each trial.
			 * 
			 * \see pbinom(double x, unsigned int n, double p)
			 * \see dbinom(double x, unsigned int n, double p)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      unsigned int
      rbinom (unsigned int n, double p)
      {
        unsigned int report;
        unsigned int count = 0;
        double hold;
          
        // Check for allowable parameters
        SCYTHE_CHECK_10(n == 0, scythe_invalid_arg, "n == 0");
        SCYTHE_CHECK_10(p < 0 || p > 1, scythe_invalid_arg, 
            "p not in [0,1]");
          
        // Loop and count successes
        for (unsigned int i = 0; i < n; i++) {
          hold = runif ();
          if (hold < p)
            ++count;
        }
        report = count;
        
        return (report);
      }

      SCYTHE_RNGMETH_MATRIX(rbinom, unsigned int, SCYTHE_ARGSET(n, p),
          unsigned int n, double p);

      /*! \brief Generate a \f$\chi^2\f$ distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * \f$\chi^2\f$distribution with \a df degress of freedom.
       *
       * \param df The degrees of freedom.
			 * 
			 * \see pchisq(double x, double df)
			 * \see dchisq(double x, double df)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rchisq (double df)
      {
        double report;
          
        // Check for allowable paramter
        SCYTHE_CHECK_10(df <= 0, scythe_invalid_arg,
            "Degrees of freedom <= 0");
      
        // Return Gamma(nu/2, 1/2) variate
        report = rgamma (df / 2, .5);
        
        return (report);
      }

      SCYTHE_RNGMETH_MATRIX(rchisq, double, df, double df);

      /*! \brief Generate an exponentially distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * exponential distribution described by the inverse scale
			 * parameter \a invscale.
       *
       * \param invscale The inverse scale parameter.
			 * 
			 * \see pexp(double x, double scale)
			 * \see dexp(double x, double scale)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rexp (double invscale)
      {
        double report;
        
        // Check for allowable parameter
        SCYTHE_CHECK_10(invscale <= 0, scythe_invalid_arg,
            "Inverse scale parameter <= 0");
        
        report = -std::log (runif ()) / invscale;
        
        return (report);
      }

      SCYTHE_RNGMETH_MATRIX(rexp, double, invscale, double invscale);
    
      /*! \brief Generate an F distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * F distribution with degress of freedom \a df1 and \a df2.
       *
       * \param df1 The positive degrees of freedom for the
			 * \f$chi^2\f$ variate in the nominator of the F statistic.
       * \param df2 The positive degrees of freedom for the
			 * \f$chi^2\f$ variate in the denominator of the F statistic.
			 *
			 * \see pf(double x, double df1, double df2)
			 * \see df(double x, double df1, double df2)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rf (double df1, double df2)
      {
        SCYTHE_CHECK_10(df1 <= 0 || df2 <= 0, scythe_invalid_arg,
            "n1 or n2 <= 0");

        return ((rchisq(df1) / df1) / (rchisq(df2) / df2));
      }

      SCYTHE_RNGMETH_MATRIX(rf, double, SCYTHE_ARGSET(df1, df2),
          double df1, double df2);

      /*! \brief Generate a gamma distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * gamma distribution with a given \a shape and \a scale.
       *
       * \param shape The strictly positive shape of the distribution.
			 * \param rate The inverse of the strictly positive scale of the distribution.  That is, 1 / scale.
			 * 
			 * \see pgamma(double x, double shape, double scale)
			 * \see dgamma(double x, double shape, double scale)
			 * \see gammafn(double x)
			 * \see lngammafn(double x)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rgamma (double shape, double rate)
      {
        double report;

        // Check for allowable parameters
        SCYTHE_CHECK_10(shape <= 0, scythe_invalid_arg, "shape <= 0");
        SCYTHE_CHECK_10(rate <= 0, scythe_invalid_arg, "rate <= 0");

        if (shape > 1)
          report = rgamma1 (shape) / rate;
        else if (shape == 1)
          report = -std::log (runif ()) / rate;
        else if (shape < 1)
          report = rgamma1 (shape + 1) 
            * std::pow (runif (), 1 / shape) / rate;

        return (report);
      }

      SCYTHE_RNGMETH_MATRIX(rgamma, double, SCYTHE_ARGSET(shape, rate),
          double shape, double rate);

      /*! \brief Generate a logistically distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * logistic distribution described by the given \a location and
			 * \a scale variables.
       *
       * \param location The location of the distribution.
			 * \param scale The scale of the distribution.
			 * 
			 * \see plogis(double x, double location, double scale)
			 * \see dlogis(double x, double location, double scale)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rlogis (double location, double scale)
      {
        double report;
        double unif;
          
        // Check for allowable paramters
        SCYTHE_CHECK_10(scale <= 0, scythe_invalid_arg, "scale <= 0");
        
        unif = runif ();
        report = location + scale * std::log (unif / (1 - unif));
        
        return (report);
      }

      SCYTHE_RNGMETH_MATRIX(rlogis, double, 
          SCYTHE_ARGSET(location, scale),
          double location, double scale);

      /*! \brief Generate a log-normal distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * log-normal distribution with given logged mean and standard
			 * deviation.
       *
       * \param logmean The logged mean of the distribtion.
			 * \param logsd The strictly positive logged standard deviation
			 * of the distribution.
			 * 
			 * \see plnorm(double x, double logmean, double logsd)
			 * \see dlnorm(double x, double logmean, double logsd)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rlnorm (double logmean, double logsd)
      {
        SCYTHE_CHECK_10(logsd < 0.0, scythe_invalid_arg,
            "standard deviation < 0");

        return std::exp(rnorm(logmean, logsd));
      }

      SCYTHE_RNGMETH_MATRIX(rlnorm, double, 
          SCYTHE_ARGSET(logmean, logsd),
          double logmean, double logsd);

			/*! \brief Generate a negative binomial distributed random
			 * variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * negative binomial distribution with given dispersion
			 * parameter and probability of success on each trial.
       *
       * \param n The strictly positive target number of successful
			 * trials (dispersion parameters).
			 * \param p The probability of success on each trial.
			 * 
			 * \see pnbinom(unsigned int x, double n, double p)
			 * \see dnbinom(unsigned int x, double n, double p)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      unsigned int
      rnbinom (double n, double p)
      {
        SCYTHE_CHECK_10(n == 0 || p <= 0 || p > 1, scythe_invalid_arg,
            "n == 0, p <= 0, or p > 1");

        return rpois(rgamma(n, (1 - p) / p));
      }

      SCYTHE_RNGMETH_MATRIX(rnbinom, unsigned int,
          SCYTHE_ARGSET(n, p), double n, double p);

      /*! \brief Generate a normally distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * normal distribution with given \a mean and \a standard
			 * distribution.
       *
       * \param mean The mean of the distribution.
			 * \param sd The standard deviation of the distribution.
			 * 
			 * \see pnorm(double x, double mean, double sd)
			 * \see dnorm(double x, double mean, double sd)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rnorm (double mean = 0, double sd = 1)
      {
        SCYTHE_CHECK_10(sd <= 0, scythe_invalid_arg, 
            "Negative standard deviation");
        
        return (mean + rnorm1 () * sd);
      }

      SCYTHE_RNGMETH_MATRIX(rnorm, double, SCYTHE_ARGSET(mean, sd),
          double mean, double sd);

      /*! \brief Generate a Poisson distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * Poisson distribution with expected number of occurrences \a
			 * lambda.
       *
       * \param lambda The strictly positive expected number of
			 * occurrences.
			 * 
			 * \see ppois(double x, double lambda)
			 * \see dpois(double x, double lambda)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      unsigned int
      rpois(double lambda)
      {
        SCYTHE_CHECK_10(lambda <= 0, scythe_invalid_arg, "lambda <= 0");
        unsigned int n;
        
        if (lambda < 33) {
          double cutoff = std::exp(-lambda);
          n = -1;
          double t = 1.0;
          do {
            ++n;
            t *= runif();
          } while (t > cutoff);    
        } else {
          bool accept = false;
          double c = 0.767 - 3.36/lambda;
          double beta = M_PI/std::sqrt(3*lambda);
          double alpha = lambda*beta;
          double k = std::log(c) - lambda - std::log(beta);
            
          while (! accept){
            double u1 = runif();
            double x = (alpha - std::log((1-u1)/u1))/beta;
            while (x <= -0.5){
              u1 = runif();
              x = (alpha - std::log((1-u1)/u1))/beta;
            } 
            n = static_cast<int>(x + 0.5);
            double u2 = runif();
            double lhs = alpha - beta*x +
              std::log(u2/std::pow(1+std::exp(alpha-beta*x),2));
            double rhs = k + n*std::log(lambda) - lnfactorial(n);
            if (lhs <= rhs)
              accept = true;
          }
        }
        
        return n;
      }

      SCYTHE_RNGMETH_MATRIX(rpois, unsigned int, lambda, double lambda);

      /* There is a naming issue here, with respect to the p- and d-
       * functions in distributions.  This is really analagous to rt1-
       * and dt1- XXX Clear up.  Also, we should probably have a
       * random number generator for both versions of the student t.
       */

      /*! \brief Generate a Student t distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * Student's t distribution with given mean \a mu, variance \a
			 * sigma2, and degrees of freedom \a nu
       *
       * \param mu The mean of the distribution.
			 * \param sigma2 The variance of the distribution.
			 * \param nu The degrees of freedom of the distribution.
			 * 
			 * \see dt1(double x, double mu, double sigma2, double nu)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rt (double mu, double sigma2, double nu)
      {
        double report;
        double x, z;
          
        // Check for allowable paramters
        SCYTHE_CHECK_10(sigma2 <= 0, scythe_invalid_arg,
            "Variance parameter sigma2 <= 0");
        SCYTHE_CHECK_10(nu <= 0, scythe_invalid_arg,
            "D.O.F parameter nu <= 0");
        
        z = rnorm1 ();
        x = rchisq (nu);
        report = mu + std::sqrt (sigma2) * z 
          * std::sqrt (nu) / std::sqrt (x);
        
        return (report);
      }

      SCYTHE_RNGMETH_MATRIX(rt1, double, SCYTHE_ARGSET(mu, sigma2, nu),
          double mu, double sigma2, double nu);

      /*! \brief Generate a Weibull distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * Weibull distribution with given \a shape and \a scale.
       *
       * \param shape The strictly positive shape of the distribution.
			 * \param scale The strictly positive scale of the distribution.
			 * 
			 * \see pweibull(double x, double shape, double scale)
			 * \see dweibull(double x, double shape, double scale)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rweibull (double shape, double scale)
      {
        SCYTHE_CHECK_10(shape <= 0 || scale <= 0, scythe_invalid_arg,
            "shape or scale <= 0");

        return scale * std::pow(-std::log(runif()), 1.0 / shape);
      }

      SCYTHE_RNGMETH_MATRIX(rweibull, double,
          SCYTHE_ARGSET(shape, scale), double shape, double scale);

			/*! \brief Generate an inverse \f$\chi^2\f$ distributed random
			 * variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * inverse \f$\chi^2\f$ distribution with \a nu degress of
			 * freedom.
       *
       * \param nu The degrees of freedom.
			 *
			 * \see rchisq(double df)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      richisq (double nu)
      {
        double report;
          
        // Check for allowable parameter
        SCYTHE_CHECK_10(nu <= 0, scythe_invalid_arg,
            "Degrees of freedom <= 0");
          
        // Return Inverse-Gamma(nu/2, 1/2) variate
        report = rigamma (nu / 2, .5);
        return (report);
      }

      SCYTHE_RNGMETH_MATRIX(richisq, double, nu, double nu);

      /*! \brief Generate an inverse gamma distributed random variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * inverse gamma distribution with given \a shape and \a scale.
       *
       * \param shape The strictly positive shape of the distribution.
			 * \param scale The strictly positive scale of the distribution.
			 * 
			 * \see rgamma(double alpha, double beta)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rigamma (double alpha, double beta)
      {
        double report;
        
        // Check for allowable parameters
        SCYTHE_CHECK_10(alpha <= 0, scythe_invalid_arg, "alpha <= 0");
        SCYTHE_CHECK_10(beta <= 0, scythe_invalid_arg, "beta <= 0");

        // Return reciprocal of gamma variate
        report = std::pow (rgamma (alpha, beta), -1);

        return (report);
      }

      SCYTHE_RNGMETH_MATRIX(rigamma, double, SCYTHE_ARGSET(alpha, beta),
          double alpha, double beta);

      /* Truncated Distributions */

			/*! \brief Generate a truncated normally distributed random
			 * variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * normal distribution with given \a mean and \a variance,
			 * truncated both above and below.  It uses the inverse CDF
			 * method.
       *
       * \param mean The mean of the distribution.
			 * \param variance The variance of the distribution.
			 * \param below The lower truncation point of the distribution.
			 * \param above The upper truncation point of the distribution.
			 * 
			 * \see rtnorm_combo(double mean, double variance, double below, double above)
			 * \see rtbnorm_slice(double mean, double variance, double below, unsigned int iter = 10)
			 * \see rtanorm_slice(double mean, double variance, double above, unsigned int iter = 10)
			 * \see rtbnorm_combo(double mean, double variance, double below, unsigned int iter = 10)
			 * \see rtanorm_combo(double mean, double variance, double above, unsigned int iter = 10)
			 * \see rnorm(double x, double mean, double sd)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double 
      rtnorm(double mean, double variance, double below, double above)
      {  
        SCYTHE_CHECK_10(below >= above, scythe_invalid_arg,
            "Truncation bound not logically consistent");
        SCYTHE_CHECK_10(variance <= 0, scythe_invalid_arg,
            "Variance <= 0");
        
        double sd = std::sqrt(variance);
        double FA = 0.0;
        double FB = 0.0;
        if ((std::fabs((above-mean)/sd) < 8.2) 
            && (std::fabs((below-mean)/sd) < 8.2)){
          FA = pnorm1((above-mean)/sd, true, false);
          FB = pnorm1((below-mean)/sd, true, false);
        }
        if ((((above-mean)/sd) < 8.2)  && (((below-mean)/sd) <= -8.2) ){ 
          FA = pnorm1((above-mean)/sd, true, false);
          FB = 0.0;
        }
        if ( (((above-mean)/sd) >= 8.2)  && (((below-mean)/sd) > -8.2) ){ 
          FA = 1.0;
          FB = pnorm1((below-mean)/sd, true, false);
        } 
        if ( (((above-mean)/sd) >= 8.2) && (((below-mean)/sd) <= -8.2)){
          FA = 1.0;
          FB = 0.0;
        }
        double term = runif()*(FA-FB)+FB;
        if (term < 5.6e-17)
          term = 5.6e-17;
        if (term > (1 - 5.6e-17))
          term = 1 - 5.6e-17;
        double draw = mean + sd * qnorm1(term);
        if (draw > above)
          draw = above;
        if (draw < below)
          draw = below;
         
        return draw;
      }

      SCYTHE_RNGMETH_MATRIX(rtnorm, double, 
          SCYTHE_ARGSET(mean, variance, above, below), double mean,
          double variance, double above, double below);

			/*! \brief Generate a truncated normally distributed random
			 * variate.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * normal distribution with given \a mean and \a variance,
			 * truncated both above and below.  It uses a combination of
			 * rejection sampling (when \a below <= mean <= \a above)
			 * sampling method of Robert and Casella (1999), pp. 288-289
			 * (when \a meam < \a below or \a mean > \a above).
       *
       * \param mean The mean of the distribution.
			 * \param variance The variance of the distribution.
			 * \param below The lower truncation point of the distribution.
			 * \param above The upper truncation point of the distribution.
			 * 
			 * \see rtnorm(double mean, double variance, double below, double above)
			 * \see rtbnorm_slice(double mean, double variance, double below, unsigned int iter = 10)
			 * \see rtanorm_slice(double mean, double variance, double above, unsigned int iter = 10)
			 * \see rtbnorm_combo(double mean, double variance, double below, unsigned int iter = 10)
			 * \see rtanorm_combo(double mean, double variance, double above, unsigned int iter = 10)
			 * \see rnorm(double x, double mean, double sd)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rtnorm_combo(double mean, double variance, double below, 
                   double above)
      {
        SCYTHE_CHECK_10(below >= above, scythe_invalid_arg,
            "Truncation bound not logically consistent");
        SCYTHE_CHECK_10(variance <= 0, scythe_invalid_arg,
            "Variance <= 0");
        
        double sd = std::sqrt(variance);
        if ((((above-mean)/sd > 0.5) && ((mean-below)/sd > 0.5))
            ||
            (((above-mean)/sd > 2.0) && ((below-mean)/sd < 0.25))
            ||
            (((mean-below)/sd > 2.0) && ((above-mean)/sd > -0.25))) { 
          double x = rnorm(mean, sd);
          while ((x > above) || (x < below))
            x = rnorm(mean,sd);
          return x;
        } else {
          // use the inverse cdf method
          double FA = 0.0;
          double FB = 0.0;
          if ((std::fabs((above-mean)/sd) < 8.2) 
              && (std::fabs((below-mean)/sd) < 8.2)){
            FA = pnorm1((above-mean)/sd, true, false);
            FB = pnorm1((below-mean)/sd, true, false);
          }
          if ((((above-mean)/sd) < 8.2)  && (((below-mean)/sd) <= -8.2) ){ 
            FA = pnorm1((above-mean)/sd, true, false);
            FB = 0.0;
          }
          if ( (((above-mean)/sd) >= 8.2)  && (((below-mean)/sd) > -8.2) ){ 
            FA = 1.0;
            FB = pnorm1((below-mean)/sd, true, false);
          } 
          if ( (((above-mean)/sd) >= 8.2) && (((below-mean)/sd) <= -8.2)){
            FA = 1.0;
            FB = 0.0;
          }
          double term = runif()*(FA-FB)+FB;
          if (term < 5.6e-17)
            term = 5.6e-17;
          if (term > (1 - 5.6e-17))
            term = 1 - 5.6e-17;
          double x = mean + sd * qnorm1(term);
          if (x > above)
            x = above;
          if (x < below)
            x = below;
          return x;
        }    
      }

      SCYTHE_RNGMETH_MATRIX(rtnorm_combo, double, 
          SCYTHE_ARGSET(mean, variance, above, below), double mean,
          double variance, double above, double below);

			/*! \brief Generate a normally distributed random variate,
			 * truncated below.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * normal distribution with given \a mean and \a variance,
			 * truncated below.  It uses the slice sampling method of
			 * Robert and Casella (1999), pp. 288-289.
       *
			 * \param mean The mean of the distribution.
			 * \param variance The variance of the distribution.
			 * \param below The lower truncation point of the distribution.
			 * \param iter The number of iterations to use.
			 * 
			 * \see rtnorm(double mean, double variance, double below, double above)
			 * \see rtnorm_combo(double mean, double variance, double below, double above)
			 * \see rtanorm_slice(double mean, double variance, double above, unsigned int iter = 10)
			 * \see rtbnorm_combo(double mean, double variance, double below, unsigned int iter = 10)
			 * \see rtanorm_combo(double mean, double variance, double above, unsigned int iter = 10)
			 * \see rnorm(double x, double mean, double sd)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rtbnorm_slice (double mean, double variance, double below,
                     unsigned int iter = 10)
      {
        SCYTHE_CHECK_10(below < mean, scythe_invalid_arg,
            "Truncation point < mean");
        SCYTHE_CHECK_10(variance <= 0, scythe_invalid_arg,
            "Variance <= 0");
         
        double z = 0;
        double x = below + .00001;
         
        for (unsigned int i=0; i<iter; ++i){
          z = runif()*std::exp(-1*std::pow((x-mean),2)/(2*variance));
          x = runif()*
            ((mean + std::sqrt(-2*variance*std::log(z))) - below) + below;
        }

        if (! finite(x)) {
          SCYTHE_WARN("Mean extremely far from truncation point. "
              << "Returning truncation point");
          return below; 
        }

        return x;
      }

      SCYTHE_RNGMETH_MATRIX(rtbnorm_slice, double, 
          SCYTHE_ARGSET(mean, variance, below, iter), double mean, 
          double variance, double below, unsigned int iter = 10);

			/*! \brief Generate a normally distributed random variate,
			 * truncated above.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * normal distribution with given \a mean and \a variance,
			 * truncated above.  It uses the slice sampling method of Robert
			 * and Casella (1999), pp. 288-289.
       *
       * \param mean The mean of the distribution.
			 * \param variance The variance of the distribution.
			 * \param above The upper truncation point of the distribution.
			 * \param iter The number of iterations to use.
			 * 
			 * \see rtnorm(double mean, double variance, double below, double above)
			 * \see rtnorm_combo(double mean, double variance, double below, double above)
			 * \see rtbnorm_slice(double mean, double variance, double below, unsigned int iter = 10)
			 * \see rtbnorm_combo(double mean, double variance, double below, unsigned int iter = 10)
			 * \see rtanorm_combo(double mean, double variance, double above, unsigned int iter = 10)
			 * \see rnorm(double x, double mean, double sd)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rtanorm_slice (double mean, double variance, double above, 
          unsigned int iter = 10)
      {
        SCYTHE_CHECK_10(above > mean, scythe_invalid_arg,
            "Truncation point > mean");
        SCYTHE_CHECK_10(variance <= 0, scythe_invalid_arg,
            "Variance <= 0");
      
        double below = -1*above;
        double newmu = -1*mean;
        double z = 0;
        double x = below + .00001;
         
        for (unsigned int i=0; i<iter; ++i){
          z = runif()*std::exp(-1*std::pow((x-newmu),2)
              /(2*variance));
          x = runif()
            *( (newmu + std::sqrt(-2*variance*std::log(z))) - below) 
            + below;
        }
        if (! finite(x)) {
          SCYTHE_WARN("Mean extremely far from truncation point. "
              << "Returning truncation point");
          return above; 
        }
        
        return -1*x;
      }

      SCYTHE_RNGMETH_MATRIX(rtanorm_slice, double, 
          SCYTHE_ARGSET(mean, variance, above, iter), double mean, 
          double variance, double above, unsigned int iter = 10);

			/*! \brief Generate a normally distributed random
			 * variate, truncated below.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * normal distribution with given \a mean and \a variance,
			 * truncated below.  It uses a combination of
			 * rejection sampling (when \a mean >= \a below) and the slice
			 * sampling method of Robert and Casella (1999), pp. 288-289
			 * (when \a mean < \a below).
       *
       * \param mean The mean of the distribution.
			 * \param variance The variance of the distribution.
			 * \param below The lower truncation point of the distribution.
			 * \param iter The number of iterations to run the slice
			 * sampler.
			 * 
			 * \see rtnorm(double mean, double variance, double below, double above)
			 * \see rtnorm_combo(double mean, double variance, double below, double above)
			 * \see rtbnorm_slice(double mean, double variance, double below, unsigned int iter = 10)
			 * \see rtanorm_slice(double mean, double variance, double above, unsigned int iter = 10)
			 * \see rtanorm_combo(double mean, double variance, double above, unsigned int iter = 10)
			 * \see rnorm(double x, double mean, double sd)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rtbnorm_combo (double mean, double variance, double below, 
          unsigned int iter = 10)
      {
        SCYTHE_CHECK_10(variance <= 0, scythe_invalid_arg,
            "Variance <= 0");
        
        double s = std::sqrt(variance);
        // do rejection sampling and return value
        //if (m >= below){
        if ((mean/s - below/s ) > -0.5){
          double x = rnorm(mean, s);
          while (x < below)
            x = rnorm(mean,s);
          return x; 
        } else if ((mean/s - below/s ) > -5.0 ){
          // use the inverse cdf method
          double above =  std::numeric_limits<double>::infinity();
          double x = rtnorm(mean, variance, below, above);
          return x;
        } else {
          // do slice sampling and return value
          double z = 0;
          double x = below + .00001;
          for (unsigned int i=0; i<iter; ++i){
            z = runif() * std::exp(-1 * std::pow((x - mean), 2)
                / (2 * variance));
            x = runif() 
              * ((mean + std::sqrt(-2 * variance * std::log(z))) 
                - below) + below;
          }
          if (! finite(x)) {
            SCYTHE_WARN("Mean extremely far from truncation point. "
                << "Returning truncation point");
            return below; 
          }
          return x;
        }
      }

      SCYTHE_RNGMETH_MATRIX(rtbnorm_combo, double, 
          SCYTHE_ARGSET(mean, variance, below, iter), double mean, 
          double variance, double below, unsigned int iter = 10);

			/*! \brief Generate a normally distributed random variate,
			 * truncated above.
       *
			 * This function returns a pseudo-random variate drawn from the
			 * normal distribution with given \a mean and \a variance,
			 * truncated above.  It uses a combination of rejection sampling
			 * (when \a mean <= \a above) and the slice sampling method of
			 * Robert and Casella (1999), pp. 288-289 (when \a mean > \a
			 * above).
       *
			 * \param mean The mean of the distribution.
			 * \param variance The variance of the distribution.
			 * \param above The upper truncation point of the distribution.
			 * \param iter The number of iterations to run the slice sampler.
			 * 
			 * \see rtnorm(double mean, double variance, double below, double above)
			 * \see rtnorm_combo(double mean, double variance, double below, double above)
			 * \see rtbnorm_slice(double mean, double variance, double below, unsigned int iter = 10)
			 * \see rtanorm_slice(double mean, double variance, double above, unsigned int iter = 10)
			 * \see rtbnorm_combo(double mean, double variance, double below, unsigned int iter = 10)
			 * \see rnorm(double x, double mean, double sd)
			 *
       * \throw scythe_invalid_arg (Level 1)
       */
      double
      rtanorm_combo (double mean, double variance, double above, 
          const unsigned int iter = 10)
      {
        SCYTHE_CHECK_10(variance <= 0, scythe_invalid_arg,
            "Variance <= 0");

        double s = std::sqrt(variance);
        // do rejection sampling and return value
        if ((mean/s - above/s ) < 0.5){ 
          double x = rnorm(mean, s);
          while (x > above)
            x = rnorm(mean,s);
          return x;
        } else if ((mean/s - above/s ) < 5.0 ){
          // use the inverse cdf method
          double below =  -std::numeric_limits<double>::infinity();
          double x = rtnorm(mean, variance, below, above);
          return x;
        } else {
          // do slice sampling and return value
          double below = -1*above;
          double newmu = -1*mean;
          double z = 0;
          double x = below + .00001;
             
          for (unsigned int i=0; i<iter; ++i){
            z = runif() * std::exp(-1 * std::pow((x-newmu), 2)
                /(2 * variance));
            x = runif() 
              * ((newmu + std::sqrt(-2 * variance * std::log(z)))
                  - below) + below;
          }
          if (! finite(x)) {
            SCYTHE_WARN("Mean extremely far from truncation point. "
                << "Returning truncation point");
            return above; 
          }
          return -1*x;
        }
      }

      SCYTHE_RNGMETH_MATRIX(rtanorm_combo, double, 
          SCYTHE_ARGSET(mean, variance, above, iter), double mean, 
          double variance, double above, unsigned int iter = 10);

      /* Multivariate Distributions */
      
      /*! \brief Generate a Wishart distributed random variate Matrix.
       *
       * This function returns a pseudo-random matrix-valued variate
       * drawn from the Wishart disribution described by the scale
       * matrix \a Sigma, with \a v degrees of freedom.
       *
       * \param v The degrees of freedom of the distribution.
			 * \param Sigma The square scale matrix of the distribution.
			 *
       * \throw scythe_invalid_arg (Level 1)
       * \throw scythe_dimension_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
      Matrix<double, O, Concrete>
      rwish(unsigned int v, const Matrix<double, O, S> &Sigma)
      {
        SCYTHE_CHECK_10(! Sigma.isSquare(), scythe_dimension_error,
            "Sigma not square");
        SCYTHE_CHECK_10(v < Sigma.rows(), scythe_invalid_arg, 
            "v < Sigma.rows()");
          
        Matrix<double,O,Concrete> 
          A(Sigma.rows(), Sigma.rows());
        Matrix<double,O,Concrete> C = cholesky<O,Concrete>(Sigma);
        Matrix<double,O,Concrete> alpha;
          
        for (unsigned int i = 0; i < v; ++i) {
          alpha = C * rnorm(Sigma.rows(), 1, 0, 1);
          A = A + (alpha * (t(alpha)));
        }

        return A;
      }

      /*! \brief Generate a Dirichlet distributed random variate Matrix.
       *
       * This function returns a pseudo-random matrix-valued variate
       * drawn from the Dirichlet disribution described by the vector
       * \a alpha.
       *
       * \param alpha A vector of non-negative reals.
			 *
       * \throw scythe_invalid_arg (Level 1)
       * \throw scythe_dimension_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
      Matrix<double, O, Concrete>
      rdirich(const Matrix<double, O, S>& alpha) 
      { 
        // Check for allowable parameters
        SCYTHE_CHECK_10(std::min(alpha) <= 0, scythe_invalid_arg,
            "alpha has elements < 0");
        SCYTHE_CHECK_10(! alpha.isColVector(), scythe_dimension_error,
            "alpha not column vector");
     
        Matrix<double, O, Concrete> y(alpha.rows(), 1);
        double ysum = 0;

        // We would use std::transform here but rgamma is a function
        // and wouldn't get inlined.
        const_matrix_forward_iterator<double,O,O,S> ait;
        const_matrix_forward_iterator<double,O,O,S> alast
          = alpha.template end_f();
        typename Matrix<double,O,Concrete>::forward_iterator yit 
          = y.begin_f();
        for (ait = alpha.begin_f(); ait != alast; ++ait) {
          *yit = rgamma(*ait, 1);
          ysum += *yit;
          ++ait;
        }

        y /= ysum;

        return y;
      }

      /*! \brief Generate a multivariate normal distributed random
       * variate Matrix.
       *
       * This function returns a pseudo-random matrix-valued variate
       * drawn from the multivariate normal disribution with means \mu
       * and variance-covariance matrix \a sigma.
       *
       * \param mu A vector containing the distribution means.
       * \param sigma The distribution variance-covariance matrix.
			 *
       * \throw scythe_invalid_arg (Level 1)
       * \throw scythe_dimension_error (Level 1)
       */
      template <matrix_order PO1, matrix_style PS1,
                matrix_order PO2, matrix_style PS2>
      Matrix<double, PO1, Concrete>
      rmvnorm(const Matrix<double, PO1, PS1>& mu, 
              const Matrix<double, PO2, PS2>& sigma)
      {  
        unsigned int dim = mu.rows();
        SCYTHE_CHECK_10(! mu.isColVector(), scythe_dimension_error,
            "mu not column vector");
        SCYTHE_CHECK_10(! sigma.isSquare(), scythe_dimension_error,
            "sigma not square");
        SCYTHE_CHECK_10(sigma.rows() != dim, scythe_conformation_error,
            "mu and sigma not conformable");
        
        return(mu + cholesky(sigma) * rnorm(dim, 1, 0, 1));
      }

      /*! \brief Generate a multivariate Student t distributed random
       * variate Matrix.
       *
       * This function returns a pseudo-random matrix-valued variate
       * drawn from the multivariate Student t disribution with
       * and variance-covariance matrix \a sigma, and degrees of
       * freedom \a nu
       *
       * \param sigma The distribution variance-covariance matrix.
       * \param nu The strictly positive degrees of freedom.
			 *
       * \throw scythe_invalid_arg (Level 1)
       * \throw scythe_dimension_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
      Matrix<double, O, Concrete>
      rmvt (const Matrix<double, O, S>& sigma, double nu)
      {
        Matrix<double, O, Concrete> result;
        SCYTHE_CHECK_10(nu <= 0, scythe_invalid_arg,
            "D.O.F parameter nu <= 0");

        result = 
          rmvnorm(Matrix<double, O>(sigma.rows(), 1, true, 0), sigma);
        result /= std::sqrt(rchisq(nu) / nu);
        return result;
      }

    protected:
      /* Default (and only) constructor */
      /*! \brief Default constructor
       *
       * Instantiate a random number generator
       */
      rng() 
        : rnorm_count_ (1) // Initialize the normal counter
      {}

      /* For Barton and Nackman trick. */
      RNGTYPE& as_derived()
      {
        return static_cast<RNGTYPE&>(*this);
      }


      /* Generate Standard Normal variates */

      /* These instance variables were static in the old
       * implementation.  Making them instance variables provides
       * thread safety, as long as two threads don't access the same
       * rng at the same time w/out precautions.  Fixes possible
       * previous issues with lecuyer.  See the similar approach in
       * rgamma1 below.
       */
      int rnorm_count_;
      double x2_;

      double
      rnorm1 ()
      {
        double nu1, nu2, rsquared, sqrt_term;
        if (rnorm_count_ == 1){ // odd numbered passses
          do {
            nu1 = -1 +2*runif();
            nu2 = -1 +2*runif();
            rsquared = ::pow(nu1,2) + ::pow(nu2,2);
          } while (rsquared >= 1 || rsquared == 0.0);
          sqrt_term = std::sqrt(-2*std::log(rsquared)/rsquared);
          x2_ = nu2*sqrt_term;
          rnorm_count_ = 2;
          return nu1*sqrt_term;
        } else { // even numbered passes
          rnorm_count_ = 1;
          return x2_;
        } 
      }

      /* Generate standard gamma variates */
      double accept_;

      double
      rgamma1 (double alpha)
      {
        int test;
        double u, v, w, x, y, z, b, c;

        // Check for allowable parameters
        SCYTHE_CHECK_10(alpha <= 1, scythe_invalid_arg, "alpha <= 1");

        // Implement Best's (1978) simulator
        b = alpha - 1;
        c = 3 * alpha - 0.75;
        test = 0;
        while (test == 0) {
          u = runif ();
          v = runif ();

          w = u * (1 - u);
          y = std::sqrt (c / w) * (u - .5);
          x = b + y;

          if (x > 0) {
            z = 64 * std::pow (v, 2) * std::pow (w, 3);
            if (z <= (1 - (2 * std::pow (y, 2) / x))) {
              test = 1;
              accept_ = x;
            } else if ((2 * (b * std::log (x / b) - y)) >= ::log (z)) {
              test = 1;
              accept_ = x;
            } else {
              test = 0;
            }
          }
        }
        
        return (accept_);
      }

  };

  
} // end namespace scythe    
#endif /* RNG_H */
