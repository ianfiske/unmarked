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
 *  scythestat/distributions.h
 *
 */

/*!  \file distributions.h
 *
 * \brief Definitions for probability density functions
 * (PDFs), cumulative distribution functions (CDFs), and some common
 * functions (gamma, beta, etc).
 *
 * This file provides functions that evaluate the PDFs and CDFs of a
 * number of probability distributions.  In addition, it includes
 * definitions for another of related functions, such as the gamma
 * and beta functions.
 *
 * The various distribution functions in this file operate on both
 * scalar quantiles and matrices of quantiles and the
 * definitions of both forms of these functions appear below.  We
 * provide explicit documentation only for the scalar versions of the
 * these functions and describe the Matrix versions in the scalar
 * calls' documents.  Much like the operators in matrix.h, we
 * implement these overloaded versions of the distribution functions
 * in terms of both generalized and default templates to allow for
 * explicit control over the template type of the returned Matrix.
 *
 * \note Doxygen does not correctly expand the macro definitions we use
 * to generate the Matrix versions of the various distribution
 * functions.  Therefore, it incorrectly substitutes the macro
 * variable
 * __VA_ARGS__ for the actual parameter values in the parameter lists
 * of each of these functions.  For example, the definitions of the
 * Matrix versions of pbeta are listed as
 * \code
 * template<matrix_order RO, matrix_style RS, matrix_order PO, matrix_style PS>
 * Matrix<double, RO, RS> scythe::pbeta (const Matrix<double, PO, PS> &X, __VA_ARGS__)
 *
 * template<matrix_order O, matrix_style S>
 * Matrix<double, O, Concrete> scythe::pbeta (const Matrix<double, O, S> &X, __VA_ARGS__)
 * \endcode
 * when they should be
 * \code
 * template<matrix_order RO, matrix_style RS, matrix_order PO, matrix_style PS>
 * Matrix<double, RO, RS> scythe::pbeta (const Matrix<double, PO, PS> &X, double a, double b)
 *
 * template<matrix_order O, matrix_style S>
 * Matrix<double, O, Concrete> scythe::pbeta (const Matrix<double, O, S> &X, double a, double b)
 * \endcode
 *
 * \par
 * Furthermore, Doxygen erroneously lists a number of variables at the
 * end of this document that are not, in fact, declared in
 * distributions.h.  Again, this error is a result of Doxygen's macro
 * parsing capabilities.
 * 
 */

/* TODO: Figure out how to get doxygen to stop listing erroneous
 * variables at the end of the doc for this file.  They stem from it
 * misreading the nested macro calls used to generate matrix procs.
 */

/* TODO: Full R-style versions of these function including arbitrary
 * recycling of matrix arguments.  This is going to have to wait for
 * variadic templates to be doable without a complete mess.  There is
 * currently a variadic patch available for g++, perhaps we can add a
 * SCYTHE_VARIADIC flag and include these as option until they become
 * part of the standard in 2009.  Something to come back to after 1.0.
 */

#ifndef SCYTHE_DISTRIBUTIONS_H
#define SCYTHE_DISTRIBUTIONS_H

#include <iostream>
#include <cmath>
#include <cfloat>
#include <climits>
#include <algorithm>
#include <limits>

#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#include "ide.h"
#include "error.h"
#else
#include "scythestat/matrix.h"
#include "scythestat/ide.h"
#include "scythestat/error.h"
#endif

/* Fill in some defs from R that aren't in math.h */
#ifndef M_PI
#define M_PI 3.141592653589793238462643383280
#endif
#define M_LN_SQRT_2PI 0.918938533204672741780329736406
#define M_LN_SQRT_PId2  0.225791352644727432363097614947
#define M_1_SQRT_2PI  0.39894228040143267793994605993
#define M_2PI   6.28318530717958647692528676655
#define M_SQRT_32 5.656854249492380195206754896838

#ifndef HAVE_TRUNC
/*! @cond */
inline double trunc(double x) throw ()
{
    if (x >= 0) 
      return std::floor(x);
    else
      return std::ceil(x);
}
/*! @endcond */
#endif

/* Many random number generators, pdfs, cdfs, and functions (gamma,
 * etc) in this file are based on code from the R Project, version
 * 1.6.0-1.7.1.  This code is available under the terms of the GNU
 * GPL.  Original copyright:
 *
 * Copyright (C) 1998      Ross Ihaka
 * Copyright (C) 2000-2002 The R Development Core Team
 * Copyright (C) 2003      The R Foundation
 */

namespace scythe {
  
  /*! @cond */
  
  /* Forward declarations */
  double gammafn (double);
  double lngammafn (double);
  double lnbetafn (double, double);
  double pgamma (double, double, double);
  double dgamma(double, double, double);
  double pnorm (double, double, double);

  /*! @endcond */

  /********************
   * Helper Functions *
   ********************/
  namespace {

    /* Evaluate a Chebysheve series at a given point */
    double
    chebyshev_eval (double x, const double *a, int n)
    {
      SCYTHE_CHECK_10(n < 1 || n > 1000, scythe_invalid_arg,
          "n not on [1, 1000]");
  
      SCYTHE_CHECK_10(x < -1.1 || x > 1.1, scythe_invalid_arg,
          "x not on [-1.1, 1.1]");
      
      double b0, b1, b2;
      b0 = b1 = b2 = 0;
  
      double twox = x * 2;
  
      for (int i = 1; i <= n; ++i) {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
      }
  
      return (b0 - b2) * 0.5;
    }

    /* Computes the log gamma correction factor for x >= 10 */
    double
    lngammacor(double x)
    {
      const double algmcs[15] = {
        +.1666389480451863247205729650822e+0,
        -.1384948176067563840732986059135e-4,
        +.9810825646924729426157171547487e-8,
        -.1809129475572494194263306266719e-10,
        +.6221098041892605227126015543416e-13,
      };
    
      SCYTHE_CHECK_10(x < 10, scythe_invalid_arg,
          "This function requires x >= 10");  
      SCYTHE_CHECK_10(x >= 3.745194030963158e306, scythe_range_error,
          "Underflow");
      
      if (x < 94906265.62425156) {
        double tmp = 10 / x;
        return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, 5) / x;
      }
      
      return 1 / (x * 12);
    }

    /* Evaluates the "deviance part" */
    double
    bd0(double x, double np)
    {
      
      if(std::fabs(x - np) < 0.1 * (x + np)) {
        double v = (x - np) / (x + np);
        double s = (x - np) * v;
        double ej = 2 * x * v;
        v = v * v;
        for (int j = 1; ; j++) {
          ej *= v;
          double s1 = s + ej / ((j << 1) + 1);
          if (s1 == s)
            return s1;
          s = s1;
        }
      }
      
      return x * std::log(x / np) + np - x;
    }
  
    /* Computes the log of the error term in Stirling's formula */
    double
    stirlerr(double n)
    {
#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
      
      /* error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0 */
      const double sferr_halves[31] = {
        0.0, /* n=0 - wrong, place holder only */
        0.1534264097200273452913848,  /* 0.5 */
        0.0810614667953272582196702,  /* 1.0 */
        0.0548141210519176538961390,  /* 1.5 */
        0.0413406959554092940938221,  /* 2.0 */
        0.03316287351993628748511048, /* 2.5 */
        0.02767792568499833914878929, /* 3.0 */
        0.02374616365629749597132920, /* 3.5 */
        0.02079067210376509311152277, /* 4.0 */
        0.01848845053267318523077934, /* 4.5 */
        0.01664469118982119216319487, /* 5.0 */
        0.01513497322191737887351255, /* 5.5 */
        0.01387612882307074799874573, /* 6.0 */
        0.01281046524292022692424986, /* 6.5 */
        0.01189670994589177009505572, /* 7.0 */
        0.01110455975820691732662991, /* 7.5 */
        0.010411265261972096497478567, /* 8.0 */
        0.009799416126158803298389475, /* 8.5 */
        0.009255462182712732917728637, /* 9.0 */
        0.008768700134139385462952823, /* 9.5 */
        0.008330563433362871256469318, /* 10.0 */
        0.007934114564314020547248100, /* 10.5 */
        0.007573675487951840794972024, /* 11.0 */
        0.007244554301320383179543912, /* 11.5 */
        0.006942840107209529865664152, /* 12.0 */
        0.006665247032707682442354394, /* 12.5 */
        0.006408994188004207068439631, /* 13.0 */
        0.006171712263039457647532867, /* 13.5 */
        0.005951370112758847735624416, /* 14.0 */
        0.005746216513010115682023589, /* 14.5 */
        0.005554733551962801371038690  /* 15.0 */
      };
      double nn;
      
      if (n <= 15.0) {
        nn = n + n;
        if (nn == (int)nn)
          return(sferr_halves[(int)nn]);
        return (scythe::lngammafn(n + 1.) - (n + 0.5) * std::log(n) + n -
            std::log(std::sqrt(2 * M_PI)));
      }
      
      nn = n*n;
      if (n > 500)
        return((S0 - S1 / nn) / n);
      if (n > 80)
        return((S0 - (S1 - S2 / nn) / nn) / n);
      if (n > 35)
        return((S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n);
      /* 15 < n <= 35 : */
      return((S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
#undef S1
#undef S2
#undef S3
#undef S4
    }


    /* Helper for dpois and dgamma */
    double
    dpois_raw (double x, double lambda)
    {
      if (lambda == 0)
        return ( (x == 0) ? 1.0 : 0.0);

      if (x == 0)
        return std::exp(-lambda);

      if (x < 0)
        return 0.0;

      return std::exp(-stirlerr(x) - bd0(x, lambda))
        / std::sqrt(2 * M_PI * x);
    }

  
    /* helper for pbeta */
    double
    pbeta_raw(double x, double pin, double qin)
    {
      double ans, c, finsum, p, ps, p1, q, term, xb, xi, y;
      int n, i, ib, swap_tail;
      
      const double eps = .5 * DBL_EPSILON;
      const double sml = DBL_MIN;
      const double lneps = std::log(eps);
      const double lnsml = std::log(eps);
      
      if (pin / (pin + qin) < x) {
        swap_tail = 1;
        y = 1 - x;
        p = qin;
        q = pin;
      } else {
        swap_tail=0;
        y = x;
        p = pin;
        q = qin;
      }
      
      if ((p + q) * y / (p + 1) < eps) {
        ans = 0;
        xb = p * std::log(std::max(y,sml)) - std::log(p) - lnbetafn(p,q);
        if (xb > lnsml && y != 0)
          ans = std::exp(xb);
        if (swap_tail)
          ans = 1-ans;
      } else {
        ps = q - std::floor(q);
        if (ps == 0)
          ps = 1;
        xb = p * std::log(y) - lnbetafn(ps, p) - std::log(p);
        ans = 0;
        if (xb >= lnsml) {
          ans = std::exp(xb);
          term = ans * p;
          if (ps != 1) {
            n = (int)std::max(lneps/std::log(y), 4.0);
            for(i = 1; i <= n; i++){
              xi = i;
              term *= (xi-ps)*y/xi;
              ans += term/(p+xi);
            }
          }
        }
        if (q > 1) {
          xb = p * std::log(y) + q * std::log(1 - y)
            - lnbetafn(p, q) - std::log(q);
          ib = (int) std::max(xb / lnsml, 0.0);
          term = std::exp(xb - ib * lnsml);
          c = 1 / (1 - y);
          p1 = q * c / (p + q - 1);
              
          finsum = 0;
          n = (int) q;
          if(q == n)
            n--;
          for (i = 1; i <= n; i++) {
            if(p1 <= 1 && term / eps <= finsum)
              break;
            xi = i;
            term = (q -xi + 1) * c * term / (p + q - xi);
            if (term > 1) {
              ib--;
              term *= sml;
            }
            if (ib == 0)
              finsum += term;
          }
          ans += finsum;
        }
        
        if(swap_tail)
          ans = 1-ans;
        ans = std::max(std::min(ans,1.),0.);
      }
      return ans;
    }
  
   /* Helper for dbinom */
    double
    dbinom_raw (double x, double n, double p, double q)
    { 
      double f, lc;

      if (p == 0)
        return((x == 0) ? 1.0 : 0.0);
      if (q == 0)
        return((x == n) ? 1.0 : 0.0);

      if (x == 0) { 
        if(n == 0)
          return 1.0;
        
        lc = (p < 0.1) ? -bd0(n, n * q) - n * p : n * std::log(q);
        return(std::exp(lc));
      }
      if (x == n) { 
        lc = (q < 0.1) ? -bd0(n,n * p) - n * q : n * std::log(p);
        return(std::exp(lc));
      }

      if (x < 0 || x > n)
        return 0.0;

      lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x,n*p) -
        bd0(n - x, n * q);
      
      f = (M_2PI * x * (n-x)) / n;

      return (std::exp(lc) / std::sqrt(f));
    }

    /* The normal probability density function implementation. */

#define SIXTEN 16
#define do_del(X)              \
    xsq = trunc(X * SIXTEN) / SIXTEN;        \
    del = (X - xsq) * (X + xsq);          \
    if(log_p) {              \
        *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + std::log(temp);  \
        if((lower && x > 0.) || (upper && x <= 0.))      \
        *ccum = ::log1p(-std::exp(-xsq * xsq * 0.5) *     \
          std::exp(-del * 0.5) * temp);    \
    }                \
    else {                \
        *cum = std::exp(-xsq * xsq * 0.5) * std::exp(-del * 0.5) * temp;  \
        *ccum = 1.0 - *cum;            \
    }

#define swap_tail            \
    if (x > 0.) {/* swap  ccum <--> cum */      \
        temp = *cum; if(lower) *cum = *ccum; *ccum = temp;  \
    }

    void
    pnorm_both(double x, double *cum, double *ccum, int i_tail,
                bool log_p)
    {
      const double a[5] = {
        2.2352520354606839287,
        161.02823106855587881,
        1067.6894854603709582,
        18154.981253343561249,
        0.065682337918207449113
      };
      const double b[4] = {
        47.20258190468824187,
        976.09855173777669322,
        10260.932208618978205,
        45507.789335026729956
      };
      const double c[9] = {
        0.39894151208813466764,
        8.8831497943883759412,
        93.506656132177855979,
        597.27027639480026226,
        2494.5375852903726711,
        6848.1904505362823326,
        11602.651437647350124,
        9842.7148383839780218,
        1.0765576773720192317e-8
      };
      const double d[8] = {
        22.266688044328115691,
        235.38790178262499861,
        1519.377599407554805,
        6485.558298266760755,
        18615.571640885098091,
        34900.952721145977266,
        38912.003286093271411,
        19685.429676859990727
      };
      const double p[6] = {
        0.21589853405795699,
        0.1274011611602473639,
        0.022235277870649807,
        0.001421619193227893466,
        2.9112874951168792e-5,
        0.02307344176494017303
      };
      const double q[5] = {
        1.28426009614491121,
        0.468238212480865118,
        0.0659881378689285515,
        0.00378239633202758244,
        7.29751555083966205e-5
      };
      
      double xden, xnum, temp, del, eps, xsq, y;
      int i, lower, upper;

      /* Consider changing these : */
      eps = DBL_EPSILON * 0.5;

      /* i_tail in {0,1,2} =^= {lower, upper, both} */
      lower = i_tail != 1;
      upper = i_tail != 0;

      y = std::fabs(x);
      if (y <= 0.67448975) {
        /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
        if (y > eps) {
          xsq = x * x;
          xnum = a[4] * xsq;
          xden = xsq;
          for (i = 0; i < 3; ++i) {
            xnum = (xnum + a[i]) * xsq;
            xden = (xden + b[i]) * xsq;
          }
        } else xnum = xden = 0.0;
        
        temp = x * (xnum + a[3]) / (xden + b[3]);
        if(lower)  *cum = 0.5 + temp;
        if(upper) *ccum = 0.5 - temp;
        if(log_p) {
          if(lower)  *cum = std::log(*cum);
          if(upper) *ccum = std::log(*ccum);
        }
      } else if (y <= M_SQRT_32) {
        /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) 
         * ~= 5.657 */

        xnum = c[8] * y;
        xden = y;
        for (i = 0; i < 7; ++i) {
          xnum = (xnum + c[i]) * y;
          xden = (xden + d[i]) * y;
        }
        temp = (xnum + c[7]) / (xden + d[7]);
        do_del(y);
        swap_tail;
      } else if (log_p
                || (lower && -37.5193 < x && x < 8.2924)
                || (upper && -8.2929 < x && x < 37.5193)
          ) {
        /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
        xsq = 1.0 / (x * x);
        xnum = p[5] * xsq;
        xden = xsq;
        for (i = 0; i < 4; ++i) {
          xnum = (xnum + p[i]) * xsq;
          xden = (xden + q[i]) * xsq;
        }
        temp = xsq * (xnum + p[4]) / (xden + q[4]);
        temp = (M_1_SQRT_2PI - temp) / y;
        do_del(x);
        swap_tail;
      } else {
        if (x > 0) {
          *cum = 1.;
          *ccum = 0.;
        } else {
          *cum = 0.;
          *ccum = 1.;
        }
        SCYTHE_THROW_10(scythe_convergence_error, "Did not converge");
      }

      return;
    }
#undef SIXTEN
#undef do_del
#undef swap_tail

    /* The standard normal distribution function */
    double
    pnorm1 (double x, bool lower_tail, bool log_p)
    {
      SCYTHE_CHECK_10(! finite(x), scythe_invalid_arg,
          "Quantile x is inifinte (+/-Inf) or NaN");

      double p, cp;
      pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

      return (lower_tail ? p : cp);
    }
  } // anonymous namespace

  /*************
   * Functions *
   *************/
  
  /* The gamma function */

  /*! \brief The gamma function.
   *
   * Computes the gamma function, evaluated at \a x.
   *
   * \param x The value to compute gamma at.
   *
   * \see lngammafn(double x)
   * \see pgamma(double x, double shape, double scale)
   * \see dgamma(double x, double shape, double scale)
   * \see rng::rgamma(double shape, double scale)
   *
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double 
  gammafn (double x)
  {
    const double gamcs[22] = {
      +.8571195590989331421920062399942e-2,
      +.4415381324841006757191315771652e-2,
      +.5685043681599363378632664588789e-1,
      -.4219835396418560501012500186624e-2,
      +.1326808181212460220584006796352e-2,
      -.1893024529798880432523947023886e-3,
      +.3606925327441245256578082217225e-4,
      -.6056761904460864218485548290365e-5,
      +.1055829546302283344731823509093e-5,
      -.1811967365542384048291855891166e-6,
      +.3117724964715322277790254593169e-7,
      -.5354219639019687140874081024347e-8,
      +.9193275519859588946887786825940e-9,
      -.1577941280288339761767423273953e-9,
      +.2707980622934954543266540433089e-10,
      -.4646818653825730144081661058933e-11,
      +.7973350192007419656460767175359e-12,
      -.1368078209830916025799499172309e-12,
      +.2347319486563800657233471771688e-13,
      -.4027432614949066932766570534699e-14,
      +.6910051747372100912138336975257e-15,
      -.1185584500221992907052387126192e-15,
    };


    double y = std::fabs(x);

    if (y <= 10) {

      /* Compute gamma(x) for -10 <= x <= 10
       * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
       * first of all. */

      int n = (int) x;
      if (x < 0)
        --n;
      
      y = x - n;/* n = floor(x)  ==>  y in [ 0, 1 ) */
      --n;
      double value = chebyshev_eval(y * 2 - 1, gamcs, 22) + .9375;
      
      if (n == 0)
        return value;/* x = 1.dddd = 1+y */

      if (n < 0) {
        /* compute gamma(x) for -10 <= x < 1 */

        /* If the argument is exactly zero or a negative integer */
        /* then return NaN. */
        SCYTHE_CHECK_10(x == 0 || (x < 0 && x == n + 2),
            scythe_range_error, "x is 0 or a negative integer");

        /* The answer is less than half precision */
        /* because x too near a negative integer. */
        SCYTHE_CHECK_10(x < -0.5 && 
            std::fabs(x - (int)(x - 0.5) / x) < 67108864.0,
            scythe_precision_error,
            "Answer < 1/2 precision because x is too near" <<
            " a negative integer");

        /* The argument is so close to 0 that the result
         * * would overflow. */
        SCYTHE_CHECK_10(y < 2.2474362225598545e-308, scythe_range_error,
            "x too close to 0");

        n = -n;

        for (int i = 0; i < n; i++)
          value /= (x + i);
        
        return value;
      } else {
        /* gamma(x) for 2 <= x <= 10 */

        for (int i = 1; i <= n; i++) {
          value *= (y + i);
        }
        return value;
      }
    } else {
      /* gamma(x) for   y = |x| > 10. */

      /* Overflow */
      SCYTHE_CHECK_10(x > 171.61447887182298, 
          scythe_range_error,"Overflow");

      /* Underflow */
      SCYTHE_CHECK_10(x < -170.5674972726612,
          scythe_range_error, "Underflow");

      double value = std::exp((y - 0.5) * std::log(y) - y 
          + M_LN_SQRT_2PI + lngammacor(y));

      if (x > 0)
        return value;

      SCYTHE_CHECK_10(std::fabs((x - (int)(x - 0.5))/x) < 67108864.0,
          scythe_precision_error, 
          "Answer < 1/2 precision because x is " <<
            "too near a negative integer");

      double sinpiy = std::sin(M_PI * y);

      /* Negative integer arg - overflow */
      SCYTHE_CHECK_10(sinpiy == 0, scythe_range_error, "Overflow");

      return -M_PI / (y * sinpiy * value);
    }
  }

  /* The natural log of the absolute value of the gamma function */
  /*! \brief The natural log of the absolute value of the gamma 
   * function.
   *
   * Computes the natural log of the absolute value of the gamma 
   * function, evaluated at \a x.
   *
   * \param x The value to compute log(abs(gamma())) at.
   *
   * \see gammafn(double x)
   * \see pgamma(double x, double shape, double scale)
   * \see dgamma(double x, double shape, double scale)
   * \see rng::rgamma(double shape, double scale)
   *
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  lngammafn(double x)
  {
    SCYTHE_CHECK_10(x <= 0 && x == (int) x, scythe_range_error,
        "x is 0 or a negative integer");

    double y = std::fabs(x);

    if (y <= 10)
      return std::log(std::fabs(gammafn(x)));

    SCYTHE_CHECK_10(y > 2.5327372760800758e+305, scythe_range_error,
        "Overflow");

    if (x > 0) /* i.e. y = x > 10 */
      return M_LN_SQRT_2PI + (x - 0.5) * std::log(x) - x
        + lngammacor(x);
    
    /* else: x < -10; y = -x */
    double sinpiy = std::fabs(std::sin(M_PI * y));

    if (sinpiy == 0) /* Negative integer argument */
      throw scythe_exception("UNEXPECTED ERROR",
           __FILE__, __func__, __LINE__,
           "ERROR:  Should never happen!");

    double ans = M_LN_SQRT_PId2 + (x - 0.5) * std::log(y) - x - std::log(sinpiy)
      - lngammacor(y);

    SCYTHE_CHECK_10(std::fabs((x - (int)(x - 0.5)) * ans / x) 
        < 1.490116119384765696e-8, scythe_precision_error, 
        "Answer < 1/2 precision because x is " 
        << "too near a negative integer");
    
    return ans;
  }

  /* The beta function */
  /*! \brief The beta function.
   *
   * Computes beta function, evaluated at (\a a, \a b).
   *
   * \param a The first parameter.
   * \param b The second parameter.
   *
   * \see lnbetafn(double a, double b)
   * \see pbeta(double x, double a, double b)
   * \see dbeta(double x, double a, double b)
   * \see rng::rbeta(double a, double b)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  betafn(double a, double b)
  {
    SCYTHE_CHECK_10(a <= 0 || b <= 0, scythe_invalid_arg, "a or b < 0");

    if (a + b < 171.61447887182298) /* ~= 171.61 for IEEE */
      return gammafn(a) * gammafn(b) / gammafn(a+b);

    double val = lnbetafn(a, b);
    SCYTHE_CHECK_10(val < -708.39641853226412, scythe_range_error,
        "Underflow");
    
    return std::exp(val);
  }

  /* The natural log of the beta function */
  /*! \brief The natural log of the beta function.
   *
   * Computes the natural log of the beta function, 
   * evaluated at (\a a, \a b).
   *
   * \param a The first parameter.
   * \param b The second parameter.
   *
   * \see betafn(double a, double b)
   * \see pbeta(double x, double a, double b)
   * \see dbeta(double x, double a, double b)
   * \see rng::rbeta(double a, double b)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  lnbetafn (double a, double b)
  {
    double p, q;

    p = q = a;
    if(b < p) p = b;/* := min(a,b) */
    if(b > q) q = b;/* := max(a,b) */

    SCYTHE_CHECK_10(p <= 0 || q <= 0,scythe_invalid_arg, "a or b <= 0");

    if (p >= 10) {
      /* p and q are big. */
      double corr = lngammacor(p) + lngammacor(q) - lngammacor(p + q);
      return std::log(q) * -0.5 + M_LN_SQRT_2PI + corr
        + (p - 0.5) * std::log(p / (p + q)) + q * std::log(1 + (-p / (p + q)));
    } else if (q >= 10) {
      /* p is small, but q is big. */
      double corr = lngammacor(q) - lngammacor(p + q);
      return lngammafn(p) + corr + p - p * std::log(p + q)
        + (q - 0.5) * std::log(1 + (-p / (p + q)));
    }
    
    /* p and q are small: p <= q > 10. */
    return std::log(gammafn(p) * (gammafn(q) / gammafn(p + q)));
  }

  /* Compute the factorial of a non-negative integer */
  /*! \brief The factorial function.
   *
   * Computes the factorial of \a n.
   *
   * \param n The non-negative integer value to compute the factorial of.
   *
   * \see lnfactorial(unsigned int n)
   *
   */
  inline int
  factorial (unsigned int n)
  {
    if (n == 0)
      return 1;

    return n * factorial(n - 1);
  }

  /* Compute the natural log of the factorial of a non-negative
   * integer
   */
  /*! \brief The log of the factorial function.
   *
   * Computes the natural log of the factorial of \a n.
   *
   * \param n The non-negative integer value to compute the natural log of the factorial of.
   *
   * \see factorial(unsigned int n)
   *
   */
  inline double
  lnfactorial (unsigned int n)
  {
    double x = n+1;
    double cof[6] = {
      76.18009172947146, -86.50532032941677,
      24.01409824083091, -1.231739572450155,
      0.1208650973866179e-2, -0.5395239384953e-5
    };
    double y = x;
    double tmp = x + 5.5 - (x + 0.5) * std::log(x + 5.5);
    double ser = 1.000000000190015;
    for (int j = 0; j <= 5; j++) {
      ser += (cof[j] / ++y);
    }
    return(std::log(2.5066282746310005 * ser / x) - tmp);
  }

  /*********************************
   * Fully Specified Distributions * 
   *********************************/

  /* These macros provide a nice shorthand for the matrix versions of
   * the pdf and cdf functions.
   */
 
#define SCYTHE_ARGSET(...) __VA_ARGS__

#define SCYTHE_DISTFUN_MATRIX(NAME, XTYPE, ARGNAMES, ...)             \
  template <matrix_order RO, matrix_style RS,                         \
            matrix_order PO, matrix_style PS>                         \
  Matrix<double, RO, RS>                                              \
  NAME (const Matrix<XTYPE, PO, PS>& X, __VA_ARGS__)                  \
  {                                                                   \
    Matrix<double, RO, Concrete> ret(X.rows(), X.cols(), false);      \
    const_matrix_forward_iterator<XTYPE,RO,PO,PS> xit;                \
    const_matrix_forward_iterator<XTYPE,RO,PO,PS> xlast               \
      = X.template end_f<RO>();                                       \
    typename Matrix<double,RO,Concrete>::forward_iterator rit         \
      = ret.begin_f();                                                \
    for (xit = X.template begin_f<RO>(); xit != xlast; ++xit) {       \
      *rit = NAME (*xit, ARGNAMES);                                   \
      ++rit;                                                          \
    }                                                                 \
    SCYTHE_VIEW_RETURN(double, RO, RS, ret)                           \
  }                                                                   \
                                                                      \
  template <matrix_order O, matrix_style S>                           \
  Matrix<double, O, Concrete>                                         \
  NAME (const Matrix<XTYPE, O, S>& X, __VA_ARGS__)                    \
  {                                                                   \
    return NAME <O, Concrete, O, S> (X, ARGNAMES);                    \
  }

  /**** The Beta Distribution ****/

  /* CDFs */

  /*! \brief The beta distribution function.
   *
   * Computes the value of the beta cumulative distribution function
   * with shape parameters \a a and \a b at the desired quantile,
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile, between 0 and 1.
   * \param a The first non-negative beta shape parameter.
   * \param b The second non-negative beta shape parameter.
   *
   * \see dbeta(double x, double a, double b)
   * \see rng::rbeta(double a, double b)
   * \see betafn(double a, double b)
   * \see lnbetafn(double a, double b)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  pbeta(double x, double a, double b)
  {
    SCYTHE_CHECK_10(a <= 0 || b <= 0,scythe_invalid_arg, "a or b <= 0");
    
    if (x <= 0)
      return 0.;
    if (x >= 1)
      return 1.;
    
    return pbeta_raw(x,a,b);
  }

  SCYTHE_DISTFUN_MATRIX(pbeta, double, SCYTHE_ARGSET(a, b), double a, double b)

  /* PDFs */
  /*! \brief The beta density function.
   *
   * Computes the value of the beta probability density function
   * with shape parameters \a a and \a b at the desired quantile,
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile, between 0 and 1.
   * \param a The first non-negative beta shape parameter.
   * \param b The second non-negative beta shape parameter.
   *
   * \see pbeta(double x, double a, double b)
   * \see rng::rbeta(double a, double b)
   * \see betafn(double a, double b)
   * \see lnbetafn(double a, double b)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  dbeta(double x, double a, double b)
  {
    SCYTHE_CHECK_10((x < 0.0) || (x > 1.0), scythe_invalid_arg,
        "x not in [0,1]");
    SCYTHE_CHECK_10(a < 0.0, scythe_invalid_arg, "a < 0");
    SCYTHE_CHECK_10(b < 0.0, scythe_invalid_arg, "b < 0");

    return (std::pow(x, (a-1.0)) * std::pow((1.0-x), (b-1.0)) )
      / betafn(a,b);
  }

  SCYTHE_DISTFUN_MATRIX(dbeta, double, SCYTHE_ARGSET(a, b), double a, double b)

  /* Returns the natural log of the ordinate of the Beta density
   * evaluated at x with Shape1 a, and Shape2 b
   */

  /*! \brief The natural log of the ordinate of the beta density
   * function.
   *
   * Computes the value of the natural log of the ordinate of the beta
   * probability density function
   * with shape parameters \a a and \a b at the desired quantile,
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile, between 0 and 1.
   * \param a The first non-negative beta shape parameter.
   * \param b The second non-negative beta shape parameter.
   *
   * \see dbeta(double x, double a, double b)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  lndbeta1(double x, double a, double b)
  { 
    SCYTHE_CHECK_10((x < 0.0) || (x > 1.0), scythe_invalid_arg,
        "x not in [0,1]");
    SCYTHE_CHECK_10(a < 0.0, scythe_invalid_arg, "a < 0");
    SCYTHE_CHECK_10(b < 0.0, scythe_invalid_arg, "b < 0");

    return (a-1.0) * std::log(x) + (b-1) * std::log(1.0-x)
      - lnbetafn(a,b);
  }

  SCYTHE_DISTFUN_MATRIX(lndbeta1, double, SCYTHE_ARGSET(a, b), double a, double b)


  /**** The Binomial Distribution ****/

  /* CDFs */

  /*! \brief The binomial distribution function.
   *
   * Computes the value of the binomial cumulative distribution function
   * with \a n trials and \a p probability of success on each trial,
   * at the desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param n The number of trials.
   * \param p The probability of success on each trial.
   *
   * \see dbinom(double x, unsigned int n, double p)
   * \see rng::rbinom(unsigned int n, double p)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  pbinom(double x, unsigned int n, double p)
  {
      
    SCYTHE_CHECK_10(p < 0 || p > 1, scythe_invalid_arg, "p not in [0,1]");
    double X = std::floor(x);
      
    if (X < 0.0)
      return 0;
    
    if (n <= X)
      return 1;
      
    return pbeta(1 - p, n - X, X + 1);
  }

  SCYTHE_DISTFUN_MATRIX(pbinom, double, SCYTHE_ARGSET(n, p), unsigned int n, double p)

  /* PDFs */
  /*! \brief The binomial density function.
   *
   * Computes the value of the binomial probability density function
   * with \a n trials and \a p probability of success on each trial,
   * at the desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param n The number of trials.
   * \param p The probability of success on each trial.
   *
   * \see pbinom(double x, unsigned int n, double p)
   * \see rng::rbinom(unsigned int n, double p)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  dbinom(double x, unsigned int n, double p)
  {
    SCYTHE_CHECK_10(p < 0 || p > 1, scythe_invalid_arg, "p not in [0, 1]");
    double X = std::floor(x + 0.5);
    return dbinom_raw(X, n, p, 1 - p);
  }

  SCYTHE_DISTFUN_MATRIX(dbinom, double, SCYTHE_ARGSET(n, p), unsigned int n, double p)

  /**** The Chi Squared Distribution ****/
  
  /* CDFs */
  /*! \brief The \f$\chi^2\f$ distribution function.
   *
   * Computes the value of the \f$\chi^2\f$ cumulative distribution
   * function with \a df degrees of freedom, at the desired quantile
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param df The degrees of freedom.

   * \see dchisq(double x, double df)
   * \see rng::rchisq(double df)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   * \throw scythe_convergence_error (Level 1)
   *
   */
  inline double
  pchisq(double x, double df)
  {
    return pgamma(x, df/2.0, 2.0);
  }

  SCYTHE_DISTFUN_MATRIX(pchisq, double, df, double df)

  /* PDFs */
  /*! \brief The \f$\chi^2\f$ density function.
   *
   * Computes the value of the \f$\chi^2\f$ probability density
   * function with \a df degrees of freedom, at the desired quantile
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param df The degrees of freedom.

   * \see pchisq(double x, double df)
   * \see rng::rchisq(double df)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   * \throw scythe_convergence_error (Level 1)
   *
   */
  inline double
  dchisq(double x, double df)
  {
    return dgamma(x, df / 2.0, 2.0);
  }

  SCYTHE_DISTFUN_MATRIX(dchisq, double, df, double df)

  /**** The Exponential Distribution ****/

  /* CDFs */
  /*! \brief The exponential distribution function.
   *
   * Computes the value of the exponential cumulative distribution
   * function with given \a scale, at the desired quantile
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param scale The positive scale of the function.
   *
   * \see dexp(double x, double scale)
   * \see rng::rexp(double scale)
   *
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  pexp(double x, double scale)
  {
    SCYTHE_CHECK_10(scale <= 0, scythe_invalid_arg, "scale <= 0");
    
    if (x <= 0)
      return 0;
    
    return (1 - std::exp(-x*scale));
  }

  SCYTHE_DISTFUN_MATRIX(pexp, double, scale, double scale)

  /* PDFs */
  /*! \brief The exponential density function.
   *
   * Computes the value of the exponential probability density
   * function with given \a scale, at the desired quantile
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param scale The positive scale of the function.
   *
   * \see pexp(double x, double scale)
   * \see rng::rexp(double scale)
   *
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  dexp(double x, double scale)
  {
    SCYTHE_CHECK_10(scale <= 0, scythe_invalid_arg, "scale <= 0");
    
    if (x < 0)
      return 0;
      
    return std::exp(-x * scale) * scale;
  }

  SCYTHE_DISTFUN_MATRIX(dexp, double, scale, double scale)

  /**** The f Distribution ****/

  /* CDFs */
  /*! \brief The F distribution function.
   *
   * Computes the value of the F cumulative distribution function with
   * \a df1 and \a df2 degrees of freedom, at the desired quantile \a
   * x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param df1 The non-negative degrees of freedom for the
   * \f$\chi^2\f$ variate in the nominator of the F statistic.
   * \param df2 The non-negative degrees of freedom for the
   * \f$\chi^2\f$ variate in the denominator of the F statistic.
   *
   *
   * \see df(double x, double df1, double df2)
   * \see rng::rf(double df1, double df2)
   * 
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   * \throw scythe_convergence_error (Level 1)
   */
  inline double
  pf(double x, double df1, double df2)
  {
    SCYTHE_CHECK_10(df1 <= 0 || df2 <= 0, scythe_invalid_arg, 
        "df1 or df2 <= 0");
  
    if (x <= 0)
      return 0;
    
    if (df2 > 4e5)
      return pchisq(x*df1,df1);
    if (df1 > 4e5)
      return 1-pchisq(df2/x,df2);
    
    return (1-pbeta(df2 / (df2 + df1 * x), df2 / 2.0, df1 / 2.0));
  }

  SCYTHE_DISTFUN_MATRIX(pf, double, SCYTHE_ARGSET(df1, df2), double df1, double df2)

  /* PDFs */

  /*! \brief The F density function.
   *
   * Computes the value of the F probability density function with
   * \a df1 and \a df2 degrees of freedom, at the desired quantile \a
   * x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param df1 The non-negative degrees of freedom for the
   * \f$\chi^2\f$ variate in the nominator of the F statistic.
   * \param df2 The non-negative degrees of freedom for the
   * \f$\chi^2\f$ variate in the denominator of the F statistic.
   *
   * \see df(double x, double df1, double df2)
   * \see rng::rf(double df1, double df2)
   * 
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   * \throw scythe_convergence_error (Level 1)
   */
  inline double
  df(double x, double df1, double df2)
  {
    double dens;
    
    SCYTHE_CHECK_10(df1 <= 0 || df2 <= 0, scythe_invalid_arg, 
        "df1 or df2 <= 0");
    
    if (x <= 0)
      return 0;
      
    double f = 1 / (df2 + x * df1);
    double q = df2 * f;
    double p = x * df1 * f;
    
    if (df1 >= 2) {
      f = df1 * q / 2;
      dens = dbinom_raw((df1 - 2) / 2,(df1 + df2 - 2) / 2, p, q);
    } else {
      f = (df1 * df1 * q) /(2 * p * (df1 + df2));
      dens = dbinom_raw(df1 / 2,(df1 + df2)/ 2, p, q);
    }
    
    return f*dens;
  }

  SCYTHE_DISTFUN_MATRIX(df, double, SCYTHE_ARGSET(df1, df2), double df1, double df2)

  /**** The Gamma Distribution ****/

  /* CDFs */
  /*! \brief The gamma distribution function.
   *
   * Computes the value of the gamma cumulative distribution
   * function with given \a shape and \a scale, at the desired quantile
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param shape The non-negative shape of the distribution.
   * \param scale The non-negative scale of the distribution.
   *
   * \see dgamma(double x, double shape, double scale)
   * \see rng::rgamma(double shape, double scale)
   * \see gammafn(double x)
   * \see lngammafn(double x)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   * \throw scythe_convergence_error (Level 1)
   */
  inline double
  pgamma (double x, double shape, double scale)
  {
    const double xbig = 1.0e+8, xlarge = 1.0e+37, 
      alphlimit = 1000.;/* normal approx. for shape > alphlimit */
      
    int lower_tail = 1;

    double pn1, pn2, pn3, pn4, pn5, pn6, arg, a, b, c, an, osum, sum;
    long n;
    int pearson;

    /* check that we have valid values for x and shape */

    SCYTHE_CHECK_10(shape <= 0. || scale <= 0., scythe_invalid_arg,
        "shape or scale <= 0");

    x /= scale;
    
    if (x <= 0.)
      return 0.0;

    /* use a normal approximation if shape > alphlimit */

    if (shape > alphlimit) {
      pn1 = std::sqrt(shape) * 3. * (std::pow(x/shape, 1./3.) + 1.
            / (9. * shape) - 1.);
      return pnorm(pn1, 0., 1.);
    }

    /* if x is extremely large __compared to shape__ then return 1 */

    if (x > xbig * shape)
      return 1.0;

    if (x <= 1. || x < shape) {
      pearson = 1;/* use pearson's series expansion. */
      arg = shape * std::log(x) - x - lngammafn(shape + 1.);
      c = 1.;
      sum = 1.;
      a = shape;
      do {
        a += 1.;
        c *= x / a;
        sum += c;
      } while (c > DBL_EPSILON);
      arg += std::log(sum);
    }
    else { /* x >= max( 1, shape) */
      pearson = 0;/* use a continued fraction expansion */
      arg = shape * std::log(x) - x - lngammafn(shape);
      a = 1. - shape;
      b = a + x + 1.;
      pn1 = 1.;
      pn2 = x;
      pn3 = x + 1.;
      pn4 = x * b;
      sum = pn3 / pn4;
      for (n = 1; ; n++) {
        a += 1.;/* =   n+1 -shape */
        b += 2.;/* = 2(n+1)-shape+x */
        an = a * n;
        pn5 = b * pn3 - an * pn1;
        pn6 = b * pn4 - an * pn2;
        if (std::fabs(pn6) > 0.) {
          osum = sum;
          sum = pn5 / pn6;
          if (std::fabs(osum - sum) <= DBL_EPSILON * std::min(1., sum))
            break;
        }
        pn1 = pn3;
        pn2 = pn4;
        pn3 = pn5;
        pn4 = pn6;
        if (std::fabs(pn5) >= xlarge) {
          /* re-scale terms in continued fraction if they are large */
          pn1 /= xlarge;
          pn2 /= xlarge;
          pn3 /= xlarge;
          pn4 /= xlarge;
        }
      }
      arg += std::log(sum);
    }

    lower_tail = (lower_tail == pearson);

    sum = std::exp(arg);

    return (lower_tail) ? sum : 1 - sum;
  }

  SCYTHE_DISTFUN_MATRIX(pgamma, double, SCYTHE_ARGSET(shape, scale), double shape, double scale)

  /* PDFs */
  /*! \brief The gamma density function.
   *
   * Computes the value of the gamma probability density
   * function with given \a shape and \a scale, at the desired quantile
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param shape The non-negative shape of the distribution.
   * \param scale The non-negative scale of the distribution.
   *
   * \see pgamma(double x, double shape, double scale)
   * \see rng::rgamma(double shape, double scale)
   * \see gammafn(double x)
   * \see lngammafn(double x)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   * \throw scythe_convergence_error (Level 1)
   */
  inline double
  dgamma(double x, double shape, double scale)
  {
    SCYTHE_CHECK_10(shape <= 0 || scale <= 0,scythe_invalid_arg,
        "shape or scale <= 0");

    if (x < 0)
      return 0.0;
    
    if (x == 0) {
      SCYTHE_CHECK_10(shape < 1,scythe_invalid_arg, 
          "x == 0 and shape < 1");
      
      if (shape > 1)
        return 0.0;
      
      return 1 / scale;
    }
    
    if (shape < 1) { 
      double pr = dpois_raw(shape, x/scale);
      return pr * shape / x;
    }
    
    /* else  shape >= 1 */
    double pr = dpois_raw(shape - 1, x / scale);
    return pr / scale;
  }

  SCYTHE_DISTFUN_MATRIX(dgamma, double, SCYTHE_ARGSET(shape, scale), double shape, double scale)

  /**** The Logistic Distribution ****/

  /* CDFs */
  /*! \brief The logistic distribution function.
   *
   * Computes the value of the logistic cumulative distribution
   * function with given \a location and \a scale, at the desired
   * quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param location The location of the distribution.
   * \param scale The positive scale of the distribution.
   *
   * \see dlogis(double x, double location, double scale)
   * \see rng::rlogis(double location, double scale)
   * 
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  plogis (double x, double location, double scale)
  {
    SCYTHE_CHECK_10(scale <= 0.0, scythe_invalid_arg, "scale <= 0");
    
    double X = (x-location) / scale;
      
    X = std::exp(-X);
      
    return 1 / (1+X);
  }

  SCYTHE_DISTFUN_MATRIX(plogis, double, SCYTHE_ARGSET(location, scale), double location, double scale)

  /* PDFs */
  /*! \brief The logistic density function.
   *
   * Computes the value of the logistic probability density
   * function with given \a location and \a scale, at the desired
   * quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param location The location of the distribution.
   * \param scale The positive scale of the distribution.
   *
   * \see plogis(double x, double location, double scale)
   * \see rng::rlogis(double location, double scale)
   * 
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  dlogis(double x, double location, double scale)
  {
    SCYTHE_CHECK_10(scale <= 0.0, scythe_invalid_arg, "scale <= 0");
    
    double X = (x - location) / scale;
    double e = std::exp(-X);
    double f = 1.0 + e;
      
    return e / (scale * f * f);
  }

  SCYTHE_DISTFUN_MATRIX(dlogis, double, SCYTHE_ARGSET(location, scale), double location, double scale)

  /**** The Log Normal Distribution ****/

  /* CDFs */
  /*! \brief The log-normal distribution function.
   *
   * Computes the value of the log-normal cumulative distribution
   * function with mean \a logmean and standard
   * deviation \a logsd, at the desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param logmean The mean of the distribution.
   * \param logsd The positive standard deviation of the distribution.
   *
   * \see dlnorm(double x, double logmean, double logsd)
   * \see rng::rlnorm(double logmean, double logsd)
   * \see pnorm(double x, double logmean, double logsd)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_convergence_error (Level 1)
   */
  inline double
  plnorm (double x, double logmean, double logsd)
  {
    SCYTHE_CHECK_10(logsd <= 0, scythe_invalid_arg, "logsd <= 0");
    
    if (x > 0)
      return pnorm(std::log(x), logmean, logsd);
    
    return 0;
  }

  SCYTHE_DISTFUN_MATRIX(plnorm, double, SCYTHE_ARGSET(logmean, logsd), double logmean, double logsd)

  /* PDFs */
  /*! \brief The log-normal density function.
   *
   * Computes the value of the log-normal probability density
   * function with mean \a logmean and standard
   * deviation \a logsd, at the desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param logmean The mean of the distribution.
   * \param logsd The positive standard deviation of the distribution.
   *
   * \see plnorm(double x, double logmean, double logsd)
   * \see rng::rlnorm(double logmean, double logsd)
   * \see dnorm(double x, double logmean, double logsd)
   *
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  dlnorm(double x, double logmean, double logsd)
  {
    SCYTHE_CHECK_10(logsd <= 0, scythe_invalid_arg, "logsd <= 0");
    
    if (x == 0)
      return 0;
    
    double y = (std::log(x) - logmean) / logsd;
    
    return (1 / (std::sqrt(2 * M_PI))) * std::exp(-0.5 * y * y) / (x * logsd);
  }

  SCYTHE_DISTFUN_MATRIX(dlnorm, double, SCYTHE_ARGSET(logmean, logsd), double logmean, double logsd)

  /**** The Negative Binomial Distribution ****/

  /* CDFs */
  /*! \brief The negative binomial distribution function.
   *
   * Computes the value of the negative binomial cumulative distribution
   * function with \a n target number of successful trials and \a p
   * probability of success on each trial, at the desired quantile \a
   * x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired non-negative, integer, quantile.
   * \param n The positive target number of successful trials
   * (dispersion parameter).
   * \param p The probability of success on each trial.
   *
   * \see dnbinom(unsigned int x, double n, double p)
   * \see rng::rnbinom(double n, double p)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  pnbinom(unsigned int x, double n, double p)
  {
    SCYTHE_CHECK_10(n == 0 || p <= 0 || p >= 1, scythe_invalid_arg,
        "n == 0 or p not in (0,1)");
    
    return pbeta(p, n, x + 1);
  }

  SCYTHE_DISTFUN_MATRIX(pnbinom, unsigned int, SCYTHE_ARGSET(n, p), double n, double p)

  /* PDFs */
  /*! \brief The negative binomial density function.
   *
   * Computes the value of the negative binomial probability density
   * function with \a n target number of successful trials and \a p
   * probability of success on each trial, at the desired quantile \a
   * x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired non-negative, integer, quantile.
   * \param n The positive target number of successful trials
   * (dispersion parameter).
   * \param p The probability of success on each trial.
   *
   * \see dnbinom(unsigned int x, double n, double p)
   * \see rng::rnbinom(double n, double p)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  dnbinom(unsigned int x, double n, double p)
  {
    SCYTHE_CHECK_10(n == 0 || p <= 0 || p >= 1, scythe_invalid_arg,
        "n == 0 or p not in (0,1)");
    
    double prob = dbinom_raw(n, x + n, p, 1 - p);
    double P = (double) n / (n + x);
    
    return P * prob;
  }

  SCYTHE_DISTFUN_MATRIX(dnbinom, unsigned int, SCYTHE_ARGSET(n, p), double n, double p)

  /**** The Normal Distribution ****/
  
  /* CDFs */
  /*! \brief The normal distribution function.
   *
   * Computes the value of the normal cumulative distribution
   * function with given \a mean and standard deviation \a sd, at the
   * desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param mean The mean of the distribution.
   * \param sd The positive standard deviation of the distribution.
   *
   * \see dnorm(double x, double mean, double sd)
   * \see rng::rnorm(double mean, double sd)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_convergence_error (Level 1)
   */
  inline double
  pnorm (double x, double mean, double sd)
  
  {
    SCYTHE_CHECK_10(sd <= 0, scythe_invalid_arg,
        "negative standard deviation");

    return pnorm1((x - mean) / sd, true, false);
  }

  SCYTHE_DISTFUN_MATRIX(pnorm, double, SCYTHE_ARGSET(mean, sd), double mean, double sd)
  

  /* PDFs */
  /*! \brief The normal density function.
   *
   * Computes the value of the normal probability density
   * function with given \a mean and standard deviation \a sd, at the
   * desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param mean The mean of the distribution.
   * \param sd The positive standard deviation of the distribution.
   *
   * \see pnorm(double x, double mean, double sd)
   * \see rng::rnorm(double mean, double sd)
   *
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  dnorm(double x, double mean, double sd)
  {
    SCYTHE_CHECK_10(sd <= 0, scythe_invalid_arg,
        "negative standard deviation");
    
    double X = (x - mean) / sd;
    
    return (M_1_SQRT_2PI * std::exp(-0.5 * X * X) / sd);
  }

  SCYTHE_DISTFUN_MATRIX(dnorm, double, SCYTHE_ARGSET(mean, sd), double mean, double sd)
 

  /* Return the natural log of the normal PDF */
  /*! \brief The natural log of normal density function.
   *
   * Computes the value of the natural log of the normal probability
   * density function with given \a mean and standard deviation \a sd,
   * at the desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param mean The mean of the distribution.
   * \param sd The positive standard deviation of the distribution.
   *
   * \see dnorm(double x, double mean, double sd)
   * \see pnorm(double x, double mean, double sd)
   * \see rng::rnorm(double mean, double sd)
   *
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  lndnorm (double x, double mean, double sd)
  {
    SCYTHE_CHECK_10(sd <= 0, scythe_invalid_arg,
        "negative standard deviation");
    
    double X = (x - mean) / sd;
    
    return -(M_LN_SQRT_2PI  +  0.5 * X * X + std::log(sd));
  }

  SCYTHE_DISTFUN_MATRIX(lndnorm, double, SCYTHE_ARGSET(mean, sd), double mean, double sd)

  /* Quantile functions */
  /*! \brief The standard normal quantile function.
   *
   * Computes the value of the standard normal quantile function
   * at the desired probability \a in_p.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param in_p The desired probability.
   *
   * \see pnorm(double x, double mean, double sd)
   * \see dnorm(double x, double mean, double sd)
   * \see rng::rnorm(double mean, double sd)
   *
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  qnorm1 (double in_p)
  {
    double p0 = -0.322232431088;
    double q0 = 0.0993484626060;
    double p1 = -1.0;
    double q1 = 0.588581570495;
    double p2 = -0.342242088547;
    double q2 = 0.531103462366;
    double p3 = -0.0204231210245;
    double q3 = 0.103537752850;
    double p4 = -0.453642210148e-4;
    double q4 = 0.38560700634e-2;
    double xp = 0.0;
    double p = in_p;
      
    if (p > 0.5)
      p = 1 - p;
        
    SCYTHE_CHECK_10(p < 10e-20, scythe_range_error,
        "p outside accuracy limit");
      
    if (p == 0.5)
      return xp;
      
    double y = std::sqrt (std::log (1.0 / std::pow (p, 2)));
    xp = y + ((((y * p4 + p3) * y + p2) * y + p1) * y + p0) /
      ((((y * q4 + q3) * y + q2) * y + q1) * y + q0);
      
    if (in_p < 0.5)
      xp = -1 * xp;
    
    return xp;
  }

  SCYTHE_DISTFUN_MATRIX(qnorm1, double, in_p, double in_p)

  /**** The Poisson Distribution ****/

  /* CDFs */
  /*! \brief The Poisson distribution function.
   *
   * Computes the value of the Poisson cumulative distribution
   * function with expected number of occurrences \a lambda, at the
   * desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired integer quantile.
   * \param lambda The expected positive number of occurrences.
   *
   * \see dpois(unsigned int x, double lambda)
   * \see rng::rpois(double lambda)
   *
   * \throws scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   * \throw scythe_convergence_error (Level 1)
   */
  inline double
  ppois(unsigned int x, double lambda)
  {
    SCYTHE_CHECK_10(lambda<=0.0, scythe_invalid_arg, "lambda <= 0");
    
    if (lambda == 1)
      return 1;
    
    return 1 - pgamma(lambda, x + 1, 1.0);
  }

  SCYTHE_DISTFUN_MATRIX(ppois, unsigned int, lambda, double lambda)

  /* PDFs */
  /*! \brief The Poisson density function.
   *
   * Computes the value of the Poisson probability density
   * function with expected number of occurrences \a lambda, at the
   * desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired integer quantile.
   * \param lambda The expected positive number of occurrences.
   *
   * \see ppois(unsigned int x, double lambda)
   * \see rng::rpois(double lambda)
   *
   * \throws scythe_invalid_arg (Level 1)
   */
  inline double
  dpois(unsigned int x, double lambda)
  {
    SCYTHE_CHECK_10(lambda<=0.0, scythe_invalid_arg, "lambda <= 0");
    
    // compute log(x!)
    double xx = x+1;
    double cof[6] = {
      76.18009172947146, -86.50532032941677,
      24.01409824083091, -1.231739572450155,
      0.1208650973866179e-2, -0.5395239384953e-5
    };
    double y = xx;
    double tmp = xx + 5.5 - (xx + 0.5) * std::log(xx + 5.5);
    double ser = 1.000000000190015;
    for (int j = 0; j <= 5; j++) {
      ser += (cof[j] / ++y);
    }
    double lnfactx = std::log(2.5066282746310005 * ser / xx) - tmp;
      
    return (std::exp( -1*lnfactx + x * std::log(lambda) - lambda));
  }

  SCYTHE_DISTFUN_MATRIX(dpois, unsigned int, lambda, double lambda)

  /**** The t Distribution ****/

  /* CDFs */
  /*! \brief The Student t distribution function.
   *
   * Computes the value of the Student t cumulative distribution
   * function with \a n degrees of freedom, at the desired quantile
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param n The positive degrees of freedom of the distribution.
   *
   * \see dt(double x, bool b1, bool b2)
   * \see rng::rt1(double mu, double sigma2, double nu)
   * 
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_convergence_error (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  pt(double x, double n)
  {
    double val;
    
    SCYTHE_CHECK_10(n <= 0, scythe_invalid_arg, "n <= 0");
    
    if (n > 4e5) {
      val = 1/(4*n);
      return pnorm1(x * (1 - val) / ::sqrt(1 + x * x * 2. * val), 
          true, false);
    }
    
    val = pbeta(n / (n + x * x), n / 2.0, 0.5);
    
    val /= 2;
    
    if (x <= 0)
      return val;
    else
      return 1 - val;
  }

  SCYTHE_DISTFUN_MATRIX(pt, double, n, double n)
  
  /* PDFs */
  /*! \brief The Student t distribution function.
   *
   * Computes the value of the Student t cumulative distribution
   * function with \a n degrees of freedom, at the desired quantile
   * \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param n The positive degrees of freedom of the distribution.
   *
   * \see pt(double x, bool b1, bool b2)
   * \see rng::rt1(double mu, double sigma2, double nu)
   * 
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  dt(double x, double n)
  {
    double u;

    SCYTHE_CHECK_10(n <= 0, scythe_invalid_arg, "n <= 0");
    
    double t = -bd0(n/2., (n + 1) / 2.)
      + stirlerr((n + 1) / 2.)
      - stirlerr(n / 2.);
    if(x*x > 0.2*n)
      u = std::log(1+x*x/n)*n/2;
    else
      u = -bd0(n/2., (n+x*x)/2.) + x*x/2;
    
    return std::exp(t-u)/std::sqrt(2*M_PI*(1+x*x/n));
  }

  SCYTHE_DISTFUN_MATRIX(dt, double, n, double n)
  
  /* Returns the univariate Student-t density evaluated at x 
   * with mean mu, scale sigma^2, and nu degrees of freedom.  
   *
   * TODO:  Do we want a pt for this distribution?
   */

  /*! \brief The univariate Student t density function.
   *
   * Computes the value of the univariate Student t probability
   * density function with mean \a mu, variance \a sigma2,
   * and degrees of freedom \a nu, at the desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param mu The mean of the distribution.
   * \param sigma2 The variance of the distribution.
   * \param nu The degrees of freedom of the distribution.
   *
   * \see rng::rt1(double mu, double sigma2, double nu)
   * \see dt(double x, bool b1, bool b2)
   * \see pt(double x, bool b1, bool b2)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double
  dt1(double x, double mu, double sigma2, double nu)
  {
    double logdens =   lngammafn((nu + 1.0) /2.0)
      - std::log(std::sqrt(nu * M_PI))
      - lngammafn(nu / 2.0) - std::log(std::sqrt(sigma2))
      - (nu + 1.0) / 2.0 * std::log(1.0 + (std::pow((x - mu), 2.0))
            / (nu * sigma2));
    
    return(std::exp(logdens));
  }

  SCYTHE_DISTFUN_MATRIX(dt1, double, SCYTHE_ARGSET(mu, sigma2, nu), double mu, double sigma2, double nu)

  /* Returns the natural log of the univariate Student-t density 
   * evaluated at x with mean mu, scale sigma^2, and nu 
   * degrees of freedom
   */
   
  /*! \brief The natural log of the univariate Student t density
   * function.
   *
   * Computes the value of the natural log of the univariate Student t
   * probability density function with mean \a mu, variance \a sigma2,
   * and degrees of freedom \a nu, at the desired quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param mu The mean of the distribution.
   * \param sigma2 The variance of the distribution.
   * \param nu The degrees of freedom of the distribution.
   *
   * \see rng::rt1(double mu, double sigma2, double nu)
   * \see dt(double x, bool b1, bool b2)
   * \see pt(double x, bool b1, bool b2)
   *
   * \throw scythe_invalid_arg (Level 1)
   * \throw scythe_range_error (Level 1)
   * \throw scythe_precision_error (Level 1)
   */
  inline double 
  lndt1(double x, double mu, double sigma2, double nu)
  {
    double logdens = lngammafn((nu+1.0)/2.0)
      - std::log(std::sqrt(nu*M_PI))
      - lngammafn(nu/2.0) - std::log(std::sqrt(sigma2))
      - (nu+1.0)/2.0 * std::log(1.0 + (std::pow((x-mu),2.0))
        /(nu * sigma2));
    
    return(logdens);
  }

  SCYTHE_DISTFUN_MATRIX(lndt1, double, SCYTHE_ARGSET(mu, sigma2, nu), double mu, double sigma2, double nu)

  /**** The Uniform Distribution ****/

  /* CDFs */
  /*! \brief The uniform distribution function.
   *
   * Computes the value of the uniform cumulative distribution
   * function evaluated on the interval [\a a, \a b], at the desired
   * quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile x.
   * \param a The lower end-point of the distribution.
   * \param b The upper end-point of the distribution.
   *
   * \see dunif(double x, double a, double b)
   * \see rng::runif()
   *
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  punif(double x, double a, double b)
  {
    SCYTHE_CHECK_10(b <= a, scythe_invalid_arg, "b <= a");
      
    if (x <= a)
      return 0.0;
        
    if (x >= b)
      return 1.0;
      
    return (x - a) / (b - a);
  }

  SCYTHE_DISTFUN_MATRIX(punif, double, SCYTHE_ARGSET(a, b), double a, double b)

  /* PDFs */
  /*! \brief The uniform density function.
   *
   * Computes the value of the uniform probability density
   * function evaluated on the interval [\a a, \a b], at the desired
   * quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile x.
   * \param a The lower end-point of the distribution.
   * \param b The upper end-point of the distribution.
   *
   * \see punif(double x, double a, double b)
   * \see rng::runif()
   *
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  dunif(double x, double a, double b)
  {
    SCYTHE_CHECK_10(b <= a, scythe_invalid_arg, "b <= a");
    
    if (a <= x && x <= b)
      return 1.0 / (b - a);
    
    return 0.0;
  }

  SCYTHE_DISTFUN_MATRIX(dunif, double, SCYTHE_ARGSET(a, b), double a, double b)

  /**** The Weibull Distribution ****/

  /* CDFs */
  /*! \brief The Weibull distribution function.
   *
   * Computes the value of the Weibull cumulative distribution
   * function with given \a shape and \a scale, at the desired
   * quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param shape The positive shape of the distribution.
   * \param scale The positive scale of the distribution.
   *
   * \see dweibull(double x, double shape, double scale)
   * \see rng::rweibull(double shape, double scale)
   *
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  pweibull(double x, double shape, double scale)
  {
    SCYTHE_CHECK_10(shape <= 0 || scale <= 0, scythe_invalid_arg,
        "shape or scale <= 0");
    
    if (x <= 0)
      return 0.0;
    
    return 1 - std::exp(-std::pow(x / scale, shape));
  }

  SCYTHE_DISTFUN_MATRIX(pweibull, double, SCYTHE_ARGSET(shape, scale), double shape, double scale)

  /* PDFs */
  /*! \brief The Weibull density function.
   *
   * Computes the value of the Weibull probability density
   * function with given \a shape and \a scale, at the desired
   * quantile \a x.
   *
   * It is also possible to call this function with a Matrix of
   * doubles as its first argument.  In this case the function will
   * return a Matrix of doubles of the same dimension as \a x,
   * containing the result of evaluating this function at each value
   * in \a x, given the remaining fixed parameters.  By default, the
   * returned Matrix will be concrete and have the same matrix_order
   * as \a x, but you may invoke a generalized version of the function
   * with an explicit template call.
   *
   * \param x The desired quantile.
   * \param shape The positive shape of the distribution.
   * \param scale The positive scale of the distribution.
   *
   * \see pweibull(double x, double shape, double scale)
   * \see rng::rweibull(double shape, double scale)
   *
   * \throw scythe_invalid_arg (Level 1)
   */
  inline double
  dweibull(double x, double shape, double scale)
  {
    SCYTHE_CHECK_10(shape <= 0 || scale <= 0, scythe_invalid_arg,
        "shape or scale <= 0");

    if (x < 0)
      return 0.;
      
    double tmp1 = std::pow(x / scale, shape - 1);
    double tmp2 = tmp1*(x / scale);
      
    return shape * tmp1 * std::exp(-tmp2) / scale;
  }

  SCYTHE_DISTFUN_MATRIX(dweibull, double, SCYTHE_ARGSET(shape, scale), double shape, double scale)

  /* Multivariate Normal */

  // TODO: distribution function.  Plain old (non-logged) dmvnorm.

  /*! \brief The natural log of the multivariate normal density
   * function.
   *
   * Computes the value of the natural log of the multivariate normal
   * probability density function with vector of mean \a mu and
   * variance-covariance matrix \a Sigma, at the vector of desired
   * quantiles \a x.
   *
   * \param x The vector of desired quantiles.
   * \param mu The vector of means.
   * \param Sigma The variance-covariance matrix.
   *
   * \see rng:rmvnorm(const Matrix<double, PO1, PS1>& mu, const Matrix<double, PO2, PS2>& sigma)
   *
   * \throw scythe_dimension_error (Level 1)
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_null_error (Level 1)
   */
  template <matrix_order O1, matrix_style S1,
            matrix_order O2, matrix_style S2,
            matrix_order O3, matrix_style S3>
  double lndmvn (const Matrix<double, O1, S1>& x,
                 const Matrix<double, O2, S2>& mu,
                 const Matrix<double, O3, S3>& Sigma)
  {
    SCYTHE_CHECK_10(! x.isColVector(), scythe_dimension_error,
        "x is not a column vector");
    SCYTHE_CHECK_10(! mu.isColVector(), scythe_dimension_error,
        "mu is not a column vector");
    SCYTHE_CHECK_10(! Sigma.isSquare(), scythe_dimension_error,
        "Sigma is not square");
    SCYTHE_CHECK_10(mu.rows()!=Sigma.rows() || x.rows()!=Sigma.rows(), 
                    scythe_conformation_error,
                    "mu, x and Sigma have mismatched row lengths")
    int k = (int) mu.rows();
    return ( (-k/2.0)*std::log(2*M_PI) -0.5 * std::log(det(Sigma)) 
       -0.5 * (t(x - mu)) * invpd(Sigma) * (x-mu) )[0];
  }

} // end namespace scythe


#endif /* SCYTHE_DISTRIBUTIONS_H */
