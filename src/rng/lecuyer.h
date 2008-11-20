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
 * scythestat/rng/lecuyer.h
 *
 * Provides the class definition for the L'Ecuyer random number
 * generator, a rng capable of generating many independent substreams.
 * This class extends the abstract rng class by implementing runif().
 * Based on RngStream.cpp, by Pierre L'Ecuyer.
 *
 * Pierre L'Ecuyer agreed to the following dual-licensing terms in an
 * email received 7 August 2004.  This dual-license was prompted by
 * the Debian maintainers of R and MCMCpack. 
 *
 * This software is Copyright (C) 2004 Pierre L'Ecuyer.
 *
 * License: this code can be used freely for personal, academic, or
 * non-commercial purposes.  For commercial licensing, please contact
 * P. L'Ecuyer at lecuyer@iro.umontreal.ca.
 *
 * This code may also be redistributed and modified under the terms of
 * the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA.
 *
 */
/*! \file lecuyer.h
 * \brief The L'Ecuyer random number generator.
 *
 * This file contains the lecuyer class, a class that extends Scythe's
 * base random number generation class (scythe::rng) by providing an
 * implementation of scythe::rng::runif(), using L'Ecuyer's algorithm.
 *
 */
#ifndef SCYTHE_LECUYER_H
#define SCYTHE_LECUYER_H

#include<cstdlib>
#include<iostream>
#include<string>

#ifdef SCYTHE_COMPILE_DIRECT
#include "rng.h"
#else
#include "scythestat/rng.h"
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

	SCYTHE_MINGW32_STATIC const double m1   = 4294967087.0;
	SCYTHE_MINGW32_STATIC const double m2   = 4294944443.0;
	SCYTHE_MINGW32_STATIC const double norm = 1.0 / (m1 + 1.0);
	SCYTHE_MINGW32_STATIC const double a12  = 1403580.0;
	SCYTHE_MINGW32_STATIC const double a13n = 810728.0;
	SCYTHE_MINGW32_STATIC const double a21  = 527612.0;
	SCYTHE_MINGW32_STATIC const double a23n = 1370589.0;
	SCYTHE_MINGW32_STATIC const double two17 =131072.0;
	SCYTHE_MINGW32_STATIC const double two53 =9007199254740992.0;
  /* 1/2^24 */
	SCYTHE_MINGW32_STATIC const double fact = 5.9604644775390625e-8;

	// The following are the transition matrices of the two MRG
	// components (in matrix form), raised to the powers -1, 1, 2^76,
	// and 2^127, resp.

	SCYTHE_MINGW32_STATIC const double InvA1[3][3] = { // Inverse of A1p0
				 { 184888585.0,   0.0,  1945170933.0 },
				 {         1.0,   0.0,           0.0 },
				 {         0.0,   1.0,           0.0 } };

	SCYTHE_MINGW32_STATIC const double InvA2[3][3] = { // Inverse of A2p0
				 {      0.0,  360363334.0,  4225571728.0 },
				 {      1.0,          0.0,           0.0 },
				 {      0.0,          1.0,           0.0 } };

	SCYTHE_MINGW32_STATIC const double A1p0[3][3] = {
				 {       0.0,        1.0,       0.0 },
				 {       0.0,        0.0,       1.0 },
				 { -810728.0,  1403580.0,       0.0 } };

	SCYTHE_MINGW32_STATIC const double A2p0[3][3] = {
				 {        0.0,        1.0,       0.0 },
				 {        0.0,        0.0,       1.0 },
				 { -1370589.0,        0.0,  527612.0 } };

	SCYTHE_MINGW32_STATIC const double A1p76[3][3] = {
				 {      82758667.0, 1871391091.0, 4127413238.0 },
				 {    3672831523.0,   69195019.0, 1871391091.0 },
				 {    3672091415.0, 3528743235.0,   69195019.0 } };

	SCYTHE_MINGW32_STATIC const double A2p76[3][3] = {
				 {    1511326704.0, 3759209742.0, 1610795712.0 },
				 {    4292754251.0, 1511326704.0, 3889917532.0 },
				 {    3859662829.0, 4292754251.0, 3708466080.0 } };

	SCYTHE_MINGW32_STATIC const double A1p127[3][3] = {
				 {    2427906178.0, 3580155704.0,  949770784.0 },
				 {     226153695.0, 1230515664.0, 3580155704.0 },
				 {    1988835001.0,  986791581.0, 1230515664.0 } };

	SCYTHE_MINGW32_STATIC const double A2p127[3][3] = {
				 {    1464411153.0,  277697599.0, 1610723613.0 },
				 {      32183930.0, 1464411153.0, 1022607788.0 },
				 {    2824425944.0,   32183930.0, 2093834863.0 } };

	// Return (a*s + c) MOD m; a, s, c and m must be < 2^35
	SCYTHE_MINGW32_STATIC double
	MultModM (double a, double s, double c, double m)
	{
		double v;
		long a1;

		v = a * s + c;

		if (v >= two53 || v <= -two53) {
			a1 = static_cast<long> (a / two17);    a -= a1 * two17;
			v  = a1 * s;
			a1 = static_cast<long> (v / m);     v -= a1 * m;
			v = v * two17 + a * s + c;
		}

		a1 = static_cast<long> (v / m);
		/* in case v < 0)*/
		if ((v -= a1 * m) < 0.0) return v += m;   else return v;
	}

	// Compute the vector v = A*s MOD m. Assume that -m < s[i] < m.
	// Works also when v = s.
	SCYTHE_MINGW32_STATIC void
	MatVecModM (const double A[3][3], const double s[3],
							double v[3], double m)
	{
		int i;
		double x[3];               // Necessary if v = s

		for (i = 0; i < 3; ++i) {
			x[i] = MultModM (A[i][0], s[0], 0.0, m);
			x[i] = MultModM (A[i][1], s[1], x[i], m);
			x[i] = MultModM (A[i][2], s[2], x[i], m);
		}
		for (i = 0; i < 3; ++i)
			v[i] = x[i];
	}

	// Compute the matrix C = A*B MOD m. Assume that -m < s[i] < m.
	// Note: works also if A = C or B = C or A = B = C.
	SCYTHE_MINGW32_STATIC void
	MatMatModM (const double A[3][3], const double B[3][3],
							double C[3][3], double m)
	{
		int i, j;
		double V[3], W[3][3];

		for (i = 0; i < 3; ++i) {
			for (j = 0; j < 3; ++j)
				V[j] = B[j][i];
			MatVecModM (A, V, V, m);
			for (j = 0; j < 3; ++j)
				W[j][i] = V[j];
		}
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
				C[i][j] = W[i][j];
	}

	// Compute the matrix B = (A^(2^e) Mod m);  works also if A = B. 
	SCYTHE_MINGW32_STATIC void
	MatTwoPowModM(const double A[3][3], double B[3][3],
								double m, long e)
	{
	 int i, j;

	 /* initialize: B = A */
	 if (A != B) {
		 for (i = 0; i < 3; ++i)
			 for (j = 0; j < 3; ++j)
				 B[i][j] = A[i][j];
	 }
	 /* Compute B = A^(2^e) mod m */
	 for (i = 0; i < e; i++)
		 MatMatModM (B, B, B, m);
	}

	// Compute the matrix B = (A^n Mod m);  works even if A = B.
	SCYTHE_MINGW32_STATIC void
	MatPowModM (const double A[3][3], double B[3][3], double m,
							long n)
	{
		int i, j;
		double W[3][3];

		/* initialize: W = A; B = I */
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j) {
				W[i][j] = A[i][j];
				B[i][j] = 0.0;
			}
		for (j = 0; j < 3; ++j)
			B[j][j] = 1.0;

		/* Compute B = A^n mod m using the binary decomposition of n */
		while (n > 0) {
			if (n % 2) MatMatModM (W, B, B, m);
			MatMatModM (W, W, W, m);
			n /= 2;
		}
	}

	// Check that the seeds are legitimate values. Returns 0 if legal
	// seeds, -1 otherwise.
	SCYTHE_MINGW32_STATIC int
	CheckSeed (const unsigned long seed[6])
	{
		int i;

		for (i = 0; i < 3; ++i) {
			if (seed[i] >= m1) {
      SCYTHE_THROW(scythe_randseed_error,
          "Seed[" << i << "] >= 4294967087, Seed is not set");
				return -1;
			}
		}
		for (i = 3; i < 6; ++i) {
			if (seed[i] >= m2) {
      SCYTHE_THROW(scythe_randseed_error,
          "Seed[" << i << "] >= 4294944443, Seed is not set");
				return -1;
			}
		}
		if (seed[0] == 0 && seed[1] == 0 && seed[2] == 0) {
      SCYTHE_THROW(scythe_randseed_error, "First 3 seeds = 0");
			return -1;
		}
		if (seed[3] == 0 && seed[4] == 0 && seed[5] == 0) {
      SCYTHE_THROW(scythe_randseed_error, "Last 3 seeds = 0");
			return -1;
		}

		return 0;
	}
 
#ifndef __MINGW32__
  } // end anonymous namespace
#endif

   /*! \brief The L'Ecuyer random number generator.
    *
    * This class defines a random number generator, using Pierre
    * L'Ecuyer's algorithm (2000) and source code (2001) for
    * generating multiple simultaneous streams of random uniform
    * variates.  The period of the underlying single-stream generator
    * is approximately \f$3.1 \times 10^{57}\f$.  Each individual
    * stream is implemented in terms of a sequence of substreams (see
    * L'Ecuyer et al (2000) for details).
    *
    * The lecuyer class extends Scythe's basic random number
    * generating class, scythe::rng, implementing the interface that
    * it defines.
    *
    * \see rng
    * \see mersenne
    *
    */
  class lecuyer : public rng<lecuyer>
  {
    public:

      // Constructor
      /*! \brief Constructor
       *
       * This constructor creates an object encapsulating a random
       * number stream, with an optional name.  It also sets the seed
       * of the stream to the package (default or user-specified) seed
       * if this is the first stream generated, or, otherwise, to a
       * value \f$2^{127}\f$ steps ahead of the seed of the previously
       * constructed stream.
       *
       * \param streamname The optional name for the stream.
       *
       * \see SetPackageSeed(unsigned long seed[6])
       * \see SetSeed(unsigned long seed[6])
       * \see SetAntithetic(bool)
       * \see IncreasedPrecis(bool)
       * \see name()
       */
      lecuyer (std::string streamname = "")
        : rng<lecuyer> (),
          streamname_ (streamname)
      {
        anti = false;
        incPrec = false;
        
        /* Information on a stream. The arrays {Cg, Bg, Ig} contain
         * the current state of the stream, the starting state of the
         * current SubStream, and the starting state of the stream.
         * This stream generates antithetic variates if anti = true.
         * It also generates numbers with extended precision (53 bits
         * if machine follows IEEE 754 standard) if incPrec = true.
         * nextSeed will be the seed of the next declared RngStream.
         */

        for (int i = 0; i < 6; ++i) {
          Bg[i] = Cg[i] = Ig[i] = nextSeed[i];
        }

        MatVecModM (A1p127, nextSeed, nextSeed, m1);
        MatVecModM (A2p127, &nextSeed[3], &nextSeed[3], m2);
      }

      /*! \brief Get the stream's name.
       *
       * This method returns a stream's name string.
       *
       * \see lecuyer(const char*)
       */
      std::string
      name() const
      {
        return streamname_;
      }

      /*! \brief Reset the stream.
       *
       * This method resets the stream to its initial seeded state.
       *
       * \see ResetStartSubstream()
       * \see ResetNextSubstream()
       * \see SetSeed(unsigned long seed[6])
       */
      void
      ResetStartStream ()
      {
        for (int i = 0; i < 6; ++i)
          Cg[i] = Bg[i] = Ig[i];
      }

      /*! \brief Reset the current substream.
       *
       *
       * This method resets the stream to the first state of its
       * current substream.
       *
       * \see ResetStartStream()
       * \see ResetNextSubstream()
       * \see SetSeed(unsigned long seed[6])
       * 
       */
      void
      ResetStartSubstream ()
      {
        for (int i = 0; i < 6; ++i)
          Cg[i] = Bg[i];
      }

      /*! \brief Jump to the next substream.
       *
       * This method resets the stream to the first state of its next
       * substream.
       *
       * \see ResetStartStream()
       * \see ResetStartSubstream()
       * \see SetSeed(unsigned long seed[6])
       * 
       */
      void
      ResetNextSubstream ()
      {
        MatVecModM(A1p76, Bg, Bg, m1);
        MatVecModM(A2p76, &Bg[3], &Bg[3], m2);
        for (int i = 0; i < 6; ++i)
          Cg[i] = Bg[i];
      }

      /*! \brief Set the package seed.
       *
       *  This method sets the overall package seed.  The default
       *  initial seed is (12345, 12345, 12345, 12345, 12345, 12345).
       *  The package seed is the seed used to initialize the first
       *  constructed random number stream in a given program.
       *
       *  \param seed An array of six integers to seed the package.
       *  The first three values cannot all equal 0 and must all be
       *  less than 4294967087 while the second trio of integers must
       *  all be less than 4294944443 and not all 0.
       *
       * \see SetSeed(unsigned long seed[6])
       *
       * \throw scythe_randseed_error (Level 0)
       */
      static void
      SetPackageSeed (unsigned long seed[6])
      {
         if (CheckSeed (seed)) return;
         for (int i = 0; i < 6; ++i)
            nextSeed[i] = seed[i];
      }

      /*! \brief Set the stream seed.
       *
       *  This method sets the stream seed which is used to initialize
       *  the state of the given stream.
       *
       *  \warning This method sets the stream seed in isolation and
       *  does not coordinate with any other streams.  Therefore,
       *  using this method without care can result in multiple
       *  streams that overlap in the course of their runs.
       *
       *  \param seed An array of six integers to seed the stream.
       *  The first three values cannot all equal 0 and must all be
       *  less than 4294967087 while the second trio of integers must
       *  all be less than 4294944443 and not all 0.
       *
       * \see SetPackageSeed(unsigned long seed[6])
       * \see ResetStartStream()
       * \see ResetStartSubstream()
       * \see ResetNextSubstream()
       *
       * \throw scythe_randseed_error (Level 0)
       */
      void
      SetSeed (unsigned long seed[6])
      {
        if (CheckSeed (seed)) return;
          for (int i = 0; i < 6; ++i)
            Cg[i] = Bg[i] = Ig[i] = seed[i];
      }

      // XXX: get the cases formula working!
      /*! \brief Advances the state of the stream.
       *
       * This method advances the input \f$n\f$ steps, using the rule:
       * \f[
       * n =
       * \begin{cases}
       *  2^e + c \quad if~e > 0, \\
       *  -2^{-e} + c \quad if~e < 0, \\
       *  c \quad if~e = 0.
       * \end{cases}
       * \f]
       * 
       * \param e This parameter controls state advancement.
       * \param c This parameter also controls state advancement.
       *
       * \see GetState()
       * \see ResetStartStream()
       * \see ResetStartSubstream()
       * \see ResetNextSubstream()
       */
      void
      AdvanceState (long e, long c)
      {
        double B1[3][3], C1[3][3], B2[3][3], C2[3][3];

        if (e > 0) {
          MatTwoPowModM (A1p0, B1, m1, e);
          MatTwoPowModM (A2p0, B2, m2, e);
        } else if (e < 0) {
          MatTwoPowModM (InvA1, B1, m1, -e);
          MatTwoPowModM (InvA2, B2, m2, -e);
        }

        if (c >= 0) {
          MatPowModM (A1p0, C1, m1, c);
          MatPowModM (A2p0, C2, m2, c);
        } else {
          MatPowModM (InvA1, C1, m1, -c);
          MatPowModM (InvA2, C2, m2, -c);
        }

        if (e) {
          MatMatModM (B1, C1, C1, m1);
          MatMatModM (B2, C2, C2, m2);
        }

        MatVecModM (C1, Cg, Cg, m1);
        MatVecModM (C2, &Cg[3], &Cg[3], m2);
      }

      /*! \brief Get the current state.
       *
       * This method places the current state of the stream, as
       * represented by six integers, into the array argument.  This
       * is useful for saving and restoring streams across program
       * runs.
       *
       * \param seed An array of six integers that will hold the state values on return.
       *
       * \see AdvanceState()
       */
      void
      GetState (unsigned long seed[6]) const
      {
        for (int i = 0; i < 6; ++i)
          seed[i] = static_cast<unsigned long> (Cg[i]);
      }

      /*! \brief Toggle generator precision.
       *
       * This method sets the precision level of the given stream.  By
       * default, streams generate random numbers with 32 bit
       * resolution.  If the user invokes this method with \a incp =
       * true, then the stream will begin to generate variates with
       * greater precision (53 bits on machines following the IEEE 754
       * standard).  Calling this method again with \a incp = false
       * will return the precision of generated numbers to 32 bits.
       *
       * \param incp A boolean value where true implies high (most
       * likely 53 bit) precision and false implies low (32 bit)
       * precision.
       *
       * \see SetAntithetic(bool)
       */
      void
      IncreasedPrecis (bool incp)
      {
        incPrec = incp;
      }

      /*! \brief Toggle the orientation of generated random numbers.
       *
       * This methods causes the given stream to generate antithetic
       * (1 - U, where U is the default number generated) when called
       * with \a a = true.  Calling this method with \a a = false will
       * return generated numbers to their default orientation.
       *
       * \param a A boolean value that selects regular or antithetic
       * variates.
       *
       * \see IncreasedPrecis(bool)
       */
      void
      SetAntithetic (bool a)
      {
        anti = a;
      }

      /*! \brief Generate a random uniform variate on (0, 1).
       *
       * This routine returns a random double precision floating point
       * number from the uniform distribution on the interval (0,
       * 1).  This method overloads the pure virtual method of the
       * same name in the rng base class.
       *
       * \see runif(unsigned int, unsigned int)
       * \see RandInt(long, long)
       * \see rng
       */
      double
      runif ()
      {
        if (incPrec)
          return U01d();
        else
          return U01();
      }

      /* We have to override the overloaded form of runif because
       * overloading the no-arg runif() hides the base class
       * definition; C++ stops looking once it finds the above.
       */
      /*! \brief Generate a Matrix of random uniform variates.
       *
       * This routine returns a Matrix of double precision random
       * uniform variates. on the interval (0, 1).  This method
       * overloads the virtual method of the same name in the rng base
       * class.
       *
       * This is the general template version of this method and
       * is called through explicit template instantiation.
       *
       * \param rows The number of rows in the returned Matrix.
       * \param cols The number of columns in the returned Matrix.
       * 
       * \see runif()
       * \see rng
       *
       * \note We are forced to override this overloaded method
       * because the 1-arg version of runif() hides the base class's
       * definition of this method from the compiler, although it
       * probably should not.
       */
      template <matrix_order O, matrix_style S>
      Matrix<double,O,S> runif(unsigned int rows, unsigned int cols)
      {
        return rng<lecuyer>::runif<O,S>(rows,cols);
      }

      /*! \brief Generate a Matrix of random uniform variates.
       *
       * This routine returns a Matrix of double precision random
       * uniform variates on the interval (0, 1).  This method
       * overloads the virtual method of the same name in the rng base
       * class.
       *
       * This is the default template version of this method and
       * is called through implicit template instantiation.
       *
       * \param rows The number of rows in the returned Matrix.
       * \param cols The number of columns in the returned Matrix.
       * 
       * \see runif()
       * \see rng
       *
       * \note We are forced to override this overloaded method
       * because the 1-arg version of runif() hides the base class's
       * definition of this method from the compiler, although it
       * probably should not.
       */
      Matrix<double,Col,Concrete> runif(unsigned int rows,
                                        unsigned int cols)
      {
        return rng<lecuyer>::runif<Col,Concrete>(rows, cols);
      }

      /*! \brief Generate the next random integer.
       *
       * This method generates a random integer from the discrete
       * uniform distribution on the interval [\a low, \a high].
       *
       * \param low The lower bound of the interval to evaluate.
       * \param high the upper bound of the interval to evaluate.
       *
       * \see runif()
       */
      long
      RandInt (long low, long high)
      {
        return low + static_cast<long> ((high - low + 1) * runif ());
      }

    protected:
      // Generate the next random number.
      //
      double
      U01 ()
      {
        long k;
        double p1, p2, u;

        /* Component 1 */
        p1 = a12 * Cg[1] - a13n * Cg[0];
        k = static_cast<long> (p1 / m1);
        p1 -= k * m1;
        if (p1 < 0.0) p1 += m1;
        Cg[0] = Cg[1]; Cg[1] = Cg[2]; Cg[2] = p1;

        /* Component 2 */
        p2 = a21 * Cg[5] - a23n * Cg[3];
        k = static_cast<long> (p2 / m2);
        p2 -= k * m2;
        if (p2 < 0.0) p2 += m2;
        Cg[3] = Cg[4]; Cg[4] = Cg[5]; Cg[5] = p2;

        /* Combination */
        u = ((p1 > p2) ? (p1 - p2) * norm : (p1 - p2 + m1) * norm);

        return (anti == false) ? u : (1 - u);
      }

      // Generate the next random number with extended (53 bits) precision.
      double 
      U01d ()
      {
        double u;
        u = U01();
        if (anti) {
          // Don't forget that U01() returns 1 - u in the antithetic case
          u += (U01() - 1.0) * fact;
          return (u < 0.0) ? u + 1.0 : u;
        } else {
          u += U01() * fact;
          return (u < 1.0) ? u : (u - 1.0);
        }
      }


      // Public members of the class start here
      
      // The default seed of the package; will be the seed of the first
      // declared RngStream, unless SetPackageSeed is called.
      static double nextSeed[6];

      /* Instance variables */
      double Cg[6], Bg[6], Ig[6];


      bool anti, incPrec;


      std::string streamname_;

  };

  /* Default seed definition */
  double lecuyer::nextSeed[6] = 
      {
         12345.0, 12345.0, 12345.0, 12345.0, 12345.0, 12345.0
      };

}

#endif /* SCYTHE_LECUYER_H */
