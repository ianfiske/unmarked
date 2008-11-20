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
 * scythestat/rng/mersenne.h
 *
 * Provides the class definition for the mersenne random number
 * generator.  This class extends the base rng class by providing an
 * implementation of runif() based on an implementation of the
 * mersenne twister, released under the following license:
 *
 * A C-program for MT19937, with initialization improved 2002/1/26.
 * Coded by Takuji Nishimura and Makoto Matsumoto.
 * 
 * Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above
 *    copyright
 *    notice, this list of conditions and the following disclaimer
 *    in the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The names of its contributors may not be used to endorse or
 *    promote products derived from this software without specific
 *    prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * For more information see:
 * http://www.math.keio.ac.jp/matumoto/emt.html
 *
 */

 /*! \file mersenne.h
  * \brief The Mersenne Twister random number generator.  
  *
  * This file contains the mersenne class, a class that extends
  * Scythe's base random number generation class (scythe::rng) by
  * providing an implementation of scythe::rng::runif() using the
  * Mersenne Twister algorithm.
  */

#ifndef SCYTHE_MERSENNE_H
#define SCYTHE_MERSENNE_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "rng.h"
#else
#include "scythestat/rng.h"
#endif

namespace scythe {

#ifdef __MINGW32__
	/* constant vector a */
	static const unsigned long MATRIX_A = 0x9908b0dfUL;
	
	/* most significant w-r bits */
	static const unsigned long UPPER_MASK = 0x80000000UL;
	
	/* least significant r bits */
	static const unsigned long LOWER_MASK = 0x7fffffffUL;
#else
	namespace {
		/* constant vector a */
		const unsigned long MATRIX_A = 0x9908b0dfUL;
		
		/* most significant w-r bits */
		const unsigned long UPPER_MASK = 0x80000000UL;
		
		/* least significant r bits */
		const unsigned long LOWER_MASK = 0x7fffffffUL;
	}
#endif

   /*! \brief The Mersenne Twister random number generator.
    *
    * This class defines a random number generator, using the Mersenne
    * Twister algorithm developed and implemented by Makoto Matsumoto
    * and Takuji Nishimura (1997, 2002).  The period of this random
    * number generator is \f$2^{19937} - 1\f$.
    *
    * The mersenne class extends Scythe's basic random number
    * generating class, scythe::rng, implementing the interface that it
    * defines.
    *
    * \see rng
    * \see lecuyer
    * 
    */
	class mersenne: public rng<mersenne>
	{
		public:

      /*! \brief Default constructor
       *
       * This constructor generates an unseeded and uninitialized
       * mersenne object.  It is most useful for creating arrays of
       * random number generators.  An uninitialized mersenne object
       * will be seeded with the default seed (5489UL) automatically
       * upon use.
       *
       * \see mersenne(const mersenne &m)
       * \see initialize(unsigned long s)
       */
			mersenne ()
        :	rng<mersenne> (),
          mti (N + 1)
      {}

      /*! \brief Copy constructor
       *
       * This constructor makes a copy of an existing mersenne
       * object, duplicating its seed and current state exactly.
       *
       * \param m An existing mersenne random number generator.
       *
       * \see mersenne()
       */
			mersenne (const mersenne &m)
        : rng<mersenne> (),
          mti (m.mti)
      {
      }

      /*! \brief Sets the seed.
       *
       * This method sets the seed of the random number generator and
       * readies it to begin generating random numbers.  Calling this
       * function on a mersenne object that is already in use is
       * supported, although not suggested unless you know what you
       * are doing.
       *
       * \param s A long integer seed.
       *
       * \see mersenne()
       */
			void initialize (unsigned long s)
      {
        mt[0]= s & 0xffffffffUL;
        for (mti=1; mti<N; mti++) {
          mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30))
              + mti); 
          /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
          /* In the previous versions, MSBs of the seed affect   */
          /* only MSBs of the array mt[].                        */
          /* 2002/01/09 modified by Makoto Matsumoto             */
          mt[mti] &= 0xffffffffUL;
          /* for >32 bit machines */
        }
      }
			
      /*! \brief Generate a random uniform variate on (0, 1).
       *
       * This routine returns a random double precision floating point
       * number from the uniform distribution on the interval (0,
       * 1).  This method overloads the pure virtual method of the
       * same name in the rng base class.
       *
       * \see runif(unsigned int, unsigned int)
       * \see genrand_int32()
       * \see rng
       */
      inline double runif()
			{
				return (((double) genrand_int32()) + 0.5) * 
					(1.0 / 4294967296.0);
			}

			/* We have to override the overloaded forms of runif because
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
			inline Matrix<double,O,S> runif(unsigned int rows, 
                                      unsigned int cols)
			{
				return rng<mersenne>::runif<O,S>(rows, cols);
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
        return rng<mersenne>::runif<Col,Concrete>(rows, cols);
      }

      /* generates a random number on [0,0xffffffff]-interval */
      /*! \brief Generate a random long integer.
       *
       * This method generates a random integer, drawn from the
       * discrete uniform distribution on the interval [0,0xffffffff].
       *
       * \see runif()
       * \see initialize(unsigned long s)
       */
			unsigned long genrand_int32()
      {
        unsigned long y;
        static unsigned long mag01[2]={0x0UL, MATRIX_A};
        /* mag01[x] = x * MATRIX_A  for x=0,1 */

        if (mti >= N) { /* generate N words at one time */
          int kk;

          if (mti == N+1)   // if init_genrand() has not been called,
            this->initialize(5489UL); // a default initial seed is used

          for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
          }
          for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
          }
          y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
          mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

          mti = 0;
        }
      
        y = mt[mti++];

        /* Tempering */
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);

        return y;
      }
		
		protected:
			/* Period parameters */
			static const int N = 624;
			static const int M = 398;
		
			/* the array for the state vector  */
			unsigned long mt[N];

			/* mti==N+1 means mt[N] is not initialized */
			int mti;
	};

}

#endif /* SCYTHE_MERSENNE_H */
