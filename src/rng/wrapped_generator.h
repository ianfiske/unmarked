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
 * scythestat/rng/wrapped_generator.h
 *
 * Provides a class definition that allows users to adapt non-Scythe
 * pseudo-random number generators to Scythe's rng interface.
 * Specifically, wraps any functor that generators uniform variates on
 * (0, 1).
 *
 */

 /*! \file wrapped_generator.h
  * \brief Adaptor for non-Scythe quasi-random number generators.
  *
  * This file contains the wrapped_generator class, a class that 
  * extends Scythe's base random number generation class (scythe::rng)
  * by allowing an arbitrary random uniform number generator to act as
  * the engine for random number generation in Scythe.
  */

#ifndef SCYTHE_WRAPPED_GENERATOR_H
#define SCYTHE_WRAPPED_GENERATOR_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "rng.h"
#else
#include "scythestat/rng.h"
#endif

namespace scythe {

   /*! \brief Adaptor for non-Scythe quasi-random number generators.
    *
    * This class defines a wrapper for arbitrary random uniform number
    * generators, allowing them to act as the underlying engine for
    * random number generation in Scythe.  Specifically, any function
    * object that overloads the function call operator to return
    * random uniform deviates on the interval (0, 1).
    *
    * The wrapped_generator class extends Scythe's basic random number
    * generating class, scythe::rng, implementing the interface that it
    * defines.
    *
    * \see rng
    * \see lecuyer
    * 
    */
  template <typename ENGINE>
  class wrapped_generator: public rng<wrapped_generator<ENGINE> >
	{
		public:

      /*! \brief Default constructor
       *
       * This constructor wraps the provided random uniform number
       * generating function object, creating an object suitable for
       * random number generation in Scythe.  Note that the function
       * object is passed by reference and is not copied on
       * construction.
       *
       * \param e A function object that returns uniform random
       * numbers on (0,1) when invoked.
       *
       * \see wrapped_generator(const wrapped_generator& wg)
       */
      wrapped_generator (ENGINE& e)
        : rng<wrapped_generator<ENGINE> > (),
          engine (e)
      {}

      /*! \brief Copy constructor
       *
       * This constructor makes a copy of an existing
       * wrapped_generator object, duplicating its seed and current
       * state exactly.  Note that this will create a copy of the
       * underlying function object using the function objects copy
       * construction semantics.
       *
       * \param wg An existing wrapped_generator object.
       *
       * \see wrapped_generator(ENGINE& e)
       */
      wrapped_generator(const wrapped_generator& wg)
        : rng<wrapped_generator<ENGINE> > (),
          engine (wg.engine)
      {}

      /*! \brief Generate a random uniform variate on (0, 1).
       *
       * This routine returns a random double precision floating point
       * number from the uniform distribution on the interval (0,
       * 1).  This method overloads the pure virtual method of the
       * same name in the rng base class.
       *
       * \see runif(unsigned int, unsigned int)
       * \see rng
       */
      inline double runif()
      {
        return engine();
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
      inline Matrix<double,O> runif(unsigned int rows,
                                    unsigned int cols)
      {
        return rng<wrapped_generator<ENGINE> >::runif<O, S>(rows, cols);
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
      Matrix<double,Col,Concrete> runif (unsigned int rows,
                                         unsigned int cols)
      {
        return rng<wrapped_generator<ENGINE> >::runif<Col,Concrete>(rows, 
            cols);
      }

		protected:
      ENGINE& engine; // The wrapped runif engine
	};
} // end namespace scythe

#endif /* SCYTHE_WRAPPED_GENERATOR_H */
