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
 * scythestat/defs.h
 */

/*!  \file defs.h 
 * \brief Global Scythe definitions.
 *
 * This file provides a variety of global definitions used throughout
 * the Scythe library.
 *
 * The majority of these definitions are used only within the library
 * itself.  Those definitions that are part of the public interface
 * are documented.
 *
 */

/* Doxygen main page text */
/*! \mainpage Scythe Statistical Library: Application Programmers' Interface
 *
 * \section intro Introduction
 *
 * The Scythe Statistical Library is an open source C++ library for
 * statistical computation, written by Daniel Pemstein (University of
 * Illinois), Kevin M. Quinn (Harvard University), and Andrew D.
 * Martin (Washington University). It includes a suite of matrix
 * manipulation functions, a suite of pseudo-random number generators,
 * and a suite of numerical optimization routines. Programs written
 * using Scythe are generally much faster than those written in
 * commonly used interpreted languages, such as R and MATLAB, and can
 * be compiled on any system with the GNU GCC compiler (and perhaps
 * with other C++ compilers). One of the primary design goals of the
 * Scythe developers has been ease of use for non-expert C++
 * programmers. We provide ease of use through three primary
 * mechanisms: (1) operator and function over-loading, (2) numerous
 * pre-fabricated utility functions, and (3) clear documentation and
 * example programs. Additionally, Scythe is quite flexible and
 * entirely extensible because the source code is available to all
 * users under the GNU General Public License.
 *
 * \section thisdoc About This Document
 *
 * This document is the application programmer's interface (API) to
 * Scythe.  It provides documentation for every class, function, and
 * object in Scythe that is part of the library's public interface.
 * In addition, the sections below explain how to obtain, install, and
 * compile the library.
 *
 * \section obtain Obtaining Scythe
 *
 * The most recent version of Scythe is available for download at
 * http://scythe.wustl.edu.
 *
 * \section install Installation
 * 
 * Scythe installs as a header-only C++ library.  After uncompressing,
 * simply follow the instructions in the INSTALL file included with
 * Scythe to install the library.  Alternatively, you may copy the
 * source files in scythestat and scythestat/rng into your project
 * directory and compile directly, using the SCYTHE_COMPILE_DIRECT
 * pre-processor flag.
 *
 * \section compile Compilation
 *
 * Scythe should work with the GNU GCC compiler, version 4.0.1 and
 * greater.  Scythe has not been tested with other compilers.  Scythe
 * provides a number of pre-processor flags.  The
 * SCYTHE_COMPILE_DIRECT allows the user to compile Scythe sources
 * directly.  The SCYTHE_VIEW_ASSIGNMENT_FLAG turns on R-style
 * recycling in Matrix::operator=() for view matrices. 
 *
 * The SCYTHE_DEBUG_LEVEL controls the amount of error trapping in
 * Scythe.  This level ranges from 0 (virtually no checking) to 3 (all
 * checking, including Matrix bounds checking, turned on).  By
 * default, the level is set to 3.  Here's an example of how to
 * compile a program with only basic error checking:
 *
 * \verbatim $ g++ myprog.cc -DSCYTHE_DEBUG_LEVEL=1 \endverbatim
 *
 * Finally, the SCYTHE_LAPACK flag enables LAPACK/BLAS.  You must have
 * the LAPACK and BLAS libraries installed on your system and compile
 * your program with the appropriate linker flags for this to work.  For
 * example, on linux you can enable LAPACK/BLAS support like this:
 *
 * \verbatim $ g++ myprog.cc -DSCYTHE_LAPACK -llapack -lblas -pthread \endverbatim
 * 
 * \section copy Copyright
 *
 * Scythe Statistical Library Copyright (C) 2000-2002 Andrew D. Martin
 * and Kevin M.  Quinn; 2002-present Andrew D. Martin, Kevin M. Quinn,
 * and Daniel Pemstein.  All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify under the terms of the GNU General Public License as
 * published by Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.  See the text files
 * COPYING and LICENSE, distributed with library's source code, for
 * further information.
 *
 * \section acknowledge Acknowledgments
 *
 * We gratefully acknowledge support from the United States National
 * Science Foundation (Grants SES-0350646 and SES-0350613), the
 * Department of Political Science, the Weidenbaum Center, and the
 * Center for Empirical Research in the Law at Washington University,
 * and the Department of Government and The Institute for Quantitative
 * Social Science at Harvard University. Neither the foundation,
 * Washington University, nor Harvard University bear any
 * responsibility for this software.
 *
 * We'd also like to thank the research assistants who have helped us
 * with Scythe: Matthew Fasman, Steve Haptonstahl, Kate Jensen, Laura
 * Keys, Kevin Rompala, Joe Sheehan, and Jean Yang.
 */

#ifndef SCYTHE_DEFS_H
#define SCYTHE_DEFS_H

/* In many functions returning matrices, we want to allow the user to
 * get a matrix of any style, but want to work with concretes inside
 * the function, for efficiency.  This macro originally contained the
 * code:
 * 
 * if (_STYLE_ == View)                                                \
 *   return Matrix<_TYPE_,_ORDER_,View>(_MATRIX_);                     \
 * else                                                                \
 *   return _MATRIX_;
 *
 * to convert to View before return if necessary.  Of course, this is
 * completely redundant, since the copy constructor gets called on
 * return anyway, so the body of the macro was replaced with a simple
 * return.  If we change our minds down the road about how to handle
 * these returns, code changes will be centered on this macro.
 */
#define SCYTHE_VIEW_RETURN(_TYPE_, _ORDER_, _STYLE_, _MATRIX_)        \
    return _MATRIX_;

/* Some macros to do bounds checking for iterator accesses.  The first
 * two are only called by the [] operator in the random access
 * iterator.  The third macro handles everything for checks on simple
 * current iterator location accesses.
 */
#define SCYTHE_ITER_CHECK_POINTER_BOUNDS(POINTER)                     \
{                                                                     \
	SCYTHE_CHECK_30(POINTER >= start_ + size_ || POINTER < start_,      \
		scythe_bounds_error, "Iterator access (offset "                   \
		<< offset_ << ") out of matrix bounds")                           \
}

#define SCYTHE_ITER_CHECK_OFFSET_BOUNDS(OFFSET)                       \
{                                                                     \
	SCYTHE_CHECK_30(OFFSET >= size_, scythe_bounds_error,      	        \
		"Iterator access (offset " << offset_ << ") out of matrix bounds")\
}

#define SCYTHE_ITER_CHECK_BOUNDS()                                    \
{                                                                     \
	if (M_STYLE != Concrete || M_ORDER != ORDER) {                      \
		SCYTHE_ITER_CHECK_OFFSET_BOUNDS(offset_);                         \
  } else {                                                            \
		SCYTHE_ITER_CHECK_POINTER_BOUNDS(pos_);													  \
	}                                                                   \
}

/*! \namespace scythe
 * \brief The Scythe library namespace.
 *
 * All Scythe library declarations are defined within the scythe
 * namespace.  This prevents name clashing with other libraries'
 * members or with declarations in users' program code.
 */
namespace scythe {

  /*! 
   * \brief Matrix order enumerator.
   *
   * Matrix templates may be either column-major or row-major ordered
   * and this enumerator is used to differentiate between the two
   * types.
   *
   * The enumerator provides two values: Concrete and View.
   *
   * \see Matrix
   */
  enum matrix_order { Col, Row };

  /*! 
   * \brief Matrix style enumerator.
   *
   * Matrix templates may be either concrete matrices or views and
   * this enumerator is used to differentiate between the two types.
   *
   * Concrete matrices provide direct access to an underlying array of
   * matrix data, while views offer a more general interface to data
   * arrays, with potentially many views referencing the same
   * underlying data structure.
   *
   * The enum provides two values: Col and Row.
   *
   * \see Matrix
   */
  enum matrix_style { Concrete, View };

  /*!
   * \brief A convenient marker for vector submatrix access.
   
   * Passing an all_elements object to a two-argument Matrix submatrix
   * method allows the caller to access a full vector submatrix.  We
   * further define an instance of all_elements named "_" in the
   * scythe namespace to allow users to easily reference entire
   * vectors within matrices.
   *
   * \see Matrix::operator()(const all_elements, uint)
   * \see Matrix::operator()(const all_elements, uint) const
   * \see Matrix::operator()(uint, const all_elements)
   * \see Matrix::operator()(uint, const all_elements) const
   *
   */

	struct all_elements {
	} const _ = {};

  // A little helper method to see if more col-order or row-order.
  // Tie breaks towards col.
  template <matrix_order o1, matrix_order o2, matrix_order o3>
  bool maj_col()
  {
    if ((o1 == Col && o2 == Col) ||
        (o1 == Col && o3 == Col) ||
        (o2 == Col && o3 == Col))
      return true;
    return false;
  }

  template <matrix_order o1, matrix_order o2, matrix_order o3,
            matrix_order o4>
  bool maj_col()
  {
    if ((o1 == Col && o2 == Col) ||
        (o1 == Col && o3 == Col) ||
        (o1 == Col && o4 == Col) ||
        (o2 == Col && o3 == Col) ||
        (o2 == Col && o4 == Col) ||
        (o3 == Col && o4 == Col))
      return true;
    return false;
  }
}  // end namespace scythe

#endif /* SCYTHE_ERROR_H */
