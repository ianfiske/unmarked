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
 *  scythestat/matrix_random_access_iterator.h
 *
 * Random access iterators for the matrix class.
 *
 */

/*! \file matrix_random_access_iterator.h
 * \brief Definitions of STL-compliant random access iterators for
 * the Matrix class.
 *
 * Contains definitions of const_matrix_random_access_iterator,
 * matrix_random_access_iterator, and related operators.  See a
 * Standard Template Library reference, such as Josuttis (1999), for a
 * full description of the capabilities of random access iterators.
 *
 * These iterators are templated on the type, order and style of the
 * Matrix they iterate over and their own order, which need not match
 * the iterated-over matrix.  Same-order iteration over concrete
 * matrices is extremely fast.  Cross-grain concrete and/or view
 * iteration is slower.  
 */

#ifndef SCYTHE_MATRIX_RANDOM_ACCESS_ITERATOR_H
#define SCYTHE_MATRIX_RANDOM_ACCESS_ITERATOR_H

#include <iterator>

#ifdef SCYTHE_COMPILE_DIRECT
#include "defs.h"
#include "error.h"
#include "matrix.h"
#else
#include "scythestat/defs.h"
#include "scythestat/error.h"
#include "scythestat/matrix.h"
#endif

/* The const_matrix_iterator and matrix_iterator classes are
 * essentially identical, except for the return types of the *, ->,
 * and [] operators.  matrix_iterator extends const_matrix_iterator,
 * overriding most of its members. */

/* TODO Current setup uses template argument based branches to
 * handle views and cross-grained orderings differently than simple
 * in-order concrete matrices.  The work for this gets done at
 * compile time, but we end with a few unused instance variables in
 * the concrete case.  It might be better to specialize the entire
 * class, although this will lead to a lot of code duplication.  We
 * should bench the difference down the road and see if it is worth
 * the maintenance hassle.
 *
 * At the moment this is looking like it won't be worth it.
 * Iterator-based operations on concretes provide comparable
 * performance to element-access based routines in previous versions
 * of the library, indicating little performance penalty.
 */

namespace scythe {
	/* convenience typedefs */
  namespace { // local to this file
    typedef unsigned int uint;
  }

	/* forward declaration of the matrix class */
	template <typename T_type, matrix_order ORDER, matrix_style STYLE>
	class Matrix;

  /*! \brief An STL-compliant const random access iterator for Matrix.
   *
   * Provides random access iteration over const Matrix objects.  See
   * Josuttis (1999), or some other STL reference, for a full
   * description of the random access iterator interface.
   *
   * \see Matrix
   * \see matrix_random_access_iterator
   * \see const_matrix_forward_iterator
   * \see matrix_forward_iterator
   * \see const_matrix_bidirectional_iterator
   * \see matrix_bidirectional_iterator
   */

  template <typename T_type, matrix_order ORDER, matrix_order M_ORDER,
            matrix_style M_STYLE>
  class const_matrix_random_access_iterator
    : public std::iterator<std::random_access_iterator_tag, T_type>
  {
		public:
			/**** TYPEDEFS ***/
			typedef const_matrix_random_access_iterator<T_type, ORDER, 
              M_ORDER, M_STYLE> self;

			/* These are a little formal, but useful */
			typedef typename std::iterator_traits<self>::value_type
				value_type;
			typedef typename std::iterator_traits<self>::iterator_category
				iterator_category;
			typedef typename std::iterator_traits<self>::difference_type
				difference_type;
			typedef typename std::iterator_traits<self>::pointer pointer;
			typedef typename std::iterator_traits<self>::reference reference;

		
			/**** CONSTRUCTORS ****/
			
			/* Default constructor */
			const_matrix_random_access_iterator ()
			{}

			/* Standard constructor */
			const_matrix_random_access_iterator
        ( const Matrix<value_type, M_ORDER, M_STYLE> &M)
        : start_ (M.getArray())
      {
        SCYTHE_CHECK_30 (start_ == 0, scythe_null_error,
            "Requesting iterator to NULL matrix");
        pos_ = start_;

        /* The basic story is: when M_STYLE == Concrete and ORDER ==
         * M_ORDER, we only need pos_ and start_ and iteration will be
         * as fast as possible.  All other types of iteration need
         * more variables to keep track of things and are slower.
         */

        if (M_STYLE != Concrete || M_ORDER != ORDER) {
          offset_ = 0;

          if (ORDER == Col) {
            lead_length_ = M.rows();
            lead_inc_ = M.rowstride();
            trail_inc_ = M.colstride();
          } else {
            lead_length_ = M.cols();
            lead_inc_ = M.colstride();
            trail_inc_ = M.rowstride();
          }
          jump_ = trail_inc_ + (1 - lead_length_) * lead_inc_;
        }

#if SCYTHE_DEBUG > 2
				size_ = M.size();
#endif
      }

      /* Copy constructor */
      const_matrix_random_access_iterator (const self &mi)
        : start_ (mi.start_),
          pos_ (mi.pos_)
      {
        if (M_STYLE != Concrete || M_ORDER != ORDER) {
          offset_ = mi.offset_;
          lead_length_ = mi.lead_length_;
          lead_inc_ = mi.lead_inc_;
          trail_inc_ = mi.trail_inc_;
          jump_ = mi.jump_;
        }
#if SCYTHE_DEBUG > 2
				size_ = mi.size_;
#endif
      }

      /**** FORWARD ITERATOR FACILITIES ****/

      inline self& operator= (const self& mi)
      {
        start_ = mi.start_;
        pos_ = mi.pos_;

        if (M_STYLE != Concrete || M_ORDER != ORDER) {
          offset_ = mi.offset_;
          lead_length_ = mi.lead_length_;
          lead_inc_ = mi.lead_inc_;
          trail_inc_ = mi.trail_inc_;
          jump_ = mi.jump_;
        }
#if SCYTHE_DEBUG > 2
				size_ = mi.size_;
#endif

        return *this;
      }

      inline const reference operator* () const
      {
				SCYTHE_ITER_CHECK_BOUNDS();
        return *pos_;
      }

      inline const pointer operator-> () const
      {
				SCYTHE_ITER_CHECK_BOUNDS();
        return pos_;
      }

      inline self& operator++ ()
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER)
          ++pos_;
        else if (++offset_ % lead_length_ == 0)
          pos_ += jump_;
        else
          pos_ += lead_inc_;

        return *this;
      }

      inline self operator++ (int)
      {
        self tmp = *this;
        ++(*this);
        return tmp;
      }

      /* == is only defined for iterators of the same template type
       * that point to the same matrix.  Behavior for any other
       * comparison is undefined.
       *
       * Note that we have to be careful about iterator comparisons
       * when working with views and cross-grain iterators.
       * Specifically, we always have to rely on the offset value.
       * Obviously, with <> checks pos_ can jump all over the place in
       * cross-grain iterators, but also end iterators point to the
       * value after the last in the matrix.  In some cases, the
       * equation in += and -= will actually put pos_ inside the
       * matrix (often in an early position) in this case.
       */
      inline bool operator== (const self& x) const
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
          return pos_ == x.pos_;
        } else {
          return offset_ == x.offset_;
        }
      }

      /* Again, != is only officially defined for iterators over the
       * same matrix although the test will be trivially true for
       * matrices that don't view the same data, by implementation.
       */
      inline bool operator!= (const self &x) const
      {
        return !(*this == x);
      }

      /**** BIDIRECTIONAL ITERATOR FACILITIES ****/
        
      inline self& operator-- ()
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER)
          --pos_;
        else if (offset_-- % lead_length_ == 0)
          pos_ -= jump_;
        else
          pos_ -= lead_inc_;

        return *this;
      }

      inline self operator-- (int)
      {
        self tmp = *this;
        --(*this);
        return tmp;
      }

      /**** RANDOM ACCESS ITERATOR FACILITIES ****/

      inline const reference operator[] (difference_type n) const
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
					SCYTHE_ITER_CHECK_OFFSET_BOUNDS(start_ + n);
          return *(start_ + n);
        } else {
          uint trailing = n / lead_length_;
          uint leading = n % lead_length_;

          T_type* place = start_ + leading * lead_inc_
                                 + trailing * trail_inc_;

					SCYTHE_ITER_CHECK_POINTER_BOUNDS(place);
          return *place;
        }
      }

      inline self& operator+= (difference_type n)
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
          pos_ += n;
        } else {
          offset_ += n;
          uint trailing = offset_ / lead_length_;
          uint leading = offset_ % lead_length_;

          pos_ = start_ + leading * lead_inc_ 
                        + trailing * trail_inc_;
        }

        return *this;
      }

      inline self& operator-= (difference_type n)
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
          pos_ -= n;
        } else {
          offset_ -= n;
          uint trailing = offset_ / lead_length_;
          uint leading = offset_ % lead_length_;

          pos_ = start_ + leading * lead_inc_ 
                        + trailing * trail_inc_;
        }

        return *this;
      }

      /* + and - difference operators are outside the class */

      inline difference_type operator- (const self& x) const
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
          return pos_ - x.pos_;
        } else {
          return offset_ - x.offset_;
        }
      }

      inline difference_type operator< (const self& x) const
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
          return pos_ < x.pos_;
        } else {
          return offset_ < x.offset_;
        }
      }

      inline difference_type operator> (const self& x) const
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
          return pos_ > x.pos_;
        } else {
          return offset_ > x.offset_;
        }
      }

      inline difference_type operator<= (const self& x) const
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
          return pos_ <= x.pos_;
        } else {
          return offset_ <= x.offset_;
        }
      }

      inline difference_type operator>= (const self& x) const
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
          return pos_ >= x.pos_;
        } else {
          return offset_ >= x.offset_;
        }
      }

    protected:

      /**** INSTANCE VARIABLES ****/
      T_type* start_; // pointer to beginning of data array
      T_type* pos_;         // pointer to current position in array

      uint offset_;  // Logical offset into matrix

      // TODO Some of these can probably be uints
      int lead_length_;  // Logical length of leading dimension
      int lead_inc_;  // Memory distance between vectors in ldim
      int trail_inc_; // Memory distance between vectors in tdim
      int jump_; // Memory distance between end of one ldim vector and
                 // begin of next

			// Size variable for range checking
#if SCYTHE_DEBUG > 2
			uint size_;  // Logical matrix size
#endif
      
 };

  /*! \brief An STL-compliant random access iterator for Matrix.
   *
   * Provides random access iteration over Matrix objects.  See
   * Josuttis (1999), or some other STL reference, for a full
   * description of the random access iterator interface.
   *
   * \see Matrix
   * \see const_matrix_random_access_iterator
   * \see const_matrix_forward_iterator
   * \see matrix_forward_iterator
   * \see const_matrix_bidirectional_iterator
   * \see matrix_bidirectional_iterator
   */
	template <typename T_type, matrix_order ORDER, matrix_order M_ORDER,
            matrix_style M_STYLE>
	class matrix_random_access_iterator
		: public const_matrix_random_access_iterator<T_type, ORDER, 
                                                 M_ORDER, M_STYLE>
	{
			/**** TYPEDEFS ***/
			typedef matrix_random_access_iterator<T_type, ORDER, M_ORDER, 
                                            M_STYLE> self;
			typedef const_matrix_random_access_iterator<T_type, ORDER, 
                                                 M_ORDER, M_STYLE> Base;
		
		public:
			/* These are a little formal, but useful */
			typedef typename std::iterator_traits<Base>::value_type
				value_type;
			typedef typename std::iterator_traits<Base>::iterator_category
				iterator_category;
			typedef typename std::iterator_traits<Base>::difference_type
				difference_type;
			typedef typename std::iterator_traits<Base>::pointer pointer;
			typedef typename std::iterator_traits<Base>::reference reference;

		
			/**** CONSTRUCTORS ****/
			
			/* Default constructor */
			matrix_random_access_iterator ()
				: Base () 
			{}

			/* Standard constructor */
			matrix_random_access_iterator (const Matrix<value_type, M_ORDER, 
                                                  M_STYLE> &M)
				:	Base(M)
			{}

      /* Copy constructor */
			matrix_random_access_iterator (const self &mi)
				:	Base (mi)
			{}

			/**** FORWARD ITERATOR FACILITIES ****/

			/* We have to override a lot of these to get return values
			 * right.*/
      inline self& operator= (const self& mi)
      {
        start_ = mi.start_;
        pos_ = mi.pos_;

        if (M_STYLE != Concrete || M_ORDER != ORDER) {
          offset_ = mi.offset_;
          lead_length_ = mi.lead_length_;
          lead_inc_ = mi.lead_inc_;
          trail_inc_ = mi.trail_inc_;
          jump_ = mi.jump_;
        }
#if SCYTHE_DEBUG > 2
				size_ = mi.size_;
#endif

        return *this;
      }

			inline reference operator* () const
			{
				SCYTHE_ITER_CHECK_BOUNDS();
				return *pos_;
			}

			inline pointer operator-> () const
			{
				SCYTHE_ITER_CHECK_BOUNDS();
				return pos_;
			}

			inline self& operator++ ()
			{
				Base::operator++();
				return *this;
			}

			inline self operator++ (int)
			{
				self tmp = *this;
				++(*this);
				return tmp;
			}

			/**** BIDIRECTIONAL ITERATOR FACILITIES ****/

			inline self& operator-- ()
			{
				Base::operator--();
				return *this;
			}

			inline self operator-- (int)
			{
				self tmp = *this;
				--(*this);
				return tmp;
			}

			/**** RANDOM ACCESS ITERATOR FACILITIES ****/

      inline reference operator[] (difference_type n) const
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
					SCYTHE_ITER_CHECK_POINTER_BOUNDS(start_ + n);
          return *(start_ + n);
        } else {
          uint trailing = n / lead_length_;
          uint leading = n % lead_length_;

          T_type* place = start_ + leading * lead_inc_
                                 + trailing * trail_inc_;

					SCYTHE_ITER_CHECK_POINTER_BOUNDS(place);
          return *place;
        }
      }

			inline self& operator+= (difference_type n)
			{
				Base::operator+=(n);
				return *this;
			}

			inline self& operator-= (difference_type n)
			{
				Base::operator-= (n);
				return *this;
			}

			/* +  and - difference_type operators are outside the class */

		private:
			/* Get handles to base members.  It boggles the mind */
			using Base::start_;
			using Base::pos_;
      using Base::offset_;
      using Base::lead_length_;
      using Base::lead_inc_;
      using Base::trail_inc_;
      using Base::jump_;
#if SCYTHE_DEBUG > 2
			using Base::size_;
#endif
	};

	template <class T_type, matrix_order ORDER, matrix_order M_ORDER,
            matrix_style STYLE>
	inline
  const_matrix_random_access_iterator<T_type, ORDER, M_ORDER, STYLE>
	operator+ (const_matrix_random_access_iterator<T_type, ORDER, M_ORDER,                                                 STYLE> x, int n)
	{
		x += n;
		return x;
	}
	
	template <class T_type, matrix_order ORDER, matrix_order M_ORDER,
            matrix_style STYLE>
	inline
  const_matrix_random_access_iterator<T_type, ORDER, M_ORDER, STYLE>
	operator+ (int n, const_matrix_random_access_iterator<T_type, ORDER,
                                                        M_ORDER, 
                                                        STYLE> x)
	{
		x += n;
		return x;
	}

	template <class T_type, matrix_order ORDER, matrix_order M_ORDER,
            matrix_style STYLE>
	inline
  const_matrix_random_access_iterator<T_type, ORDER, M_ORDER, STYLE>
	operator- (const_matrix_random_access_iterator<T_type, ORDER, M_ORDER,                                                 STYLE> x, int n)
	{
		x -= n;
		return x;
	}

	template <class T_type, matrix_order ORDER, matrix_order M_ORDER,
            matrix_style STYLE>
	inline matrix_random_access_iterator<T_type, ORDER, M_ORDER, STYLE>
	operator+ (matrix_random_access_iterator<T_type, ORDER, M_ORDER,
                                           STYLE> x, int n)
	{
		x += n;
		return x;
	}
	
	template <class T_type, matrix_order ORDER, matrix_order M_ORDER,
            matrix_style STYLE>
	inline matrix_random_access_iterator<T_type, ORDER, M_ORDER, STYLE>
	operator+ (int n, matrix_random_access_iterator<T_type, ORDER,
                                                  M_ORDER, STYLE> x)
	{
		x += n;
		return x;
	}

	template <class T_type, matrix_order ORDER, matrix_order M_ORDER,
            matrix_style STYLE>
	inline matrix_random_access_iterator<T_type, ORDER, M_ORDER, STYLE>
	operator- (matrix_random_access_iterator<T_type, ORDER, M_ORDER,
                                           STYLE> x, int n)
	{
		x -= n;
		return x;
	}

} // namespace scythe

#endif /* SCYTHE_MATRIX_ITERATOR_H */
