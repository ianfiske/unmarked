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
 *  scythestat/matrix_forward_iterator.h
 *
 * Forward iterators for the matrix class.
 *
 */

/*! \file matrix_forward_iterator.h
 * \brief Definitions of STL-compliant forward iterators for the
 * Matrix class.
 *
 * Contains definitions of const_matrix_forward_iterator,
 * matrix_forward_iterator, and related operators.  See a Standard
 * Template Library reference, such as Josuttis (1999), for a full
 * description of the capabilities of forward iterators.
 *
 * These iterators are templated on the type, order and style of the
 * Matrix they iterate over and their own order, which need not match
 * the iterated-over matrix.  Same-order iteration over concrete
 * matrices is extremely fast.  Cross-grain concrete and/or view
 * iteration is slower.  
 */

#ifndef SCYTHE_MATRIX_FORWARD_ITERATOR_H
#define SCYTHE_MATRIX_FORWARD_ITERATOR_H

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

namespace scythe {
	/* convenience typedefs */
  namespace { // local to this file
	  typedef unsigned int uint;
  }
  
	/* forward declaration of the matrix class */
	template <typename T_type, matrix_order ORDER, matrix_style STYLE>
	class Matrix;

  /*! \brief An STL-compliant const forward iterator for Matrix.
   *
   * Provides forward iteration over const Matrix objects.  See
   * Josuttis (1999), or some other STL reference, for a full
   * description of the forward iterator interface.
   *
   * \see Matrix
   * \see matrix_forward_iterator
   * \see const_matrix_random_access_iterator
   * \see matrix_random_access_iterator
   * \see const_matrix_bidirectional_iterator
   * \see matrix_bidirectional_iterator
   */
  
  template <typename T_type, matrix_order ORDER, matrix_order M_ORDER,
            matrix_style M_STYLE>
  class const_matrix_forward_iterator
    : public std::iterator<std::forward_iterator_tag, T_type>
  {
		public:
			/**** TYPEDEFS ***/
			typedef const_matrix_forward_iterator<T_type, ORDER, 
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
			const_matrix_forward_iterator ()
			{}

			/* Standard constructor */
			const_matrix_forward_iterator
        (const Matrix<value_type, M_ORDER, M_STYLE> &M)
        : pos_ (M.getArray()),
          matrix_ (&M)
      {
        SCYTHE_CHECK_30 (pos_ == 0, scythe_null_error,
            "Requesting iterator to NULL matrix");

        /* The basic story is: when M_STYLE == Concrete and ORDER ==
         * M_ORDER, we only need pos_ and iteration will be as fast as
         * possible.  All other types of iteration need more variables
         * to keep track of things and are slower.
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
          vend_ = pos_ + (lead_length_ - 1) * lead_inc_;
        }
#if SCYTHE_DEBUG > 2
				size_ = M.size();
        start_ = pos_;
#endif
      }

      /* Copy constructor */
      const_matrix_forward_iterator (const self &mi)
        : pos_ (mi.pos_),
          matrix_ (mi.matrix_)
      {
        if (M_STYLE != Concrete || M_ORDER != ORDER) {
          offset_ = mi.offset_;
          lead_length_ = mi.lead_length_;
          lead_inc_ = mi.lead_inc_;
          trail_inc_ = mi.trail_inc_;
          vend_ = mi.vend_;
          jump_ = mi.jump_;
        }
#if SCYTHE_DEBUG > 2
				size_ = mi.size_;
        start_ = mi.start_;
#endif
      }

      /**** EXTRA MODIFIER ****/

      /* This function lets us grab an end iterator quickly, for both
       * concrete and view matrices.  The view code is a bit of a
       * kludge, but it works.
       */
      inline self& set_end ()
      {
        if (M_STYLE == Concrete && ORDER == M_ORDER) {
          pos_ = matrix_->getArray() + matrix_->size();
        } else { 
          offset_ = matrix_->size();
        }

        return *this;
      }

      /* Returns the current index (in logical matrix terms) of the
       * iterator.
       */
      unsigned int get_index () const
      {
        return offset_;
      }

      /**** FORWARD ITERATOR FACILITIES ****/

      inline self& operator= (const self& mi)
      {
        pos_ = mi.pos_;
        matrix_ = mi.matrix_;

        if (M_STYLE != Concrete || M_ORDER != ORDER) {
          offset_ = mi.offset_;
          lead_length_ = mi.lead_length_;
          lead_inc_ = mi.lead_inc_;
          trail_inc_ = mi.trail_inc_;
          vend_ = mi.vend_;
          jump_ = mi.jump_;
        }
#if SCYTHE_DEBUG > 2
				size_ = mi.size_;
        start_ = mi.start_;
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
        else {
          if (pos_ == vend_) {
            vend_ += trail_inc_;
            pos_ += jump_;
          } else {
            pos_ += lead_inc_;
          }
          ++offset_;
        }

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

    protected:

      /**** INSTANCE VARIABLES ****/
      T_type* pos_;         // pointer to current position in array
      T_type *vend_;        // pointer to end of current vector

      uint offset_;         // logical offset into matrix

      // TODO Some of these can probably be uints
      int lead_length_;  // Logical length of leading dimension
      int lead_inc_;  // Memory distance between vectors in ldim
      int trail_inc_; // Memory distance between vectors in tdim
      int jump_; // Memory distance between end of one ldim vector and
                 // begin of next
      // Pointer to the matrix we're iterating over.  This is really
      // only needed to get variables necessary to set the end.
      // TODO Handle this more cleanly.
      const Matrix<T_type, M_ORDER, M_STYLE>* matrix_;
			// Size variable for range checking
#if SCYTHE_DEBUG > 2
			uint size_;  // Logical matrix size
      T_type* start_; // Not normally needed, but used for bound check
#endif
 };

  /*! \brief An STL-compliant forward iterator for Matrix.
   *
   * Provides forward iteration over Matrix objects.  See
   * Josuttis (1999), or some other STL reference, for a full
   * description of the forward iterator interface.
   *
   * \see Matrix
   * \see const_matrix_forward_iterator
   * \see const_matrix_random_access_iterator
   * \see matrix_random_access_iterator
   * \see const_matrix_bidirectional_iterator
   * \see matrix_bidirectional_iterator
   */
	template <typename T_type, matrix_order ORDER, matrix_order M_ORDER,
            matrix_style M_STYLE>
	class matrix_forward_iterator
		: public const_matrix_forward_iterator<T_type, ORDER, 
                                           M_ORDER, M_STYLE>
	{
			/**** TYPEDEFS ***/
			typedef matrix_forward_iterator<T_type, ORDER, M_ORDER, 
                                      M_STYLE> self;
			typedef const_matrix_forward_iterator<T_type, ORDER, 
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
			matrix_forward_iterator ()
				: Base () 
			{}

			/* Standard constructor */
			matrix_forward_iterator (const Matrix<value_type, M_ORDER, 
                                                  M_STYLE> &M)
				:	Base(M)
			{}

      /* Copy constructor */
			matrix_forward_iterator (const self &mi)
				:	Base (mi)
			{}

      /**** EXTRA MODIFIER ****/
      inline self& set_end ()
      {
        Base::set_end();
        return *this;
      }

			/**** FORWARD ITERATOR FACILITIES ****/

			/* We have to override a lot of these to get return values
			 * right.*/
      inline self& operator= (const self& mi)
      {
        pos_ = mi.pos_;
        matrix_ = mi.matrix_;

        if (M_STYLE != Concrete || M_ORDER != ORDER) {
          offset_ = mi.offset_;
          lead_length_ = mi.lead_length_;
          lead_inc_ = mi.lead_inc_;
          trail_inc_ = mi.trail_inc_;
          vend_ = mi.vend_;
          jump_ = mi.jump_;
        }
#if SCYTHE_DEBUG > 2
				size_ = mi.size_;
        start_ = mi.start_;
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

		private:
			/* Get handles to base members.  It boggles the mind */
			using Base::pos_;
      using Base::vend_;
      using Base::offset_;
      using Base::lead_length_;
      using Base::lead_inc_;
      using Base::trail_inc_;
      using Base::jump_;
      using Base::matrix_;
#if SCYTHE_DEBUG > 2
			using Base::size_;
      using Base::start_;
#endif
	};

} // namespace scythe

#endif /* SCYTHE_MATRIX_ITERATOR_H */
