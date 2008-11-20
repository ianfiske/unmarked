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
 *  scythe's/matrix.h
 *
 */

/*!
 * \file matrix.h
 * \brief Definitions of Matrix and related classes and functions.
 *
 * This file contains the definitions of the Matrix, Matrix_base, and
 * associated classes.  It also contains a number of external
 * functions that operate on Matrix objects, such as mathematical
 * operators.
 *
 * Many of the arithmetic and logical operators in this file are
 * implemented in terms of overloaded template definitions to provide
 * both generic and default versions of each operation.  Generic
 * templates allow (and force) the user to fully specify the 
 * template type of the returned matrix object (row or column order,
 * concrete or view) while default templates automatically return
 * concrete matrices with the ordering of the first or only Matrix
 * argument to the function.  Furthermore, we overload binary
 * functions to provide scalar by Matrix operations, in addition to
 * basic Matrix by Matrix arithmetic and logic.  Therefore,
 * definitions for multiple versions appear in the functions list
 * below.  We adopt the convention of providing explicit documentation
 * for only the most generic Matrix by Matrix version of each of these
 * operators and describe the behavior of the various overloaded
 * versions in these documents.
 */


#ifndef SCYTHE_MATRIX_H
#define SCYTHE_MATRIX_H

//#include <climits>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
//#include <numeric>
#include <functional>
#include <list>

#ifdef SCYTHE_COMPILE_DIRECT
#include "defs.h"
#include "algorithm.h"
#include "error.h"
#include "datablock.h"
#include "matrix_random_access_iterator.h"
#include "matrix_forward_iterator.h"
#include "matrix_bidirectional_iterator.h"
#ifdef SCYTHE_LAPACK
#include "lapack.h"
#endif
#else
#include "scythestat/defs.h"
#include "scythestat/algorithm.h"
#include "scythestat/error.h"
#include "scythestat/datablock.h"
#include "scythestat/matrix_random_access_iterator.h"
#include "scythestat/matrix_forward_iterator.h"
#include "scythestat/matrix_bidirectional_iterator.h"
#ifdef SCYTHE_LAPACK
#include "scythestat/lapack.h"
#endif
#endif

namespace scythe {

  namespace { // make the uint typedef local to this file
	  /* Convenience typedefs */
	  typedef unsigned int uint;
  }

  /* Forward declare the matrix multiplication functions for *= to use
   * within Matrix proper .
   */

  template <typename T_type, matrix_order ORDER, matrix_style STYLE,
            matrix_order L_ORDER, matrix_style L_STYLE,
            matrix_order R_ORDER, matrix_style R_STYLE>
  Matrix<T_type, ORDER, STYLE>
  operator* (const Matrix<T_type, L_ORDER, L_STYLE>& lhs,
             const Matrix<T_type, R_ORDER, R_STYLE>& rhs);


  template <matrix_order L_ORDER, matrix_style L_STYLE,
            matrix_order R_ORDER, matrix_style R_STYLE, typename T_type>
  Matrix<T_type, L_ORDER, Concrete>
  operator* (const Matrix<T_type, L_ORDER, L_STYLE>& lhs,
             const Matrix<T_type, R_ORDER, R_STYLE>& rhs);

	/* forward declaration of the matrix class */
	template <typename T_type, matrix_order ORDER, matrix_style STYLE>
	class Matrix;

  /*!  \brief A helper class for list-wise initialization.  
   *
   * This class gets used behind the scenes to provide listwise
   * initialization for Matrix objects.  This documentation is mostly
   * intended for developers.
   *
   * The Matrix class's assignment operator returns a ListInitializer
   * object when passed a scalar.  The assignment operator binds before
   * the comma operator, so this happens first, no matter if there is
   * one scalar, or a list of scalars on the right hand side of the
   * assignment sign.  The ListInitializer constructor keeps an iterator
   * to the Matrix that created it and places the initial item at the
   * head of a list.  Then the ListInitializer comma operator gets
   * called 0 or more times, appending items to the list. At this
   * point the ListInitializer object gets destructed because the
   * expression is done and it is just a temporary.  All the action is
   * in the destructor where the list is copied into the Matrix with
   * R-style recycling.
   *
   * To handle chained assignments, such as A = B = C = 1.2 where A,
   * B, and C are matrices, correctly, we encapsulate the Matrix
   * population sequence that is typically called by the destructor in
   * the private function populate, and make Matrix a friend of this
   * class.  The Matrix class contains an assignment operator for
   * ListInitializer objects that calls this function.  When a call
   * like "A = B = C = 1.2" occurs the compiler first evaluates C =
   * 1.2 and returns a ListInitializer object.  Then, the
   * ListInitializer assignment operator in the Matrix class (being
   * called on B = (C = 1.2)) forces the ListInitializer object to
   * populate C early (it would otherwise not occur until destruction
   * at the end of th entire call) by calling the populate method and
   * then does a simple Matrix assignment of B = C and the populated C
   * and the chain of assignment proceeds from there in the usual
   * fashion.
   *
   * Based on code in Blitz++ (http://www.oonumerics.org/blitz/) by
   * Todd Veldhuizen.  Blitz++ is distributed under the terms of the
   * GNU GPL.
   */

  template<typename T_elem, typename T_iter, 
           matrix_order O, matrix_style S>
  class ListInitializer {
    // An unbound friend
    template <class T, matrix_order OO, matrix_style SS>
    friend class Matrix;
    
    public:
      ListInitializer (T_elem val, T_iter begin, T_iter end, 
                       Matrix<T_elem,O,S>* matrix)
        : vals_ (),
          iter_ (begin),
          end_ (end),
          matrix_ (matrix),
          populated_ (false)
      {
        vals_.push_back(val);
      }

      ~ListInitializer ()
      {
        if (! populated_)
          populate();
      }

      ListInitializer &operator, (T_elem x)
      {
        vals_.push_back(x);
        return *this;
      }

    private:
      void populate ()
      {
        typename std::list<T_elem>::iterator vi = vals_.begin();

        while (iter_ < end_) {
          if (vi == vals_.end())
            vi = vals_.begin();
          *iter_ = *vi;
          ++iter_; ++vi;
        }

        populated_ = true;
      }

      std::list<T_elem> vals_;
      T_iter iter_;
      T_iter end_;
      Matrix<T_elem, O, S>* matrix_;
      bool populated_;
  };
	
  /*! \brief Matrix superclass.
   *
   * The Matrix_base class handles Matrix functionality that doesn't
   * need to be templatized with respect to data type.  This helps
   * reduce code bloat by reducing replication of code for member
   * functions that don't rely on templating.  Furthermore, it
   * hides all of the implementation details of indexing.  The
   * constructors of this class are protected and end-users should
   * always work with full-fledged Matrix objects.
   *
   * The public functions in this class generally provide Matrix
   * metadata like information on dimensionality and size.
	 */

  template <matrix_order ORDER = Col, matrix_style STYLE = Concrete>
  class Matrix_base
  {
    protected:
      /**** CONSTRUCTORS ****/

      /* Default constructor */
      Matrix_base ()
        : rows_ (0),
          cols_ (0),
          rowstride_ (0),
          colstride_ (0),
          storeorder_ (ORDER)
      {}

      /* Standard constructor */
      Matrix_base (uint rows, uint cols)
        : rows_ (rows),
          cols_ (cols),
          storeorder_ (ORDER)
      {
        if (ORDER == Col) {
          rowstride_ = 1;
          colstride_ = rows;
        } else {
          rowstride_ = cols;
          colstride_ = 1;
        }
      }

      /* Copy constructors 
       *
       * The first version handles matrices of the same order and
       * style.  The second handles matrices with different
       * orders/styles.  The same- templates are more specific,
       * so they will always catch same- cases.
       *
       */

      Matrix_base (const Matrix_base &m)
        : rows_ (m.rows()),
          cols_ (m.cols()),
          rowstride_ (m.rowstride()),
          colstride_ (m.colstride())
      {
        if (STYLE == View)
          storeorder_ = m.storeorder();
        else
          storeorder_ = ORDER;
      }

      template <matrix_order O, matrix_style S>
      Matrix_base (const Matrix_base<O, S> &m)
        : rows_ (m.rows()),
          cols_ (m.cols())
      {
        if (STYLE == View) {
          storeorder_ = m.storeorder();
          rowstride_ = m.rowstride();
          colstride_ = m.colstride();
         } else {
          storeorder_ = ORDER;
          if (ORDER == Col) {
            rowstride_ = 1;
            colstride_ = rows_;
          } else {
            rowstride_ = cols_;
            colstride_ = 1;
          }
         }
      }

      /* Submatrix constructor */
      template <matrix_order O, matrix_style S>
      Matrix_base (const Matrix_base<O, S> &m,
          uint x1, uint y1, uint x2, uint y2)
        : rows_ (x2 - x1 + 1),
          cols_ (y2 - y1 + 1),
          rowstride_ (m.rowstride()),
          colstride_ (m.colstride()),
          storeorder_ (m.storeorder())
      {
        /* Submatrices always have to be views, but the whole
         * concrete-view thing is just a policy maintained by the
         * software.  Therefore, we need to ensure this constructor
         * only returns views.  There's no neat way to do it but this
         * is still a compile time check, even if it only reports at
         * run-time.  Of course, we should never get this far.
         */
        if (STYLE == Concrete) {
          SCYTHE_THROW(scythe_style_error, 
              "Tried to construct a concrete submatrix (Matrix_base)!");
        }
      }


      /**** DESTRUCTOR ****/

      ~Matrix_base ()
      {}

      /**** OPERRATORS ****/

      // I'm defining one just to make sure we don't get any subtle
      // bugs from defaults being called.
      Matrix_base& operator=(const Matrix_base& m)
      {
        SCYTHE_THROW(scythe_unexpected_default_error,
          "Unexpected call to Matrix_base default assignment operator");
      }

      /**** MODIFIERS ****/

			/* Make this Matrix_base an exact copy of the matrix passed in.
			 * Just like an assignment operator but can be called from a derived
			 * class w/out making the = optor public and doing something
			 * like
			 * *(dynamic_cast<Matrix_base *> (this)) = M;
			 * in the derived class.
       *
       * Works across styles, but should be used with care
			 */
      template <matrix_order O, matrix_style S>
			inline void mimic (const Matrix_base<O, S> &m)
			{
				rows_ = m.rows();
				cols_ = m.cols();
				rowstride_ = m.rowstride();
				colstride_ = m.colstride();
        storeorder_ = m.storeorder();
			}

      /* Reset the dimensions of this Matrix_base.
       *
       * TODO This function is a bit of an interface weakness.  It
       * assumes a resize always means a fresh matrix (concrete or
       * view) with no slack between dims and strides. This happens to
       * always be the case when it is called, but it tightly couples
       * Matrix_base and extending classes.  Not a big issue (since
       * Matrix is probably the only class that will ever extend this)
       * but something to maybe fix down the road.
       */
			inline void resize (uint rows, uint cols)
			{
				rows_ = rows;
				cols_ = cols;

        if (ORDER == Col) {
          rowstride_ = 1;
          colstride_ = rows;
        } else {
          rowstride_ = cols;
          colstride_ = 1;
        }

        storeorder_ = ORDER;
			}
			
		public:
			/**** ACCESSORS ****/

      /*! \brief Returns the total number of elements in the Matrix.
       *
       * \see rows()
       * \see cols()
       * \see max_size()
       */
			inline uint size () const
			{
				return (rows() * cols());
			}

			/*! \brief Returns the  number of rows in the Matrix.
       *
       * \see size()
       * \see cols()
       */
			inline uint rows() const
			{
				return rows_;
			}

			/*! \brief Returns the number of columns in the Matrix.
       *
       * \see size()
       * \see rows()
       */
			inline uint cols () const
			{
				return cols_;
			}

      /*! \brief Check matrix ordering.
       *
       * This method returns the matrix_order of this Matrix.
       *
       * \see style()
       * \see storeorder()
       */
      inline matrix_order order() const
      {
        return ORDER;
      }

      /*! \brief Check matrix style.
       *
       * This method returns the matrix_style of this Matrix.
       *
       * \see order()
       * \see storeorder()
       */
      inline matrix_style style() const
      {
        return STYLE;
      }

      /*! \brief Returns the storage order of the underlying
       * DataBlock.  
       *
       * In view matrices, the storage order of the data may not be the
       * same as the template ORDER.
       * 
       *
       * \see rowstride()
       * \see colstride()
       * \see order()
       * \see style()
       */
      inline matrix_order storeorder () const
      {
        return storeorder_;
      }

      /*! \brief Returns the in-memory distance between elements in
       * successive rows of the matrix.
       *
       * \see colstride()
       * \see storeorder()
			 */
			inline uint rowstride () const
			{
				return rowstride_;
			}
			
      /*! \brief Returns the in-memory distance between elements in
       * successive columns of the matrix.
       *
       * \see rowstride()
       * \see storeorder()
			 */
      inline uint colstride () const
			{
				return colstride_;
			}

      /*! \brief Returns the maximum possible matrix size.
       *
       * Maximum matrix size is simply the highest available unsigned
       * int on your system.
       *
       * \see size()
			 */
			inline uint max_size () const
			{
				return UINT_MAX;
			}

			/*! \brief Returns true if this Matrix is 1x1.
       *
       * \see isNull()
       */
			inline bool isScalar () const
			{
				return (rows() == 1 && cols() == 1);
			}

			/*! \brief Returns true if this Matrix is 1xm.
       *
       * \see isColVector()
       * \see isVector()
       */
			inline bool isRowVector () const
			{
				return (rows() == 1);
			}
			
			/*! \brief Returns true if this Matrix is nx1.
       *
       * \see isRowVector()
       * \see isVector()
       */
			inline bool isColVector () const
			{
				return (cols() == 1);
			}

			/*! \brief Returns true if this Matrix is nx1 or 1xn.
       *
       * \see isRowVector()
       * \see isColVector()
       */
			inline bool isVector () const
			{
				return (cols() == 1 || rows() == 1);
			}
			
			/*! \brief Returns true if this Matrix is nxn.
       *
       * Null and scalar matrices are both considered square.
       *
       * \see isNull()
       * \see isScalar()
       */
			inline bool isSquare () const
			{
				return (cols() == rows());
			}

      /*! \brief Returns true if this Matrix has 0 elements.
       *  
       *  \see empty()
       *  \see isScalar()
       */
			inline bool isNull () const
			{
				return (rows() == 0);
			}

      /*! \brief Returns true if this Matrix has 0 elements.
       *
       * This function is identical to isNull() but conforms to STL
       * container class conventions.
       *
       * \see isNull()
       */
			inline bool empty () const
			{
				return (rows() == 0);
			}
			

			/**** HELPERS ****/

			/*! \brief Check if an index is in bounds.
       *
       * This function takes a single-argument index into a Matrix and
       * returns true iff that index is within the bounds of the
       * Matrix.  This function is equivalent to the expression:
       * \code
       * i < M.size()
       * \endcode
       * for a given Matrix M.
       *
       * \param i The element index to check.
       *
       * \see inRange(uint, uint)
       */
			inline bool inRange (uint i) const
			{
				return (i < size());
			}

			/*! \brief Check if an index is in bounds.
       *
       * This function takes a two-argument index into a Matrix and
       * returns true iff that index is within the bounds of the
       * Matrix.  This function is equivalent to the expression:
       * \code
       * i < M.rows() && j < M.cols()
       * \endcode
       * for a given Matrix M.
       *
       * \param i The row value of the index to check.
       * \param j The column value of the index to check.
       *
       * \see inRange(uint)
       */
			inline bool inRange (uint i, uint j) const
			{
				return (i < rows() && j < cols());
			}

    protected:
			/* These methods take offsets into a matrix and convert them
			 * into that actual position in the referenced memory block,
			 * taking stride into account.  Protection is debatable.  They
			 * could be useful to outside functions in the library but most
			 * callers should rely on the public () operator in the derived
			 * class or iterators.
       *
       * Note that these are very fast for concrete matrices but not
       * so great for views.  Of course, the type checks are done at
       * compile time.
			 */
			
			/* Turn an index into a true offset into the data. */
			inline uint index (uint i) const
      {
        if (STYLE == View) {
          if (ORDER == Col) {
            uint col = i / rows();
            uint row = i % rows();
            return (index(row, col));
          } else {
            uint row = i / cols();
            uint col = i % cols();
            return (index(row, col));
          }
        } else
          return(i);
      }

			/* Turn an i, j into an index. */
      inline uint index (uint row, uint col) const
      {
        if (STYLE == Concrete) {
          if (ORDER == Col)
				    return (col * rows() + row);
          else
            return (row * cols() + col);
        } else { // view
          if (storeorder_ == Col)
            return (col * colstride() + row);
          else
            return (row * rowstride() + col);
        }
      }

    
    /**** INSTANCE VARIABLES ****/
    protected:
      uint rows_;   // # of rows
      uint cols_;   // # of cols

    private:
      /* The derived class shouldn't have to worry about this
       * implementation detail.
       */
      uint rowstride_;   // the in-memory number of elements from the
      uint colstride_;   // beginning of one column(row) to the next
      matrix_order storeorder_; // The in-memory storage order of this
                                // matrix.  This will always be the
                                // same as ORDER for concrete
                                // matrices but views can look at
                                // matrices with storage orders that
                                // differ from their own.
                                // TODO storeorder is always known at
                                // compile time, so we could probably
                                // add a third template param to deal
                                // with this.  That would speed up
                                // views a touch.  Bit messy maybe.
  };

	/**** MATRIX CLASS ****/

  /*!  \brief An STL-compliant matrix container class.
   *
   * The Matrix class provides a matrix object with an interface similar
   * to standard mathematical notation.  The class provides a number
   * of unary and binary operators for manipulating matrices.
   * Operators provide such functionality as addition, multiplication,
   * and access to specific elements within the Matrix.  One can test
   * two matrices for equality or use provided methods to test the
   * size, shape, or symmetry of a given Matrix.  In addition, we
   * provide a number of facilities for saving, loading, and printing
   * matrices.  Other portions of the library provide functions for
   * manipulating matrices.  Most notably, la.h provides definitions
   * of common linear algebra routines and ide.h defines functions
   * that perform inversion and decomposition.
   * 
   * This Matrix data structure sits at the core of the library.  In
   * addition to standard matrix operations, this class allows
   * multiple matrices to view and modify the same underlying data.
   * This ability provides an elegant way in which to access and
   * modify submatrices such as isolated row vectors and greatly
   * increases the overall flexibility of the class.  In addition, we
   * provide iterators (defined in matrix_random_access_iterator.h,
   * matrix_forward_iterator.h, and matrix_bidirectional_iterator.h)
   * that allow Matrix objects to interact seamlessly with the generic
   * algorithms provided by the Standard Template Library (STL).
   *
   * The Matrix class uses template parameters to define multiple
   * behaviors.  Matrices are templated on data type, matrix_order,
   * and matrix_style.
   *
   * Matrix objects can contain elements of any type.  For the most
   * part, uses will wish to fill their matrices with single or double
   * precision floating point numbers, but matrices of integers,
   * boolean values, complex numbers, and user-defined types are all
   * possible and useful.  Although the basic book-keeping methods in
   * the Matrix class will support virtually any type, certain
   * operators require that one or more mathematical operator be
   * defined for the given type and many of the functions in the wider
   * Scythe library expect, or even demand, matrices containing floating
   * point numbers.
   *
   * There are two possible Matrix element orderings, row- or
   * column-major.  Differences in matrix ordering will be most
   * noticeable at construction time.  Constructors that build matrices
   * from streams or other list-like structures will place elements
   * into the matrix in its given order.  In general, any method that
   * processes a matrix in order will use the given matrix_order.  For
   * the most part, matrices of both orderings should exhibit the same
   * performance, but when a trade-off must be made, we err on the
   * side of column-major ordering.  In one respect, this bias is very
   * strong.  If you enable LAPACK/BLAS support in with the
   * SCYTHE_LAPACK compiler flag, the library will use these optimized
   * fortran routines to perform a number of operations on column
   * major matrices; we provide no LAPACK/BLAS support for row-major
   * matrices.  Operations on matrices with mismatched ordering are
   * legal and supported, but not guaranteed to be as fast as
   * same-order operations, especially when SCYTHE_LAPACK is enabled.
   *
   * There are also two possible styles of Matrix template, concrete
   * and view.  These two types of matrix provide distinct ways in
   * which to interact with an underlying block of data.
   * 
   * Concrete matrices behave like matrices in previous
   * Scythe releases.  They directly encapsulate a block of data and
   * always process it in the same order as it is stored (their
   * matrix_order always matches the underlying storage order).
   * All copy constructions and assignments on
   * concrete matrices make deep copies and it is not possible to use
   * the reference() method to make a concrete Matrix a view of
   * another Matrix.  Furthermore, concrete matrices are guaranteed to
   * have unit stride (That is, adjacent Matrix elements are stored
   * adjacently in memory).  
   *
   * Views, on the other hand, provide references to data blocks.
   * More than one view can look at the same underlying block of data,
   * possibly at different portions of the data at the same time.
   * Furthermore, a view may look at the data block of a concrete
   * matrix, perhaps accessing a single row vector or a small
   * submatrix of a larger matrix.  When you copy construct
   * a view a deep copy is not made, rather the view simply provides
   * access to the extant block of data underlying the copied object.  
   * Furthermore, when
   * you assign to a view, you overwrite the data the view is
   * currently pointing to, rather than generating a new data block.
   * Together, these behaviors allow
   * for matrices that view portions of other matrices
   * (submatrices) and submatrix assignment.  Views do not guarantee
   * unit stride and may even logically access their elements in a
   * different order than they are stored in memory.  Copying between
   * concretes and views is fully supported and often transparent to
   * the user.
   *
   * There is a fundamental trade-off between concrete matrices and
   * views.  Concrete matrices are simpler to work with, but not
   * as flexible as views.  Because they always have unit stride,
   * concrete matrices
   * have fast iterators and access operators but, because they must
   * always be copied deeply, they provide slow copies (for example,
   * copy construction of a Matrix returned from a function wastes
   * cycles).  Views are more flexible but also somewhat more
   * complicated to program with.  Furthermore, because they cannot
   * guarantee unit stride, their iterators and access operations are
   * somewhat slower than those for concrete matrices.  On the other
   * hand, they provide very fast copies.  The average Scythe user may
   * find that she virtually never works with views directly (although
   * they can be quite useful in certain situations) but they provide
   * a variety of functionality underneath the hood of the library and
   * underpin many common operations.
   *
   * \note
   * The Matrix interface is split between two classes: this Matrix
   * class and Matrix_base, which Matrix extends.  Matrix_base
   * includes a range of accessors that provide the programmer with
   * information about such things as the dimensionality of Matrix
   * objects.
   */

	template <typename T_type = double, matrix_order ORDER = Col, 
            matrix_style STYLE = Concrete>
	class Matrix : public Matrix_base<ORDER, STYLE>,
								 public DataBlockReference<T_type>
	{
		public:
			/**** TYPEDEFS ****/

			/* Iterator types */

      /*! \brief Random Access Iterator type.
       *
       * This typedef for matrix_random_access_iterator provides a
       * convenient shorthand for the default, and most general,
       * Matrix iterator type.
       *
       * \see const_iterator
       * \see reverse_iterator
       * \see const_reverse_iterator
       * \see forward_iterator
       * \see const_forward_iterator
       * \see reverse_forward_iterator
       * \see const_reverse_forward_iterator
       * \see bidirectional_iterator
       * \see const_bidirectional_iterator
       * \see reverse_bidirectional_iterator
       * \see const_reverse_bidirectional_iterator
       */
			typedef matrix_random_access_iterator<T_type, ORDER, ORDER, STYLE>
        iterator;

      /*! \brief Const Random Access Iterator type.
       *
       * This typedef for const_matrix_random_access_iterator provides
       * a convenient shorthand for the default, and most general,
       * Matrix const iterator type.
       *
       * \see iterator
       * \see reverse_iterator
       * \see const_reverse_iterator
       * \see forward_iterator
       * \see const_forward_iterator
       * \see reverse_forward_iterator
       * \see const_reverse_forward_iterator
       * \see bidirectional_iterator
       * \see const_bidirectional_iterator
       * \see reverse_bidirectional_iterator
       * \see const_reverse_bidirectional_iterator
       */
			typedef const_matrix_random_access_iterator<T_type, ORDER, ORDER,
                                                  STYLE> const_iterator;

      /*! \brief Reverse Random Access Iterator type.
       *
       * This typedef uses std::reverse_iterator to describe a
       * reversed matrix_random_access_iterator type.  This is the
       * reverse of iterator.
       *
       * \see iterator
       * \see const_iterator
       * \see const_reverse_iterator
       * \see forward_iterator
       * \see const_forward_iterator
       * \see reverse_forward_iterator
       * \see const_reverse_forward_iterator
       * \see bidirectional_iterator
       * \see const_bidirectional_iterator
       * \see reverse_bidirectional_iterator
       * \see const_reverse_bidirectional_iterator
       */
			typedef typename 
        std::reverse_iterator<matrix_random_access_iterator<T_type, 
                              ORDER, ORDER, STYLE> > reverse_iterator;

      /*! \brief Reverse Const Random Access Iterator type.
       *
       * This typedef uses std::reverse_iterator to describe a
       * reversed const_matrix_random_access_iterator type.  This is
       * the reverse of const_iterator.
       *
       * \see iterator
       * \see const_iterator
       * \see reverse_iterator
       * \see forward_iterator
       * \see const_forward_iterator
       * \see reverse_forward_iterator
       * \see const_reverse_forward_iterator
       * \see bidirectional_iterator
       * \see const_bidirectional_iterator
       * \see reverse_bidirectional_iterator
       * \see const_reverse_bidirectional_iterator
       */
			typedef typename
				std::reverse_iterator<const_matrix_random_access_iterator
                              <T_type, ORDER, ORDER, STYLE> >
				const_reverse_iterator;

      /*! \brief Forward Iterator type.
       *
       * This typedef for matrix_forward_iterator provides
       * a convenient shorthand for a fast (when compared to
       * matrix_random_access_iterator) Matrix iterator type.
       *
       * \see iterator
       * \see const_iterator
       * \see reverse_iterator
       * \see const_reverse_iterator
       * \see const_forward_iterator
       * \see reverse_forward_iterator
       * \see const_reverse_forward_iterator
       * \see bidirectional_iterator
       * \see const_bidirectional_iterator
       * \see reverse_bidirectional_iterator
       * \see const_reverse_bidirectional_iterator
       */
      typedef matrix_forward_iterator<T_type, ORDER, ORDER, STYLE>
        forward_iterator;

      /*! \brief Const Forward Iterator type.
       *
       * This typedef for const_matrix_forward_iterator provides a
       * convenient shorthand for a fast (when compared to
       * const_matrix_random_access_iterator) const Matrix iterator
       * type.
       *
       * \see iterator
       * \see const_iterator
       * \see reverse_iterator
       * \see const_reverse_iterator
       * \see forward_iterator
       * \see reverse_forward_iterator
       * \see const_reverse_forward_iterator
       * \see bidirectional_iterator
       * \see const_bidirectional_iterator
       * \see reverse_bidirectional_iterator
       * \see const_reverse_bidirectional_iterator
       */
      typedef const_matrix_forward_iterator<T_type, ORDER, ORDER, STYLE>
        const_forward_iterator;

      /*! \brief Bidirectional Iterator type.
       *
       * This typedef for matrix_bidirectional_iterator provides
       * a convenient shorthand for a compromise (with speed and
       * flexibility between matrix_random_access_iterator and
       * matrix_forward_iterator) Matrix iterator type.
       *
       * \see iterator
       * \see const_iterator
       * \see reverse_iterator
       * \see const_reverse_iterator
       * \see forward_iterator
       * \see const_forward_iterator
       * \see reverse_forward_iterator
       * \see const_reverse_forward_iterator
       * \see const_bidirectional_iterator
       * \see reverse_bidirectional_iterator
       * \see const_reverse_bidirectional_iterator
       */
			typedef matrix_bidirectional_iterator<T_type, ORDER, ORDER, STYLE>
        bidirectional_iterator;

      /*! \brief Const Bidirectional Iterator type.
       *
       * This typedef for const_matrix_bidirectional_iterator provides
       * a convenient shorthand for a compromise (with speed and
       * flexibility between const_matrix_random_access_iterator and
       * const_matrix_forward_iterator) const Matrix iterator type.
       * 
       * \see iterator
       * \see const_iterator
       * \see reverse_iterator
       * \see const_reverse_iterator
       * \see forward_iterator
       * \see const_forward_iterator
       * \see reverse_forward_iterator
       * \see const_reverse_forward_iterator
       * \see bidirectional_iterator
       * \see reverse_bidirectional_iterator
       * \see const_reverse_bidirectional_iterator
       */
			typedef const_matrix_bidirectional_iterator<T_type, ORDER, ORDER,
                                  STYLE> const_bidirectional_iterator;

      /*! \brief Const Bidirectional Iterator type.
       *
       * This typedef uses std::reverse_iterator to describe a
       * reversed matrix_bidirectional_iterator type.  This is
       * the reverse of bidirectional_iterator.
       *
       * \see iterator
       * \see const_iterator
       * \see reverse_iterator
       * \see const_reverse_iterator
       * \see forward_iterator
       * \see const_forward_iterator
       * \see reverse_forward_iterator
       * \see const_reverse_forward_iterator
       * \see bidirectional_iterator
       * \see const_bidirectional_iterator
       * \see const_reverse_bidirectional_iterator
       */
			typedef typename 
        std::reverse_iterator<matrix_bidirectional_iterator<T_type, 
                ORDER, ORDER, STYLE> > reverse_bidirectional_iterator;

      /*! \brief Reverse Const Bidirectional Iterator type.
       *
       * This typedef uses std::reverse_iterator to describe a
       * reversed const_matrix_bidirectional_iterator type.  This is
       * the reverse of const_bidirectional_iterator.
       *
       * \see iterator
       * \see const_iterator
       * \see reverse_iterator
       * \see const_reverse_iterator
       * \see forward_iterator
       * \see const_forward_iterator
       * \see reverse_forward_iterator
       * \see const_reverse_forward_iterator
       * \see bidirectional_iterator
       * \see const_bidirectional_iterator
       * \see reverse_bidirectional_iterator
       */
			typedef typename
				std::reverse_iterator<const_matrix_bidirectional_iterator
                              <T_type, ORDER, ORDER, STYLE> >
				const_reverse_bidirectional_iterator;

      /*!\brief The Matrix' element type.
       *
       * This typedef describes the element type (T_type) of this
       * Matrix.
       */
      typedef T_type ttype;
		
		private:
			/* Some convenience typedefs */
			typedef DataBlockReference<T_type> DBRef;
			typedef Matrix_base<ORDER, STYLE> Base;
			
		public:
			/**** CONSTRUCTORS ****/

			/*! \brief Default constructor.
       *
       * The default constructor creates an empty/null matrix.  Using
       * null matrices in operations will typically cause errors; this
       * constructor exists primarily for initialization within
       * aggregate types.
       *
       * \see Matrix(T_type)
       * \see Matrix(uint, uint, bool, T_type)
       * \see Matrix(uint, uint, T_iterator)
       * \see Matrix(const std::string&)
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see Matrix(const Matrix<T_type, O, S>&, uint, uint, uint, uint)
       *
       * \b Example:
       * \include example.matrix.constructor.default.cc
       */
			Matrix ()
				:	Base (),
					DBRef ()
			{
			}

			/*! \brief Parameterized type constructor.
       *
       * Creates a 1x1 matrix (scalar).
       *
       * \param element The scalar value of the constructed Matrix.
       *
       * \see Matrix()
       * \see Matrix(uint, uint, bool, T_type)
       * \see Matrix(uint, uint, T_iterator)
       * \see Matrix(const std::string&)
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see Matrix(const Matrix<T_type, O, S>&, uint, uint, uint, uint)
       *
       * \throw scythe_alloc_error (Level 1)
       *
       * \b Example:
       * \include example.matrix.constructor.ptype.cc
			 */
			Matrix (T_type element)
				:	Base (1, 1),
					DBRef (1)
			{
				data_[Base::index(0)] = element;  // ALWAYS use index()
			}

      /*! \brief Standard constructor.
       *
       * The standard constructor creates a rowsXcols Matrix, filled
       * with zeros by default.  Optionally, you can leave the Matrix
       * uninitialized, or choose a different fill value.
       * 
       * \param rows The number of rows in the Matrix.
       * \param cols The number of columns in the Matrix.
       * \param fill Indicates whether or not the Matrix should be
       * initialized.
       * \param fill_value The scalar value to fill the Matrix with
       * when fill == true.
       *
       * \see Matrix()
       * \see Matrix(T_type)
       * \see Matrix(uint, uint, T_iterator)
       * \see Matrix(const std::string&)
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see Matrix(const Matrix<T_type, O, S>&, uint, uint, uint, uint)
       *
       * \throw scythe_alloc_error (Level 1)
       *
       * \b Example:
       * \include example.matrix.constructor.standard.cc
       */
			Matrix (uint rows, uint cols, bool fill = true,
					T_type fill_value = 0)
				:	Base (rows, cols),
					DBRef (rows * cols)
			{
        // TODO Might use iterator here for abstraction.
				if (fill)
					for (uint i = 0; i < Base::size(); ++i)
						data_[Base::index(i)] = fill_value; // we know data contig
			}

      /*! \brief Iterator constructor.
       *
			 * Creates a \a rows X \a cols matrix, filling it sequentially
			 * (based on this template's matrix_order) with values
			 * referenced by the input iterator \a it.  Pointers are a form
			 * of input iterator, so one can use this constructor to
			 * initialize a matrix object from a c-style array.  The caller
			 * is responsible for supplying an iterator that won't be
			 * exhausted too soon.
       *
       * \param rows The number of rows in the Matrix.
       * \param cols The number of columns in the Matrix.
       * \param it The input iterator to read from.
       *
       * \see Matrix()
       * \see Matrix(T_type)
       * \see Matrix(uint, uint, bool, T_type)
       * \see Matrix(const std::string&)
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see Matrix(const Matrix<T_type, O, S>&, uint, uint, uint, uint)
       *
       * \throw scythe_alloc_error (Level 1)
       *
       * \b Example:
       * \include example.matrix.constructor.iterator.cc
       */
			template <typename T_iterator>
			Matrix (uint rows, uint cols, T_iterator it)
				:	Base (rows, cols),
					DBRef (rows * cols)
			{
        // TODO again, should probably use iterator
				for (uint i = 0; i < Base::size(); ++i) {
					data_[Base::index(i)] = *it; // we know data_ is contig
					++it;
				}
			}

      /*! \brief File constructor.
       *
       * Constructs a matrix from the contents of a file.  The
       * standard input file format is a simple rectangular text file
       * with one matrix row per line and spaces delimiting values in
       * a row.  Optionally, one can also use Scythe's old file format
       * which is a space-delimited, row-major ordered, list of values
       * with row and column lengths in the first two slots.
       *
       * \param path The path of the input file.
       * \param oldstyle Whether or not to use Scythe's old file format.
       *
       * \see Matrix()
       * \see Matrix(T_type)
       * \see Matrix(uint, uint, bool, T_type)
       * \see Matrix(uint, uint, T_iterator)
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see Matrix(const Matrix<T_type, O, S>&, uint, uint, uint, uint)
       * \see save(const std::string&)
       *
       * \throw scythe_alloc_error (Level 1)
       * \throw scythe_file_error (Level 1)
       * \throw scythe_bounds_error (Level 3)
       *
       * \b Example:
       * \include example.matrix.constructor.file.cc
       */
      Matrix (const std::string& path, bool oldstyle=false)
        : Base (),
          DBRef ()
      {
        std::ifstream file(path.c_str());
        SCYTHE_CHECK_10(! file, scythe_file_error,
            "Could not open file " << path);

        if (oldstyle) {
          uint rows, cols;
          file >> rows >> cols;
          resize(rows, cols);
          std::copy(std::istream_iterator<T_type> (file), 
                    std::istream_iterator<T_type>(), begin_f<Row>());
        } else {
          std::string row;

          unsigned int cols = -1;
          std::vector<std::vector<T_type> > vals;
          unsigned int rows = 0;
          while (std::getline(file, row)) {
            std::vector<T_type> column;
            std::istringstream rowstream(row);
            std::copy(std::istream_iterator<T_type> (rowstream),
                 std::istream_iterator<T_type>(),
                 std::back_inserter(column));

            if (cols == -1)
              cols = (unsigned int) column.size();

            SCYTHE_CHECK_10(cols != column.size(), scythe_file_error,
                "Row " << (rows + 1) << " of input file has "
                << column.size() << " elements, but should have " << cols);

            vals.push_back(column);
            rows++;
          }

          resize(rows, cols);
          for (unsigned int i = 0; i < rows; ++i)
            operator()(i, _) = Matrix<T_type>(1, cols, vals[i].begin());
        }
      }

      /* Copy constructors. Uses template args to set up correct
       * behavior for both concrete and view matrices.  The branches
       * are no-ops and get optimized away at compile time.
       *
       * We have to define this twice because we must explicitly
       * override the default copy constructor; otherwise it is the
       * most specific template in a lot of cases and causes ugliness.
       */

      /*! \brief Default copy constructor.
       *
       * Copy constructing a concrete Matrix makes an exact copy of M
       * in a new data block.  On the other hand, copy constructing a
       * view Matrix generates a new Matrix object that references (or
       * views) M's existing data block.
       *
       * \param M The Matrix to copy or make a view of.
       *
       * \see Matrix()
       * \see Matrix(T_type)
       * \see Matrix(uint, uint, bool, T_type)
       * \see Matrix(uint, uint, T_iterator)
       * \see Matrix(const std::string&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see Matrix(const Matrix<T_type, O, S>&, uint, uint, uint, uint)
       * \see copy()
       * \see copy(const Matrix<T_type, O, S> &)
       * \see reference(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_alloc_error (Level 1)
       *
       * \b Example:
       * \include example.matrix.constructor.copy.cc
       */

      Matrix (const Matrix& M)
				:	Base (M), // this call deals with concrete-view conversions
					DBRef ()
			{
        if (STYLE == Concrete) {
          referenceNew(M.size());
          scythe::copy<ORDER,ORDER>(M, *this);
        } else // STYLE == View
          referenceOther(M);
			}

      /*! \brief Cross order and/or style copy constructor.
       *
       * Copy constructing a concrete Matrix makes an exact copy of M
       * in a new data block.  On the other hand, copy constructing a
       * view Matrix generates a new Matrix object that references (or
       * views) M's existing data block.
       *
       * This version of the copy constructor extends Matrix(const
       * Matrix &) by allowing the user to make concrete copies and
       * views of matrices that have matrix_order or matrix_style that
       * does not match that of the constructed Matrix.  That is, this
       * constructor makes it possible to create views of concrete
       * matrices and concrete copies of views, row-major copies of
       * col-major matrices, and so on.
       *
       * \param M The Matrix to copy or make a view of.
       *
       * \see Matrix()
       * \see Matrix(T_type)
       * \see Matrix(uint, uint, bool, T_type)
       * \see Matrix(uint, uint, T_iterator)
       * \see Matrix(const std::string&)
       * \see Matrix(const Matrix&)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see Matrix(const Matrix<T_type, O, S>&, uint, uint, uint, uint)
       * \see copy()
       * \see copy(const Matrix<T_type, O, S> &)
       * \see reference(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_alloc_error (Level 1)
       *
       * \b Example:
       * \include example.matrix.constructor.crosscopy.cc
       */

      template <matrix_order O, matrix_style S>
			Matrix (const Matrix<T_type, O, S> &M)
				:	Base (M), // this call deals with concrete-view conversions
					DBRef ()
			{
        if (STYLE == Concrete) {
          referenceNew(M.size());
          scythe::copy<ORDER,ORDER> (M, *this);
        } else // STYLE == View
          referenceOther(M);
			}

      /*! \brief Cross type copy constructor
       *
       * The type conversion copy constructor takes a reference to an
       * existing matrix containing elements of a different type than
       * the constructed matrix and creates a copy. This constructor
       * will only work if it is possible to cast elements from the
       * copied matrix to the type of elements in the constructed
       * matrix.
       *
       * This constructor always creates a deep copy of the existing
       * matrix, even if the constructed matrix is a view. It is
       * impossible for a matrix view with one element type to
       * reference the data block of a matrix containing elements of a
       * different type. 
       * 
       * \param M The Matrix to copy.
       *
       * \see Matrix()
       * \see Matrix(T_type)
       * \see Matrix(uint, uint, bool, T_type)
       * \see Matrix(uint, uint, T_iterator)
       * \see Matrix(const std::string&)
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<T_type, O, S>&, uint, uint, uint, uint)
       *
       * \throw scythe_alloc_error (Level 1)
       *
       * \b Example:
       * \include example.matrix.constructor.convcopy.cc
       */
			template <typename S_type, matrix_order O, matrix_style S>
			Matrix (const Matrix<S_type, O, S> &M)
				:	Base(M), // this call deal with concrete-view conversions
					DBRef (M.size())
			{
        scythe::copy<ORDER,ORDER> (M, *this);
			}

      /*! \brief Submatrix constructor
       *
       * The submatrix constructor takes a reference to an existing
       * matrix and a set of indices, and generates a new Matrix
       * object referencing the submatrix described by the indices.
       * One can only construct a submatrix with a view template and
       * this constructor will throw an error if one tries to use it
       * to construct a concrete matrix.
       *
       * \note
       * The submatrix-returning operators provide the same
       * functionality as this constructor with less obtuse syntax.
       * Users should generally employ these methods instead of this
       * constructor.
       *
       * \param M  The Matrix to view.
       * \param x1 The first row coordinate, \a x1 <= \a x2.
       * \param y1 The first column coordinate, \a y1 <= \a y2.
       * \param x2 The second row coordinate, \a x2 > \a x1.
       * \param y2 The second column coordinate, \a y2 > \a y1.
       *
       * \see Matrix()
       * \see Matrix(T_type)
       * \see Matrix(uint, uint, bool, T_type)
       * \see Matrix(uint, uint, T_iterator)
       * \see Matrix(const std::string&)
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see operator()(uint, uint, uint, uint)
       * \see operator()(uint, uint, uint, uint) const
       * \see operator()(all_elements, uint)
       * \see operator()(all_elements, uint) const
       * \see operator()(uint, all_elements)
       * \see operator()(uint, all_elements) const
       * \see reference(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_style_error (Level 0)
       * \throw scythe_alloc_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
      Matrix (const Matrix<T_type, O, S> &M,
          uint x1, uint y1, uint x2, uint y2)
        : Base(M, x1, y1, x2, y2),
          DBRef(M, Base::index(x1, y1))
      {
        /* Submatrices always have to be views, but the whole
         * concrete-view thing is just a policy maintained by the
         * software.  Therefore, we need to ensure this constructor
         * only returns views.  There's no neat way to do it but this
         * is still a compile time check, even if it only reports at
         * run-time.
         */
        if (STYLE == Concrete) {
          SCYTHE_THROW(scythe_style_error, 
              "Tried to construct a concrete submatrix (Matrix)!");
        }
      }

    public:
      /**** DESTRUCTOR ****/

      /*!\brief Destructor. 
       */
      ~Matrix() {}

      /**** COPY/REFERENCE METHODS ****/

			/* Make this matrix a view of another's data. If this matrix's
			 * previous datablock is not viewed by any other object it is
			 * deallocated.  Concrete matrices cannot be turned into views
       * at run-time!  Therefore, we generate an error here if *this
       * is concrete.
			 */

      /*!\brief View another Matrix's data.
       *
       * This modifier makes this matrix a view of another's data.
       * The action detaches the Matrix from its current view; if no
       * other Matrix views the detached DataBlock, it will be
       * deallocated.  
       *
       * Concrete matrices cannot convert into views at
       * run time.  Therefore, it is an error to invoke this method on
       * a concrete Matrix.
       *
       * \param M The Matrix to view.
       *
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see copy()
       * \see copy(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_style_error (Level 0)
       *
       * \b Example:
       * \include example.matrix.reference.cc
       */
      template <matrix_order O, matrix_style S>
			inline void reference (const Matrix<T_type, O, S> &M)
			{
        if (STYLE == Concrete) {
          SCYTHE_THROW(scythe_style_error, 
              "Concrete matrices cannot reference other matrices");
        } else {
          referenceOther(M);
          mimic(M);
        }
			}

      /*!\brief Create a copy of this matrix.
       *
       * Creates a deep copy of this Matrix.  The returned concrete
       * matrix references a newly created DataBlock that contains
       * values that are identical to, but distinct from, the values
       * contained in the original Matrix.
       *
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see copy(const Matrix<T_type, O, S> &)
       * \see reference(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_alloc_error (Level 1)
       *
       * \b Example:
       * \include example.matrix.copy.cc
       */
			inline Matrix<T_type, ORDER, Concrete> copy () const
			{
				Matrix<T_type, ORDER> res (Base::rows(), Base::cols(), false);
				std::copy(begin_f(), end_f(), res.begin_f());

				return res;
			}

			/* Make this matrix a copy of another. The matrix retains its
       * own order and style in this case, because that can't change
       * at run time.
       */
      /*!\brief Make this Matrix a copy of another.
       *
       * Converts this Matrix into a deep copy of another Matrix.
       * This Matrix retains its own matrix_order and matrix_style but
       * contains copies of M's elements and becomes the same size and
       * shape as M.  Calling this method automatically detaches this
       * Matrix from its previous DataBlock before copying.
       *
       * \param M The Matrix to copy.
       *
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see copy()
       * \see reference(const Matrix<T_type, O, S> &)
       * \see detach()
       *
       * \throw scythe_alloc_error (Level 1)
       *
       * \b Example:
       * \include example.matrix.copyother.cc
       */
      template <matrix_order O, matrix_style S>
			inline void copy (const Matrix<T_type, O, S>& M)
			{
				resize2Match(M);
        scythe::copy<ORDER,ORDER> (M, *this);
      }

			/**** INDEXING OPERATORS ****/
			
      /*! \brief Access or modify an element in this Matrix.
       *
       * This indexing operator allows the caller to access or modify
       * the ith (indexed in this Matrix's matrix_order) element of
       * this Matrix, indexed from 0 to n - 1, where n is the number
       * of elements in the Matrix.
       *
       * \param i The index of the element to access/modify.
       *
       * \see operator[](uint) const
       * \see operator()(uint)
       * \see operator()(uint) const
       * \see operator()(uint, uint)
       * \see operator()(uint, uint) const
       * 
       * \throw scythe_bounds_error (Level 3)
       */
			inline T_type& operator[] (uint i)
			{
				SCYTHE_CHECK_30 (! Base::inRange(i), scythe_bounds_error,
						"Index " << i << " out of range");

				return data_[Base::index(i)];
			}
			
      /*! \brief Access an element in this Matrix.
       *
       * This indexing operator allows the caller to access 
       * the ith (indexed in this Matrix's matrix_order) element of
       * this Matrix, indexed from 0 to n - 1, where n is the number
       * of elements in the Matrix.
       *
       * \param i The index of the element to access.
       *
       * \see operator[](uint)
       * \see operator()(uint)
       * \see operator()(uint) const
       * \see operator()(uint, uint)
       * \see operator()(uint, uint) const
       * 
       * \throw scythe_bounds_error (Level 3)
       */
			inline T_type& operator[] (uint i) const
			{
				SCYTHE_CHECK_30 (! Base::inRange(i), scythe_bounds_error,
						"Index " << i << " out of range");

				return data_[Base::index(i)];
			}

      /*! \brief Access or modify an element in this Matrix.
       *
       * This indexing operator allows the caller to access or modify
       * the ith (indexed in this Matrix's matrix_order) element of
       * this Matrix, indexed from 0 to n - 1, where n is the number
       * of elements in the Matrix.
       *
       * \param i The index of the element to access/modify.
       *
       * \see operator[](uint)
       * \see operator[](uint) const
       * \see operator()(uint) const
       * \see operator()(uint, uint)
       * \see operator()(uint, uint) const
       * 
       * \throw scythe_bounds_error (Level 3)
       */
			inline T_type& operator() (uint i)
			{
				SCYTHE_CHECK_30 (! Base::inRange(i), scythe_bounds_error,
						"Index " << i << " out of range");

				return data_[Base::index(i)];
			}
			
      /*! \brief Access an element in this Matrix.
       *
       * This indexing operator allows the caller to access 
       * the ith (indexed in this Matrix's matrix_order) element of
       * this Matrix, indexed from 0 to n - 1, where n is the number
       * of elements in the Matrix.
       *
       * \param i The index of the element to access.
       *
       * \see operator[](uint)
       * \see operator[](uint) const
       * \see operator()(uint)
       * \see operator()(uint, uint)
       * \see operator()(uint, uint) const
       * 
       * \throw scythe_bounds_error (Level 3)
       */
			inline T_type& operator() (uint i) const
			{
				SCYTHE_CHECK_30 (! Base::inRange(i), scythe_bounds_error,
					"Index " << i << " out of range");

				return data_[Base::index(i)];
			}

      /*! \brief Access or modify an element in this Matrix.
       *
       * This indexing operator allows the caller to access or modify
       * the (i, j)th element of
       * this Matrix, where i is an element of 0, 1, ..., rows - 1 and
       * j is an element of 0, 1, ..., columns - 1.
       *
       * \param i The row index of the element to access/modify.
       * \param j The column index of the element to access/modify.
       *
       * \see operator[](uint)
       * \see operator[](uint) const
       * \see operator()(uint)
       * \see operator()(uint) const
       * \see operator()(uint, uint) const
       * 
       * \throw scythe_bounds_error (Level 3)
       */
			inline T_type& operator() (uint i, uint j)
			{
				SCYTHE_CHECK_30 (! Base::inRange(i, j), scythe_bounds_error,
						"Index (" << i << ", " << j << ") out of range");

				return data_[Base::index(i, j)];
			}
				
      /*! \brief Access an element in this Matrix.
       *
       * This indexing operator allows the caller to access 
       * the (i, j)th element of
       * this Matrix, where i is an element of 0, 1, ..., rows - 1 and
       * j is an element of 0, 1, ..., columns - 1.
       *
       * \param i The row index of the element to access.
       * \param j The column index of the element to access.
       *
       * \see operator[](uint)
       * \see operator[](uint) const
       * \see operator()(uint)
       * \see operator()(uint) const
       * \see operator() (uint, uint)
       * 
       * \throw scythe_bounds_error (Level 3)
       */
			inline T_type& operator() (uint i, uint j) const
			{
				SCYTHE_CHECK_30 (! Base::inRange(i, j), scythe_bounds_error,
						"Index (" << i << ", " << j << ") out of range");

				return data_[Base::index(i, j)];
			}

      /**** SUBMATRIX OPERATORS ****/


      /* Submatrices are always views.  An extra (but relatively
       * cheap) copy constructor call is made when mixing and matching
       * orders like
       *
       * Matrix<> A;
       * ...
       * Matrix<double, Row> B = A(2, 2, 4, 4);
       *
       * It is technically possible to get around this, by providing
       * templates of each function of the form
       * template <matrix_order O>
       * Matrix<T_type, O, View> operator() (...)
       *
       * but the syntax to call them (crappy return type inference):
       *
       * Matrix<double, Row> B = A.template operator()<Row>(2, 2, 4, 4)
       *
       * is such complete gibberish that I don't think this is worth
       * the slight optimization.
       */
      
      /*! \brief Returns a view of a submatrix.
       *
       * This operator returns a rectangular submatrix view of this
       * Matrix with its upper left corner at (x1, y1) and its lower
       * right corner at (x2, y2).
       *
       * \param x1 The upper row of the submatrix.
       * \param y1 The leftmost column of the submatrix.
       * \param x2 The lowest row of the submatrix.
       * \param y2 The rightmost column of the submatrix.
       *
       * \see operator()(uint, uint, uint, uint) const
       * \see operator()(all_elements, uint)
       * \see operator()(all_elements, uint) const
       * \see operator()(uint, all_elements)
       * \see operator()(uint, all_elements) const
       *
       * \throw scythe_bounds_error (Level 2)
       *
       * \b Example:
       * \include example.matrix.submatrix.cc
       */
      inline Matrix<T_type, ORDER, View> 
			operator() (uint x1, uint y1, uint x2, uint y2)
			{
				SCYTHE_CHECK_20 (! Base::inRange(x1, y1) 
            || ! Base::inRange(x2, y2)
						|| x1 > x2 || y1 > y2,
						scythe_bounds_error,
						"Submatrix (" << x1 << ", " << y1 << ") ; ("
						<< x2 << ", " << y2 << ") out of range or ill-formed");
				
				return (Matrix<T_type, ORDER, View>(*this, x1, y1, x2, y2));
			}
			
      /*! \brief Returns a view of a submatrix.
       *
       * This operator returns a rectangular submatrix view of this
       * Matrix with its upper left corner at (x1, y1) and its lower
       * right corner at (x2, y2).
       *
       * \param x1 The upper row of the submatrix.
       * \param y1 The leftmost column of the submatrix.
       * \param x2 The lowest row of the submatrix.
       * \param y2 The rightmost column of the submatrix.
       *
       * \see operator()(uint, uint, uint, uint)
       * \see operator()(all_elements, uint)
       * \see operator()(all_elements, uint) const
       * \see operator()(uint, all_elements)
       * \see operator()(uint, all_elements) const
       *
       * \throw scythe_bounds_error (Level 2)
       */
      inline Matrix<T_type, ORDER, View> 
      operator() (uint x1, uint y1, uint x2, uint y2) const
			{
				SCYTHE_CHECK_20 (! Base::inRange(x1, y1) 
            || ! Base::inRange(x2, y2)
						|| x1 > x2 || y1 > y2,
						scythe_bounds_error,
						"Submatrix (" << x1 << ", " << y1 << ") ; ("
						<< x2 << ", " << y2 << ") out of range or ill-formed");

				return (Matrix<T_type, ORDER, View>(*this, x1, y1, x2, y2));
			}

      /*! \brief Returns a view of a column vector.
       *
       * This operator returns a vector view of column j in this Matrix.
       *
       * \param a An all_elements object signifying whole vector access.
       * \param j The column to view.
       *
       * \see operator()(uint, uint, uint, uint)
       * \see operator()(uint, uint, uint, uint) const
       * \see operator()(all_elements, uint) const
       * \see operator()(uint, all_elements)
       * \see operator()(uint, all_elements) const
       *
       * \throw scythe_bounds_error (Level 2)
       *
       * \b Example:
       * \include example.matrix.vector.cc
       */
      inline Matrix<T_type, ORDER, View> 
			operator() (const all_elements a, uint j)
			{
				SCYTHE_CHECK_20 (j >= Base::cols(), scythe_bounds_error,
						"Column vector index " << j << " out of range");

				return (Matrix<T_type, ORDER, View>
           (*this, 0, j, Base::rows() - 1, j));
			}
			
      /*! \brief Returns a view of a column vector.
       *
       * This operator returns a vector view of column j in this Matrix.
       *
       * \param a An all_elements object signifying whole vector access.
       * \param j The column to view.
       *
       * \see operator()(uint, uint, uint, uint)
       * \see operator()(uint, uint, uint, uint) const
       * \see operator()(all_elements, uint)
       * \see operator()(uint, all_elements)
       * \see operator()(uint, all_elements) const
       *
       * \throw scythe_bounds_error (Level 2)
       */
      inline Matrix<T_type, ORDER, View> 
			operator() (const all_elements a, uint j) const
			{
				SCYTHE_CHECK_20 (j >= Base::cols(), scythe_bounds_error,
						"Column vector index " << j << " out of range");

				return (Matrix<T_type, ORDER, View>
           (*this, 0, j, Base::rows() - 1, j));
			}

      /*! \brief Returns a view of a row vector.
       *
       * This operator returns a vector view of row i in this Matrix.
       *
       * \param i The row to view.
       * \param b An all_elements object signifying whole vector access.
       *
       * \see operator()(uint, uint, uint, uint)
       * \see operator()(uint, uint, uint, uint) const
       * \see operator()(all_elements, uint)
       * \see operator()(all_elements, uint) const
       * \see operator()(uint, all_elements) const
       *
       * \throw scythe_bounds_error (Level 2)
       *
       * \b Example:
       * \include example.matrix.vector.cc
       */
      inline Matrix<T_type, ORDER, View> 
			operator() (uint i, const all_elements b)
			{
				SCYTHE_CHECK_20 (i >= Base::rows(), scythe_bounds_error,
						"Row vector index " << i << " out of range");

				return (Matrix<T_type, ORDER, View>
            (*this, i, 0, i, Base::cols() - 1));
			}
			
      /*! \brief Returns a view of a row vector.
       *
       * This operator returns a vector view of row i in this Matrix.
       *
       * \param i The row to view.
       * \param b An all_elements object signifying whole vector access.
       *
       * \see operator()(uint, uint, uint, uint)
       * \see operator()(uint, uint, uint, uint) const
       * \see operator()(all_elements, uint)
       * \see operator()(all_elements, uint) const
       * \see operator()(uint, all_elements)
       *
       * \throw scythe_bounds_error (Level 2)
       */
      inline Matrix<T_type, ORDER, View> 
			operator() (uint i, const all_elements b) const
			{
				SCYTHE_CHECK_20 (i >= Base::rows(), scythe_bounds_error,
						"Row vector index " << i << " out of range");
				return (Matrix<T_type, ORDER, View>
            (*this, i, 0, i, Base::cols() - 1));
			}	

      /**** ASSIGNMENT OPERATORS ****/

       /*
       * As with the copy constructor, we need to
       * explicitly define the same-order-same-style assignment
       * operator or the default operator will take over.
       *
       * TODO With views, it may be desirable to auto-grow (and,
       * technically, detach) views to the null matrix.  This means
       * you can write something like:
       *
       * Matrix<double, Col, View> X;
       * X = ...
       *
       * and not run into trouble because you didn't presize.  Still,
       * not sure this won't encourage silly mistakes...need to think
       * about it.
       */

      /*! \brief Assign the contents of one Matrix to another.
       *
       * Like copy construction, assignment works differently for
       * concrete matrices than it does for views.  When you assign to
       * a concrete Matrix it resizes itself to match the right hand
       * side Matrix and copies over the values.  Like all resizes,
       * this causes this Matrix to detach() from its original
       * DataBlock.  This means that any views attached to this Matrix
       * will no longer view this Matrix's data after the assignment;
       * they will continue to view this Matrix's previous DataBlock.
       * When you assign to a view it first checks that
       * the right hand side conforms to its dimensions (by default,
       * see below), and then copies the right hand side values over
       * into its current DataBlock, overwriting the current contents.
       *
       * Scythe also supports a slightly different model of view
       * assignment.  If the user compiled her program with the
       * SCYTHE_VIEW_ASSIGNMENT_RECYCLE flag set then it is possible
       * to copy into a view that is not of the same size as the
       * Matrix on the right hand side of the equation.  In this case,
       * the operator copies elements from the right hand side
       * object into this matrix until either this matrix runs out of
       * room, or the right hand side one does.  In the latter case,
       * the operator starts over at the beginning of the right hand
       * side object, recycling its values as many times as necessary
       * to fill the left hand side object.  The
       * SCYTHE_VIEW_ASSIGNMENT_RECYCLE flag does not affect the
       * behavior of the concrete matrices in any way.
       *
       * \param M The Matrix to copy.
       *
       * \see operator=(const Matrix<T_type, O, S>&)
       * \see operator=(T_type x)
       * \see operator=(ListInitializer<T_type, ITERATOR, O, S>)
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see copy()
       * \see copy(const Matrix<T_type, O, S> &)
       * \see reference(const Matrix<T_type, O, S> &)
       * \see resize(uint, uint, bool)
       * \see detach()
       *
       * \throw scythe_conformation_error (Level 1)
       * \throw scythe_alloc_error (Level 1)
       *
       * \b Example:
       * \include example.matrix.operator.assignment.cc
       */
      Matrix& operator= (const Matrix& M)
      {
        if (STYLE == Concrete) {
          resize2Match(M);
          scythe::copy<ORDER,ORDER> (M, *this);
        } else {
#ifndef SCYTHE_VIEW_ASSIGNMENT_RECYCLE
          SCYTHE_CHECK_10 (Base::size() != M.size(),
              scythe_conformation_error,
              "LHS has dimensions (" << Base::rows() 
              << ", " << Base::cols()
              << ") while RHS has dimensions (" << M.rows() << ", "
              << M.cols() << ")");

          scythe::copy<ORDER,ORDER> (M, *this);
#else
          copy_recycle<ORDER,ORDER>(M, *this);
#endif
        }

        return *this;
      }
      
      /*! \brief Assign the contents of one Matrix to another.
       *
       * Like copy construction, assignment works differently for
       * concrete matrices than it does for views.  When you assign to
       * a concrete Matrix it resizes itself to match the right hand
       * side Matrix and copies over the values.  Like all resizes,
       * this causes this Matrix to detach() from its original
       * DataBlock.  When you assign to a view it first checks that
       * the right hand side conforms to its dimensions, and then
       * copies the right hand side values over into its current
       * DataBlock, overwriting the current contents.
       *
       * Scythe also supports a slightly different model of view
       * assignment.  If the user compiled her program with the
       * SCYTHE_VIEW_ASSIGNMENT_RECYCLE flag set then it is possible
       * to copy into a view that is not of the same size as the
       * Matrix on the right hand side of the equation.  In this case,
       * the operator copies elements from the right hand side
       * object into this matrix until either this matrix runs out of
       * room, or the right hand side one does.  In the latter case,
       * the operator starts over at the beginning of the right hand
       * side object, recycling its values as many times as necessary
       * to fill the left hand side object.  The
       * SCYTHE_VIEW_ASSIGNMENT_RECYCLE flag does not affect the
       * behavior of the concrete matrices in any way.
       *
       * This version of the assignment operator handles assignments
       * between matrices of different matrix_order and/or
       * matrix_style.
       *
       * \param M The Matrix to copy.
       *
       * \see operator=(const Matrix&)
       * \see operator=(T_type x)
       * \see operator=(ListInitializer<T_type, ITERATOR, O, S>)
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see copy()
       * \see copy(const Matrix<T_type, O, S> &)
       * \see reference(const Matrix<T_type, O, S> &)
       * \see resize(uint, uint, bool)
       * \see detach()
       *
       * \throw scythe_conformation_error (Level 1)
       * \throw scythe_alloc_error (Level 1)
       *
       * \b Example:
       * \include example.matrix.operator.assignment.cc
       */
      template <matrix_order O, matrix_style S>
      Matrix &operator= (const Matrix<T_type, O, S> &M)
      {
        if (STYLE == Concrete) {
          resize2Match(M);
          scythe::copy<ORDER,ORDER> (M, *this);
        } else {
#ifndef SCYTHE_VIEW_ASSIGNMENT_RECYCLE
          SCYTHE_CHECK_10 (Base::size() != M.size(),
              scythe_conformation_error,
              "LHS has dimensions (" << Base::rows() 
              << ", " << Base::cols()
              << ") while RHS has dimensions (" << M.rows() << ", "
              << M.cols() << ")");

          scythe::copy<ORDER,ORDER> (M, *this);
#else
          copy_recycle<ORDER,ORDER>(M, *this);
#endif
        }

        return *this;
      }
      
      /* List-wise initialization behavior is a touch complicated.
       * List needs to be less than or equal to matrix in size and it
       * is copied into the matrix with R-style recycling.
       *
       * The one issue is that, if you want true assignment of a
       * scalar to a concrete matrix (resize the matrix to a scalar)
       * you need to be explicit:
       *
       * Matrix<> A(2, 2);
       * double x = 3;
       * ...
       * A = Matrix<>(x); // -> 3
       *
       * A = x; // -> 3 3
       *        //    3 3
       */

      /*! \brief Copy values in a comma-separated list into this Matrix.
       *
       * This assignment operator allows the user to copy the values in
       * a bare, comma-separated, list into this Matrix.  The list
       * should have no more elements in it than the Matrix has
       * elements.  If the list has fewer elements than the Matrix, it
       * will be recycled until the Matrix is full.
       *
       * If you wish to convert a concrete Matrix to a scalar-valued
       * Matrix object you need to explicitly promote the scalar to a
       * Matrix, using the parameterized type constructor
       * (Matrix(T_type)).
       *
       * \param x The first element in the list.
       *
       * \see operator=(const Matrix&)
       * \see operator=(const Matrix<T_type, O, S>&)
       * \see operator=(ListInitializer<T_type, ITERATOR, O, S>)
       * \see Matrix(const Matrix&)
			 * \see Matrix(const Matrix<T_type, O, S> &)
       * \see Matrix(const Matrix<S_type, O, S> &)
       * \see copy()
       * \see copy(const Matrix<T_type, O, S> &)
       * \see reference(const Matrix<T_type, O, S> &)
       * \see resize(uint, uint, bool)
       * \see detach()
       *
       * \b Example:
       * \include example.matrix.operator.assignment.cc
       */
			ListInitializer<T_type, iterator, ORDER, STYLE> 
      operator=(T_type x)
			{
				return (ListInitializer<T_type, iterator, ORDER, STYLE> 
          (x, begin(),end(), this));
			}

      /*! \brief A special assignment operator.
       *
       * This assignment operator provides the necessary glue to allow
       * chained assignments of matrices where the last assignment is
       * achieved through list initialization.  This allows users to
       * write code like
       * \code
       * Matrix<> A, B, C;
       * Matrix<> D(4, 4, false);
       * A = B = C = (D = 1, 2, 3, 4);
       * \endcode
       * where the assignment in the parentheses technically returns a
       * ListInitializer object, not a Matrix object.  The details of
       * this mechanism are not important for the average user and the
       * distinction can safely be ignored.
       *
       * \note
       * The parentheses in the code above are necessary because of
       * the precedence of the assignment operator.
       *
       * \see operator=(const Matrix&)
       * \see operator=(const Matrix<T_type, O, S>&)
       * \see operator=(T_type x)
       *
       * \b Example:
       * \include example.matrix.operator.assignment.cc
       */
      template <typename ITERATOR, matrix_order O, matrix_style S>
      Matrix &operator=(ListInitializer<T_type, ITERATOR, O, S> li)
      {
        li.populate();
				*this = *(li.matrix_);
        return *this;
      }

      /**** ARITHMETIC OPERATORS ****/

		private:
			/* Reusable chunk of code for element-wise operator
       * assignments.  Updates are done in-place except for the 1x1 by
       * nXm case, which forces a resize.
			 */
			template <typename OP, matrix_order O, matrix_style S>
			inline Matrix&
			elementWiseOperatorAssignment (const Matrix<T_type, O, S>& M, 
                                     OP op)
			{
				SCYTHE_CHECK_10 (Base::size() != 1 && M.size() != 1 && 
						(Base::rows () != M.rows() || Base::cols() != M.cols()),
						scythe_conformation_error,
						"Matrices with dimensions (" << Base::rows() 
            << ", " << Base::cols()
						<< ") and (" << M.rows() << ", " << M.cols()
						<< ") are not conformable");
				
				if (Base::size() == 1) { // 1x1 += nXm
					T_type tmp = (*this)(0);
					resize2Match(M);
          std::transform(M.begin_f<ORDER>(), M.end_f<ORDER>(), 
              begin_f(), std::bind1st(op, tmp));
				} else if (M.size() == 1) { // nXm += 1x1
					std::transform(begin_f(), end_f(), begin_f(),
							std::bind2nd(op, M(0)));
				} else { // nXm += nXm
            std::transform(begin_f(), end_f(), M.begin_f<ORDER>(), 
                begin_f(), op);
        }

				return *this;
			}

    public:
      /*! \brief Add another Matrix to this Matrix.
       *
       * This operator sums this Matrix with another and places the
       * result into this Matrix.  The two matrices must have the same
       * dimensions or one of the matrices must be 1x1.
       *
       * \param M The Matrix to add to this one.
       *
       * \see operator+=(T_type)
       * \see operator-=(const Matrix<T_type, O, S> &)
       * \see operator%=(const Matrix<T_type, O, S> &)
       * \see operator/=(const Matrix<T_type, O, S> &)
       * \see operator^=(const Matrix<T_type, O, S> &)
       * \see operator*=(const Matrix<T_type, O, S> &)
       * \see kronecker(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_conformation_error (Level 1)
       * \throw scythe_alloc_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
			inline Matrix& operator+= (const Matrix<T_type, O, S> &M)
			{
				return elementWiseOperatorAssignment(M, std::plus<T_type> ());
			}

      /*! \brief Add a scalar to this Matrix.
       *
       * This operator sums each element of this Matrix with the
       * scalar \a x and places the result into this Matrix.
       *
       * \param x The scalar to add to each element.
       *
       * \see operator+=(const Matrix<T_type, O, S> &)
       * \see operator-=(T_type)
       * \see operator%=(T_type)
       * \see operator/=(T_type)
       * \see operator^=(T_type)
       * \see kronecker(T_type)
       *
       * \throw scythe_conformation_error (Level 1)
       */
      inline Matrix& operator+= (T_type x)
      {
        return elementWiseOperatorAssignment(Matrix(x), 
            std::plus<T_type> ());
      }
			
      /*! \brief Subtract another Matrix from this Matrix.
       *
       * This operator subtracts another Matrix from this one and
       * places the result into this Matrix.  The two matrices must
       * have the same dimensions or one of the matrices must be 1x1.
       *
       * \param M The Matrix to subtract from this one.
       *
       * \see operator-=(T_type)
       * \see operator+=(const Matrix<T_type, O, S> &)
       * \see operator%=(const Matrix<T_type, O, S> &)
       * \see operator/=(const Matrix<T_type, O, S> &)
       * \see operator^=(const Matrix<T_type, O, S> &)
       * \see operator*=(const Matrix<T_type, O, S> &)
       * \see kronecker(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_conformation_error (Level 1)
       * \throw scythe_alloc_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
			inline Matrix& operator-= (const Matrix<T_type, O, S> &M)
			{
				return elementWiseOperatorAssignment(M, std::minus<T_type> ());
			}

      /*! \brief Subtract a scalar from this Matrix.
       *
       * This operator subtracts \a x from each element of this
       * Matrix and places the result into this Matrix.
       *
       * \param x The scalar to subtract from each element.
       *
       * \see operator-=(const Matrix<T_type, O, S> &)
       * \see operator+=(T_type)
       * \see operator%=(T_type)
       * \see operator/=(T_type)
       * \see operator^=(T_type)
       * \see kronecker(T_type)
       *
       * \throw scythe_conformation_error (Level 1)
       */
      inline Matrix& operator-= (T_type x)
      {
        return elementWiseOperatorAssignment(Matrix(x), 
            std::minus<T_type> ());
      }
			
      /*! \brief Multiply the elements of this Matrix with another's.
       *
       * This operator multiplies the elements of this Matrix with
       * another's and places the result into this Matrix.  The two
       * matrices must have the same dimensions, or one of the
       * matrices must be 1x1.
       *
       * This operator performs element-by-element multiplication
       * (calculates the Hadamard product), not conventional matrix
       * multiplication.
       *
       * \param M The Matrix to multiply with this one.
       *
       * \see operator%=(T_type)
       * \see operator+=(const Matrix<T_type, O, S> &)
       * \see operator-=(const Matrix<T_type, O, S> &)
       * \see operator/=(const Matrix<T_type, O, S> &)
       * \see operator^=(const Matrix<T_type, O, S> &)
       * \see operator*=(const Matrix<T_type, O, S> &)
       * \see kronecker(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_conformation_error (Level 1)
       * \throw scythe_alloc_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
			inline Matrix& operator%= (const Matrix<T_type, O, S> &M)
			{
				return elementWiseOperatorAssignment(M, 
            std::multiplies<T_type> ());
			}

      /*! \brief Multiply this Matrix by a scalar.
       *
       * This operator multiplies each element of this
       * Matrix with \a x and places the result into this Matrix.
       *
       * \param x The scalar to multiply each element by.
       *
       * \see operator%=(const Matrix<T_type, O, S> &)
       * \see operator+=(T_type)
       * \see operator-=(T_type)
       * \see operator/=(T_type)
       * \see operator^=(T_type)
       * \see kronecker(T_type)
       *
       * \throw scythe_conformation_error (Level 1)
       */
      inline Matrix& operator%= (T_type x)
      {
        return elementWiseOperatorAssignment(Matrix(x), 
            std::multiplies<T_type> ());
      }
			
      /*! \brief Divide the elements of this Matrix by another's.
       *
       * This operator divides the elements of this Matrix by
       * another's and places the result into this Matrix.  The two
       * matrices must have the same dimensions, or one of the
       * Matrices must be 1x1.
       *
       * \param M The Matrix to divide this one by.
       *
       * \see operator/=(T_type)
       * \see operator+=(const Matrix<T_type, O, S> &)
       * \see operator-=(const Matrix<T_type, O, S> &)
       * \see operator%=(const Matrix<T_type, O, S> &)
       * \see operator^=(const Matrix<T_type, O, S> &)
       * \see operator*=(const Matrix<T_type, O, S> &)
       * \see kronecker(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_conformation_error (Level 1)
       * \throw scythe_alloc_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
			inline Matrix& operator/= (const Matrix<T_type, O, S> &M)
			{
				return elementWiseOperatorAssignment(M, std::divides<T_type>());
			}

      /*! \brief Divide this Matrix by a scalar.
       *
       * This operator divides each element of this
       * Matrix by \a x and places the result into this Matrix.
       *
       * \param x The scalar to divide each element by.
       *
       * \see operator/=(const Matrix<T_type, O, S> &)
       * \see operator+=(T_type)
       * \see operator-=(T_type)
       * \see operator%=(T_type)
       * \see operator^=(T_type)
       * \see kronecker(T_type)
       *
       * \throw scythe_conformation_error (Level 1)
       */
      inline Matrix& operator/= (T_type x)
      {
        return elementWiseOperatorAssignment(Matrix(x), 
            std::divides<T_type> ());
      }

      /*! \brief Exponentiate the elements of this Matrix by another's.
       *
       * This operator exponentiates the elements of this Matrix by
       * another's and places the result into this Matrix.  The two
       * matrices must have the same dimensions, or one of the
       * Matrices must be 1x1.
       *
       * \param M The Matrix to exponentiate this one by.
       *
       * \see operator^=(T_type)
       * \see operator+=(const Matrix<T_type, O, S> &)
       * \see operator-=(const Matrix<T_type, O, S> &)
       * \see operator%=(const Matrix<T_type, O, S> &)
       * \see operator^=(const Matrix<T_type, O, S> &)
       * \see operator*=(const Matrix<T_type, O, S> &)
       * \see kronecker(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_conformation_error (Level 1)
       * \throw scythe_alloc_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
			inline Matrix& operator^= (const Matrix<T_type, O, S> &M)
			{
				return elementWiseOperatorAssignment(M, 
            exponentiate<T_type>());
			}

      /*! \brief Exponentiate this Matrix by a scalar.
       *
       * This operator exponentiates each element of this
       * Matrix by \a x and places the result into this Matrix.
       *
       * \param x The scalar to exponentiate each element by.
       *
       * \see operator^=(const Matrix<T_type, O, S> &)
       * \see operator+=(T_type)
       * \see operator-=(T_type)
       * \see operator%=(T_type)
       * \see operator/=(T_type)
       * \see kronecker(T_type)
       *
       * \throw scythe_conformation_error (Level 1)
       */
      inline Matrix& operator^= (T_type x)
      {
        return elementWiseOperatorAssignment(Matrix(x), 
            exponentiate<T_type> ());
      }

      /* Matrix mult always disengages views because it generally
       * requires a resize.  We force a disengage in the one place it
       * isn't absolutely necessary(this->size()==1), for consistency.
       */

      /*! \brief Multiply this Matrix by another.
       *
       * This operator multiplies this Matrix by another and places
       * the result into this Matrix.  The two matrices must conform;
       * this Matrix must have as many columns as the right hand side
       * Matrix has rows.
       *
       * Matrix multiplication always causes a Matrix to detach() from
       * its current view, because it generally requires a resize().
       * Even when it is not absolutely necessary to detach() the
       * Matrix, this function will do so to maintain consistency.
       *
       * Scythe will use LAPACK/BLAS routines to multiply concrete
       * column-major matrices of double-precision floating point
       * numbers if LAPACK/BLAS is available and you compile your
       * program with the SCYTHE_LAPACK flag enabled.
       *
       * \param M The Matrix to multiply this one by.
       *
       * \see operator*=(T_type)
       * \see operator+=(const Matrix<T_type, O, S> &)
       * \see operator-=(const Matrix<T_type, O, S> &)
       * \see operator%=(const Matrix<T_type, O, S> &)
       * \see operator/=(const Matrix<T_type, O, S> &)
       * \see operator^=(const Matrix<T_type, O, S> &)
       * \see kronecker(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_conformation_error (Level 1)
       * \throw scythe_alloc_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
			Matrix& operator*= (const Matrix<T_type, O, S>& M)
			{
        /* Farm out the work to the plain old * operator and make this
         * matrix a reference (the only reference) to the result.  We
         * always have to create a new matrix here, so there is no
         * speed-up from using *=.
         */
        
        /* This saves a copy over 
         * *this = (*this) * M;
         * if we're concrete
         */
        Matrix<T_type, ORDER> res = (*this) * M;
        referenceOther(res);
        mimic(res);

				return *this;
			}

      /*! \brief Multiply this Matrix by a scalar.
       *
       * This operator multiplies each element of this
       * Matrix with \a x and places the result into this Matrix.
       *
       * \note This method is identical in behavior to
       * operator%=(T_type).  It also slightly overgeneralizes matrix
       * multiplication but makes life easy on the user by allowing
       * the matrix multiplication operator to work for basic scaler
       * multiplications.
       *
       * \param x The scalar to multiply each element by.
       *
       * \see operator*=(const Matrix<T_type, O, S> &)
       * \see operator+=(T_type)
       * \see operator-=(T_type)
       * \see operator%=(T_type)
       * \see operator/=(T_type)
       * \see operator^=(T_type)
       * \see kronecker(T_type)
       *
       * \throw scythe_conformation_error (Level 1)
       */
      inline Matrix& operator*= (T_type x)
      {
        return elementWiseOperatorAssignment(Matrix(x), 
            std::multiplies<T_type> ());
      }

      /*! \brief Kronecker multiply this Matrix by another.
       *
       * This method computes the Kronecker product of this Matrix and
       * \a M, and sets the value of this Matrix to the result.
       *
       * Kronecker multiplication always causes a Matrix to detach()
       * from its current view, because it generally requires a
       * resize().
       *
       * \note This method would have been implemented as an operator
       * if we had any reasonable operator choices left.
       *
       * \param M The Matrix to Kronecker multiply this one by.
       *
       * \see kronecker(T_type)
       * \see operator+=(const Matrix<T_type, O, S> &)
       * \see operator-=(const Matrix<T_type, O, S> &)
       * \see operator%=(const Matrix<T_type, O, S> &)
       * \see operator/=(const Matrix<T_type, O, S> &)
       * \see operator^=(const Matrix<T_type, O, S> &)
       * \see operator*=(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_alloc_error (Level 1)
       */
      template <matrix_order O, matrix_style S> Matrix& kronecker
        (const Matrix<T_type, O, S>& M) { uint totalrows =
          Base::rows() * M.rows(); uint totalcols = Base::cols() *
            M.cols();
        // Even if we're a view, make this guy concrete.
        Matrix<T_type,ORDER> res(totalrows, totalcols, false);

        /* TODO: This the most natural way to write this in scythe
         * (with a small optimization based on ordering) but probably
         * not the fastest because it uses submatrix assignments.
         * Optimizations should be considered.
         */
        forward_iterator it = begin_f();
        if (ORDER == Row) {
          for (uint row = 0; row < totalrows; row += M.rows()) {
            for (uint col = 0; col < totalcols; col += M.cols()){
              res(row, col, row + M.rows() - 1, col + M.cols() - 1)
                 = (*it) * M;
              it++;
            }
          }
        } else {
          for (uint col = 0; col < totalcols; col += M.cols()) {
            for (uint row = 0; row < totalrows; row += M.rows()){
              res(row, col, row + M.rows() - 1, col + M.cols() - 1)
                = (*it) * M;
              it++;
            }
          }
        }
       
        referenceOther(res);
        mimic(res);

        return *this;
      }
        
      /*! \brief Kronecker multiply this Matrix by a scalar.
       *
       * This method Kronecker multiplies this Matrix with some scalar,
       * \a x.  This is a degenerate case of Kronecker
       * multiplication, simply multiplying every element in the
       * Matrix by \a x.
       *
       * \note This method is identical in behavior to
       * operator%=(T_type) and operator*=(T_type).
       *
       * \param x The scalar to Kronecker multiply this Matrix by.
       *
       * \see kronecker(const Matrix<T_type, O, S> &)
       * \see operator+=(T_type)
       * \see operator-=(T_type)
       * \see operator%=(T_type)
       * \see operator/=(T_type)
       * \see operator^=(T_type)
       * \see operator*=(T_type)
       *
       */
      inline Matrix& kronecker (T_type x)
      {
        return elementWiseOperatorAssignment(Matrix(x), 
            std::multiplies<T_type> ());
      }

      /* Logical assignment operators */

      /*! \brief Logically AND this Matrix with another.
       *
       * This operator computes the element-wise logical AND of this
       * Matrix and another and places the result into this Matrix.
       * That is, after the operation, an element in this Matrix will
       * evaluate to true (or the type-specific analog of true,
       * typically 1) iff the corresponding element previously
       * residing in this Matrix and the corresponding element in \a M
       * both evaluate to true.  The two matrices must have the same
       * dimensions, or one of the Matrices must be 1x1.
       *
       * \param M The Matrix to AND with this one.
       *
       * \see operator&=(T_type)
       * \see operator|=(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_conformation_error (Level 1)
       * \throw scythe_alloc_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
			inline Matrix& operator&= (const Matrix<T_type, O, S> &M)
			{
				return elementWiseOperatorAssignment(M, 
            std::logical_and<T_type>());
			}

      /*! \brief Logically AND this Matrix with a scalar.
       *
       * This operator computes the element-wise logical AND of this
       * Matrix and a scalar.  That is, after the operation, an
       * element in this Matrix will evaluate to true (or the
       * type-specific analog of true, typically 1) iff the
       * corresponding element previously residing in this Matrix and
       * \a x both evaluate to true.
       *
       * \param x The scalar to AND with each element.
       *
       * \see operator&=(const Matrix<T_type, O, S> &)
       * \see operator|=(T_type)
       *
       * \throw scythe_conformation_error (Level 1)
       */
      inline Matrix& operator&= (T_type x)
      {
        return elementWiseOperatorAssignment(Matrix(x), 
            std::logical_and<T_type> ());
      }

      /*! \brief Logically OR this Matrix with another.
       *
       * This operator computes the element-wise logical OR of this
       * Matrix and another and places the result into this Matrix.
       * That is, after the operation, an element in this Matrix will
       * evaluate to true (or the type-specific analog of true,
       * typically 1) if the corresponding element previously
       * residing in this Matrix or the corresponding element in \a M
       * evaluate to true.  The two matrices must have the same
       * dimensions, or one of the Matrices must be 1x1.
       *
       * \param M The Matrix to OR with this one.
       *
       * \see operator|=(T_type)
       * \see operator&=(const Matrix<T_type, O, S> &)
       *
       * \throw scythe_conformation_error (Level 1)
       * \throw scythe_alloc_error (Level 1)
       */
      template <matrix_order O, matrix_style S>
			inline Matrix& operator|= (const Matrix<T_type, O, S> &M)
			{
				return elementWiseOperatorAssignment(M, 
            std::logical_or<T_type>());
			}

      /*! \brief Logically OR this Matrix with a scalar.
       *
       * This operator computes the element-wise logical OR of this
       * Matrix and a scalar.  That is, after the operation, an
       * element in this Matrix will evaluate to true (or the
       * type-specific analog of true, typically 1) if the
       * corresponding element previously residing in this Matrix or
       * \a x evaluate to true.
       *
       * \param x The scalar to OR with each element.
       *
       * \see operator|=(const Matrix<T_type, O, S> &)
       * \see operator&=(T_type)
       *
       * \throw scythe_conformation_error (Level 1)
       */
      inline Matrix& operator|= (T_type x)
      {
        return elementWiseOperatorAssignment(Matrix(x), 
            std::logical_or<T_type> ());
      }

			/**** MODIFIERS ****/

			/* Resize a matrix view.  resize() takes dimensions as
			 * parameters while resize2Match() takes a matrix reference and
			 * uses its dimensions.
			 */

      /*! \brief Resize or reshape a Matrix.
       *
       * This modifier resizes this Matrix to the given dimensions.
       * Matrix contents after a resize is undefined (junk) unless the
       * preserve flag is set to true.  In this case, the old contents
       * of the Matrix remains at the same indices it occupied in the
       * old Matrix.  Any excess capacity is junk.
       *
			 * Resizing a Matrix ALWAYS disengages it from its current view,
			 * even if the dimensions passed to resize are the same as the
			 * current Matrix's dimensions.  Resized matrices point to new,
			 * uninitialized data blocks (technically, the Matrix might
			 * recycle its current block if it is the only Matrix viewing 
			 * the block, but callers cannot rely on this).  It is important
       * to realize that concrete matrices behave just like views in
       * this respect.  Any views to a concrete Matrix will be
       * pointing to a different underlying data block than the
       * concrete Matrix after the concrete Matrix is resized.
       *
       * \param rows The number of rows in the resized Matrix.
       * \param cols The number of columns in the resized Matrix.
       * \param preserve Whether or not to retain the current contents
       * of the Matrix.
       *
       * \see resize2Match(const Matrix<T_type, O, S>&, bool)
       * \see detach()
       *
       * \throw scythe_alloc_error (Level 1)
       */
			void resize (uint rows, uint cols, bool preserve=false)
      {
        if (preserve) {
          /* TODO Optimize this case.  It is rather clunky. */
          Matrix<T_type, ORDER, View> tmp(*this);
          DBRef::referenceNew(rows * cols);
          Base::resize(rows, cols);
          uint min_cols = std::min(Base::cols(), tmp.cols());
          uint min_rows = std::min(Base::rows(), tmp.rows());

          // TODO use iterators here perhaps
          if (ORDER == Col) {
            for (uint j = 0; j < min_cols; ++j)
              for (uint i = 0; i < min_rows; ++i)
                (*this)(i, j) = tmp(i, j);
          } else {
            for (uint i = 0; i < min_rows; ++i)
              for (uint j = 0; j < min_cols; ++j)
                (*this)(i, j) = tmp(i, j);
          }
        } else {
          DBRef::referenceNew(rows * cols);
          Base::resize(rows, cols);
        }
      }

      /*!\brief Resize a Matrix to match another.
       *
       * This modifier resizes this Matrix to match the dimensions of
       * the argument.  In all other respects, it behaves just like
       * resize().
       *
       * \param M The Matrix providing the dimensions to mimic.
       * \param preserve Whether or not to train the current contents
       * of the Matrix.
       *
       * \see resize(uint, uint, bool)
       * \see detach()
       *
       * \throw scythe_alloc_error (Level 1)
       */
      template <typename T, matrix_order O, matrix_style S>
			inline void resize2Match(const Matrix<T, O, S> &M,
                               bool preserve=false)
			{
				resize(M.rows(), M.cols(), preserve);
			}

			/* Copy this matrix to a new datablock in contiguous storage */
      /*! \brief Copy the contents of this Matrix to a new DataBlock.
       *
       * The detach method copies the data viewed by this Matrix to a
       * fresh DataBlock, detaches this Matrix from its old block and
       * attaches it to the new block.  The old DataBlock will be
       * deallocated if no other matrices view the block after this
       * one detaches.
       *
       * This method can be used to ensure that this Matrix is the
       * sole viewer of its DataBlock.  It also ensures that the
       * underlying data is stored contiguously in memory.
       *
       * \see copy()
       * \see resize(uint, uint, bool)
       *
       * \throw scythe_alloc_error (Level 1)
       */
			inline void detach ()
			{
				resize2Match(*this, true);
			}

      /* Swap operator: sort of a dual copy constructor.  Part of the
       * standard STL container interface. We only support swaps
       * between matrices of like order and style because things get
       * hairy otherwise.  The behavior of this for concrete matrices
       * is a little hairy in any case.
       *
       * Matrix<> A, B;
       * ... // fill in A and B
       * Matrix<double, Col, View> v1 = A(_, 1);
       * A.swap(B);
       * Matrix<double, Col, View> v2 = B(_, 1);
       * 
       * v1 == v2; // evaluates to true
       *
       */

      /*! \brief Swap this Matrix with another.
       *
       * This modifier is much like a dual copy constructor and is
       * part of the Standard Template Library (STL) 
       * interface for container objects.  It is only possible to swap
       * two matrices of the same matrix_order and matrix_style.  When
       * two matrices are swapped, they trade their underlying
       * DataBlock and dimensions.  This behavior is perfectly natural
       * for views, but my seem somewhat surprising for concrete
       * matrices.  When two concrete matrices are swapped, any views
       * that referenced either matrices' DataBlock will reference the
       * other matrices' DataBlock after the swap.
       *
       * \param M - The Matrix to swap with.
       */
			inline void swap (Matrix &M)
			{
			  Matrix tmp = *this;

        /* This are just reference() calls, but we do this explicitly
         * here to avoid throwing errors on the concrete case.  While
         * having a concrete matrix reference another matrix is
         * generally a bad idea, it is safe when the referenced matrix
         * is concrete, has the same order, and gets deallocated (or
         * redirected at another block) like here.
         */

        referenceOther(M);
        mimic(M);

        M.referenceOther(tmp);
        M.mimic(tmp);
			}

      /**** ACCESSORS ****/

      /* Accessors that don't access the data itself (that don't rely
       * on T_type) are in Matrix_base
       */

      /* Are all the elements of this Matrix == 0 */

      /*! \brief Returns true if every element in this Matrix equals 0.
       *
       * The return value of this method is undefined for null
       * matrices.
       *
       * \see empty()
       * \see isNull()
       */
      inline bool isZero () const
      {
        const_forward_iterator last = end_f();
        return (last == std::find_if(begin_f(), last, 
          std::bind1st(std::not_equal_to<T_type> (), 0)));
      }

      /* M(i,j) == 0 when i != j */
      /*! \brief Returns true if this Matrix is square and its
       * off-diagonal elements are all 0.
       *
       * The return value of this method is undefined for null
       * matrices.
       *
       * \see isSquare()
       * \see isIdentity()
       * \see isLowerTriangular()
       * \see isUpperTriangular()
       */
      inline bool isDiagonal() const
      {
        if (! Base::isSquare())
          return false;
        /* Always travel in order.  It would be nice to use iterators
         * here, but we'd need to take views and their iterators are
         * too slow at the moment.
         * TODO redo with views and iterators if optimized.
         */
        if (ORDER == Row) {
          for (uint i = 0; i < Base::rows(); ++i) {
            for (uint j = 0; j < Base::cols(); ++j) {
              if (i != j && (*this)(i, j) != 0)
                return false;
            }
          }
        } else { // ORDER == Col
          for (uint j = 0; j < Base::cols(); ++j) {
            for (uint i = 0; i < Base::rows(); ++i) {
              if (i != j && (*this)(i, j) != 0)
                return false;
            }
          }
        }
        return true;
      }

      /* M(I, j) == 0 when i!= j and 1 when i == j */
      /*! \brief Returns true if this Matrix is diagonal and its
       * diagonal elements are all 1s.
       *
       * The return value of this method is undefined for null
       * matrices.
       *
       * \see isSquare()
       * \see isDiagonal()
       * \see isLowerTriangular()
       * \see isUpperTriangular()
       */
      inline bool isIdentity () const
      {
        if (! Base::isSquare())
          return false;
        // TODO redo with views and iterator if optimized
        if (ORDER == Row) {
          for (uint i = 0; i < Base::rows(); ++i) {
            for (uint j = 0; j < Base::cols(); ++j) {
              if (i != j) {
                if ((*this)(i,j) != 0)
                  return false;
              } else if ((*this)(i,j) != 1)
                return false;
            }
          }
        } else { // ORDER == Col
          for (uint j = 0; j < Base::rows(); ++j) {
            for (uint i = 0; i < Base::cols(); ++i) {
              if (i != j) {
                if ((*this)(i,j) != 0)
                  return false;
              } else if ((*this)(i,j) != 1)
                return false;
            }
          }
        }
        return true;
      }

      /* M(i,j) == 0 when i < j */
      /*! \brief Returns true if all of this Matrix's above-diagonal
       * elements equal 0.
       *
       * The return value of this method is undefined for null
       * matrices.
       *
       * \see isDiagonal()
       * \see isUpperTriangular
       */
      inline bool isLowerTriangular () const
      {
        // TODO view+iterator if optimized
        if (ORDER == Row) {
          for (uint i = 0; i < Base::rows(); ++i)
            for (uint j = i + 1; j < Base::cols(); ++j)
              if ((*this)(i,j) != 0)
                return false;
        } else {
          for (uint j = 0; j < Base::cols(); ++j)
            for (uint i = 0; i < j; ++i)
              if ((*this)(i,j) != 0)
                return false;
       }
        return true;
      }

      /* M(i,j) == 0 when i > j */
      /*! \brief Returns true if all of this Matrix's below-diagonal
       * elements equal 0.
       *
       * The return value of this method is undefined for null
       * matrices.
       *
       * \see isDiagonal()
       * \see isLowerTriangular
       */
      inline bool isUpperTriangular () const
      {
        // TODO view+iterator if optimized
        if (ORDER == Row) {
          for (uint i = 0; i < Base::rows(); ++i)
            for (uint j = 0; j < i; ++j)
              if ((*this)(i,j) != 0)
                return false;
        } else {
          for (uint j = 0; j < Base::cols(); ++j)
            for (uint i = j + 1; i < Base::rows(); ++i)
              if ((*this)(i,j) != 0)
                return false;
       }
        return true;
      }

      /*! \brief Returns true if this Matrix is square and has no
       * inverse.
       *
       * \see isSquare()
       * \see operator~()
       */
      inline bool isSingular() const
      {
        if (! Base::isSquare() || Base::isNull())
          return false;
        if ((~(*this)) == (T_type) 0)
          return true;
        return false;
      }

      /* Square and t(M) = M(inv(M) * t(M) == I */
      /*! Returns true if this Matrix is equal to its transpose.
       *
       * A Matrix is symmetric when \f$M^T = M\f$ or, equivalently,
       * \f$M^{-1} M^T = I\f$.  In simple terms, this means that the
       * (i,j)th element of the Matrix is equal to the (j, i)th
       * element for all i, j.
       *
       * \see isSkewSymmetric()
       */
      inline bool isSymmetric () const
      {
        if (! Base::isSquare())
          return false;
        // No point in order optimizing here
        for (uint i = 1; i < Base::rows(); ++i)
          for (uint j = 0; j < i; ++j)
            if ((*this)(i, j) != (*this)(j, i))
              return false;

        return true;
      }

      /* The matrix is square and t(A) = -A */
      /*! Returns true if this Matrix is equal to its negated
       * transpose.
       *
       * A Matrix is skew symmetric when \f$-M^T = M\f$ or,
       * equivalently, \f$-M^{-1} M^T = I\f$.  In simple terms, this
       * means that the (i, j)th element of the Matrix is equal to the
       * negation of the (j, i)th element for all i, j.
       *
       * \see isSymmetric()
       */
      inline bool isSkewSymmetric () const
      {
        if (! Base::isSquare())
          return false;
        // No point in order optimizing here
        for (uint i = 1; i < Base::rows(); ++i)
          for (uint j = 0; j < i; ++j)
            if ((*this)(i, j) != 0 - (*this)(j, i))
              return false;

        return true;
      }

      /*! \brief Test Matrix equality.
       *
       * This method returns true if all of \a M's elements are equal
       * to those in this Matrix.  To be equal, two matrices must
       * be of the same dimension.  Matrices with differing
       * matrix_order or matrix_style may equal one another.
       *
       * \param M The Matrix to test equality with.
       *
       * \see equals(T_type x) const
       * \see operator==(const Matrix<T_type, L_ORDER, L_STYLE>& lhs, const Matrix<T_type, R_ORDER, R_STYLE>& rhs)
       */
      template <matrix_order O, matrix_style S>
      inline bool
      equals(const Matrix<T_type, O, S>& M) const
      {
        if (data_ == M.data_ && STYLE == Concrete && S == Concrete)
          return true;
        else if (data_ == M.data_ && Base::rows() == M.rows() 
                 && Base::cols() == M.cols()) {
          return true;
        } else if (this->isNull() && M.isNull())
          return true;
        else if (Base::rows() != M.rows() || Base::cols() != M.cols())
          return false;

        return std::equal(begin_f(), end_f(),
            M.template begin_f<ORDER>());
      }

      /*! \brief Test Matrix equality.
       *
       * This method returns true if all of the elements in this
       * Matrix are equal to \a x.
       *
       * \param x The scalar value to test equality with.
       *
       * \see equals(const Matrix<T_type, O, S>& M) const
       * \see operator==(const Matrix<T_type, L_ORDER, L_STYLE>& lhs, const Matrix<T_type, R_ORDER, R_STYLE>& rhs)
       */
      inline bool
      equals(T_type x) const
      {
        const_forward_iterator last = end_f();
        return (last == std::find_if(begin_f(), last, 
          std::bind1st(std::not_equal_to<T_type> (), x)));
      }


			/**** OTHER UTILITIES ****/

      /*! \brief Returns a pointer to this Matrix's internal data
       * array.
       *
       * This method returns a pointer to the internal data array
       * contained within the DataBlock that this Matrix references.
       * 
       * \warning It is generally a bad idea to use this method.  We
       * provide it only for convenience.  Please note that, when
       * working with views, the internal data array may not even be
       * stored in this Matrix's matrix_order.  Furthermore, data
       * encapsulated by a view will generally not be contiguous
       * within the data array.  It this is a concrete Matrix,
       * getArray() will always return a pointer to a data array
       * ordered like this Matrix and in contiguous storage.
       */
			inline T_type* getArray () const
			{
				return data_;
			}

      /*! \brief Saves a Matrix to disk.
       *
       * This method writes the contents of this Matrix to the file
       * specified by \a path.  The user can control file overwriting
       * with \a flag.  The parameter \a header controls the output
       * style.  When one sets \a header to true the Matrix is written
       * as a space-separated list of values, with the number of rows
       * and columns placed in the first two positions in the list.
       * If header is set to false, the file is written as a space
       * separated ascii block, with end-of-lines indicating ends of
       * rows.  The Matrix is always written out in row-major order.
       *
       * \param path The name of the file to write.
       * \param flag Overwrite flag taking values 'a': append, 'o':
       * overwrite, or 'n': do not replace.
       * \param header Boolean value indicating whether to write as a
       * flat list with dimension header or as a rectangular block.
       *
       * \see Matrix(const std::string& file)
       * \see operator>>(std::istream& is, Matrix<T,O,S>& M)
       *
       * \throw scythe_invalid_arg (Level 0)
       * \throw scythe_file_error (Level 0)
       */
      inline void
      save (const std::string& path, const char flag = 'n',
            const bool header = false)
      {
        std::ofstream out;
        if (flag == 'n') {
          std::fstream temp(path.c_str(), std::ios::in);
          if (! temp)
            out.open(path.c_str(), std::ios::out);
          else {
            temp.close();
            SCYTHE_THROW(scythe_file_error, "Cannot overwrite file "
                << path << " when flag = n");
          }
        } else if (flag == 'o')
          out.open(path.c_str(), std::ios::out | std::ios::trunc);
        else if (flag == 'a')
          out.open(path.c_str(), std::ios::out | std::ios::app);
        else
          SCYTHE_THROW(scythe_invalid_arg, "Invalid flag: " << flag);

        if (! out)
          SCYTHE_THROW(scythe_file_error, 
              "Could not open file " << path); 

        if (header) {
          out << Base::rows() << " " << Base::cols();
          for (uint i = 0; i < Base::size(); ++i)
            out << " " << (*this)[i];
          out << std::endl;
        } else {
          for (uint i = 0; i < Base::rows(); ++i) {
            for (uint j = 0; j < Base::cols(); ++j)
              out << (*this)(i,j) << " ";
            out << "\n";
          }
        }
        out.close();
      }


			/**** ITERATOR FACTORIES ****/

      /* TODO Write some cpp macro code to reduce this to something
       * manageable.
       */

      /* Random Access Iterator Factories */
      
      /* Generalized versions */

      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a random_access_iterator that
       * points to the first element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      matrix_random_access_iterator<T_type, I_ORDER, ORDER, STYLE>
      begin ()
      {
        return matrix_random_access_iterator<T_type, I_ORDER, ORDER,
                                             STYLE>(*this);
      }
      
      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a
       * const_random_access_iterator that
       * points to the first element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline
      const_matrix_random_access_iterator<T_type, I_ORDER, ORDER, STYLE>
      begin() const
      {
        return const_matrix_random_access_iterator<T_type, I_ORDER,
                                                   ORDER, STYLE>
          (*this);
      }

      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns a
       * matrix_random_access_iterator that
       * points to just after the last element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      matrix_random_access_iterator<T_type, I_ORDER, ORDER, STYLE>
      end ()
      {
        return (begin<I_ORDER>() + Base::size());
      }
      
      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an
       * const_matrix_random_access_iterator that
       * points to just after the last element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      const_matrix_random_access_iterator<T_type, I_ORDER, ORDER, STYLE>
      end () const
      {
        return (begin<I_ORDER>() + Base::size());
      }

      /*! \brief Get a reverse iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns a reverse
       * matrix_random_access_iterator that
       * points to the last element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline std::reverse_iterator<matrix_random_access_iterator<T_type,
                                   I_ORDER, ORDER, STYLE> >
      rbegin()
      {
        return std::reverse_iterator<matrix_random_access_iterator
                                     <T_type, I_ORDER, ORDER, STYLE> > 
               (end<I_ORDER>());
      }

      /*! \brief Get a reverse iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns a reverse
       * const_matrix_random_access_iterator that points to the last
       * element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      std::reverse_iterator<const_matrix_random_access_iterator
                            <T_type, I_ORDER, ORDER, STYLE> > 
      rbegin() const
      {
        return std::reverse_iterator<const_matrix_random_access_iterator
                                     <T_type, I_ORDER, ORDER, STYLE> > 
        (end<I_ORDER>());
      }

      /*! \brief Get a reverse iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a reverse
       * matrix_random_access_iterator
       * that points to the just before the first element in the given
       * Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline std::reverse_iterator<matrix_random_access_iterator
                                   <T_type, I_ORDER, ORDER, STYLE> >
      rend()
      {
        return std::reverse_iterator<matrix_random_access_iterator
                                     <T_type, I_ORDER, ORDER, STYLE> > 
               (begin<I_ORDER>());
      }

      /*! \brief Get a reverse iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a reverse
       * const_matrix_random_access_iterator that points to the just
       * before the first element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      std::reverse_iterator<const_matrix_random_access_iterator
                            <T_type, I_ORDER, ORDER, STYLE> > 
      rend() const
      {
        return std::reverse_iterator<const_matrix_random_access_iterator
                                     <T_type, I_ORDER, ORDER, STYLE> > 
          (begin<I_ORDER>());
      }

      /* Specific versions --- the generalized versions force you
       * choose the ordering explicitly.  These definitions set up
       * in-order iteration as a default */
      
      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a Matrix::iterator that
       * points to the first element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline iterator begin ()
      {
        return iterator(*this);
      }
      
      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a Matrix::const_iterator that
       * points to the first element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline const_iterator begin() const
      {
        return const_iterator (*this);
      }

      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an Matrix::iterator that
       * points to just after the last element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline iterator end ()
      {
        return (begin() + Base::size());
      }
      
      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an Matrix::const_iterator that
       * points to just after the last element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline 
      const_iterator end () const
      {
        return (begin() + Base::size());
      }

      /*! \brief Get a reverse iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns a Matrix::reverse_iterator that
       * points to the last element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline reverse_iterator rbegin()
      {
        return reverse_iterator (end());
      }

      /*! \brief Get a reverse iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns a
       * Matrix::const_reverse_iterator that points to the last
       * element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline const_reverse_iterator rbegin() const
      {
        return const_reverse_iterator (end());
      }

      /*! \brief Get a reverse iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a Matrix::reverse_iterator
       * that points to the just before the first element in the given
       * Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline reverse_iterator rend()
      {
        return reverse_iterator (begin());
      }

      /*! \brief Get a reverse iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a Matrix::const_reverse_iterator
       * that points to the just before the first element in the given
       * Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline const_reverse_iterator rend() const
      {
        return const_reverse_iterator (begin());
      }

      /* Forward Iterator Factories */

      /* Generalized versions */

      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a matrix_forward_iterator that
       * points to the first element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      matrix_forward_iterator<T_type, I_ORDER, ORDER, STYLE>
      begin_f ()
      {
        return matrix_forward_iterator<T_type, I_ORDER, ORDER,
                                             STYLE>(*this);
      }
      
      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a
       * const_matrix_forward_iterator that
       * points to the first element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline
      const_matrix_forward_iterator <T_type, I_ORDER, ORDER, STYLE>
      begin_f () const
      {
        return const_matrix_forward_iterator <T_type, I_ORDER,
                                                   ORDER, STYLE>
          (*this);
      }

      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an matrix_forward_iterator that
       * points to just after the last element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      matrix_forward_iterator<T_type, I_ORDER, ORDER, STYLE>
      end_f ()
      {
        return (begin_f<I_ORDER>().set_end());
      }
      
      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an
       * const_matrix_forward_iterator that points to just after the
       * last element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      const_matrix_forward_iterator<T_type, I_ORDER, ORDER, STYLE>
      end_f () const
      {
        return (begin_f<I_ORDER>().set_end());
      }

      /* Default Versions */
      
      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a Matrix::forward_iterator that
       * points to the first element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline forward_iterator begin_f ()
      {
        return forward_iterator(*this);
      }
      
      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a
       * Matrix::const_forward_iterator that points to the first
       * element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline const_forward_iterator begin_f () const
      {
        return const_forward_iterator (*this);
      }

      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an Matrix::forward_iterator that
       * points to just after the last element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline forward_iterator end_f ()
      {
        return (begin_f().set_end());
      }
      
      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an
       * Matrix::const_forward_iterator that points to just after the
       * last element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline 
      const_forward_iterator end_f () const
      {
        return (begin_f().set_end());
      }

      /* Bidirectional Iterator Factories */

      /* Generalized versions */

      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a
       * matrix_bidirectional_iterator that
       * points to the first element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      matrix_bidirectional_iterator<T_type, I_ORDER, ORDER, STYLE>
      begin_bd ()
      {
        return matrix_bidirectional_iterator<T_type, I_ORDER, ORDER,
                                             STYLE>(*this);
      }
      
      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a
       * const_matrix_bidirectional_iterator that points to the first
       * element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline
      const_matrix_bidirectional_iterator<T_type, I_ORDER, ORDER, STYLE>
      begin_bd () const
      {
        return const_matrix_bidirectional_iterator<T_type, I_ORDER,
                                                   ORDER, STYLE>
          (*this);
      }

      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an
       * matrix_bidirectional_iterator that points to just after the
       * last element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      matrix_bidirectional_iterator<T_type, I_ORDER, ORDER, STYLE>
      end_bd ()
      {
        return (begin_bd<I_ORDER>().set_end());
      }
      
      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an
       * const_matrix_bidirectional_iterator that points to just after
       * the last element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      const_matrix_bidirectional_iterator<T_type, I_ORDER, ORDER, STYLE>
      end_bd () const
      {
        return (begin_bd<I_ORDER>.set_end());
      }

      /*! \brief Get a reverse iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns a reverse
       * matrix_bidirectional_iterator that points to the last element
       * in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline std::reverse_iterator<matrix_bidirectional_iterator<T_type,
                                   I_ORDER, ORDER, STYLE> >
      rbegin_bd ()
      {
        return std::reverse_iterator<matrix_bidirectional_iterator
                                     <T_type, I_ORDER, ORDER, STYLE> > 
               (end_bd<I_ORDER>());
      }

      /*! \brief Get a reverse iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns a reverse
       * const_matrix_bidirectional_iterator that points to the last
       * element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      std::reverse_iterator<const_matrix_bidirectional_iterator
                            <T_type, I_ORDER, ORDER, STYLE> > 
      rbegin_bd () const
      {
        return std::reverse_iterator<const_matrix_bidirectional_iterator
                                     <T_type, I_ORDER, ORDER, STYLE> > 
        (end_bd<I_ORDER>());
      }

      /*! \brief Get a reverse iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a reverse
       * matrix_bidirectional_iterator that points to the just before
       * the first element in the given
       * Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline std::reverse_iterator<matrix_bidirectional_iterator
                                   <T_type, I_ORDER, ORDER, STYLE> >
      rend_bd ()
      {
        return std::reverse_iterator<matrix_bidirectional_iterator
                                     <T_type, I_ORDER, ORDER, STYLE> > 
               (begin_bd<I_ORDER>());
      }

      /*! \brief Get a reverse iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a reverse
       * const_matrix_bidirectional_iterator that points to the just
       * before the first element in the given Matrix.
       *
       * This is a general template of this function.  It allows the
       * user to generate iterators that iterate over the given Matrix
       * in any order through an explicit template instantiation.
       */
      template <matrix_order I_ORDER>
      inline 
      std::reverse_iterator<const_matrix_bidirectional_iterator
                            <T_type, I_ORDER, ORDER, STYLE> > 
      rend_bd () const
      {
        return std::reverse_iterator<const_matrix_bidirectional_iterator
                                     <T_type, I_ORDER, ORDER, STYLE> > 
          (begin_bd<I_ORDER>());
      }

      /* Specific versions --- the generalized versions force you
       * choose the ordering explicitly.  These definitions set up
       * in-order iteration as a default */
      
      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a
       * Matrix::bidirectional_iterator that points to the first
       * element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline bidirectional_iterator begin_bd ()
      {
        return bidirectional_iterator(*this);
      }
      
      /*! \brief Get an iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a
       * Matrix::const_bidirectional_iterator that points to the first
       * element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline const_bidirectional_iterator begin_bd() const
      {
        return const_bidirectional_iterator (*this);
      }

      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an
       * Matrix::bidirectional_iterator that points to just after the
       * last element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline bidirectional_iterator end_bd ()
      {
        return (begin_bd().set_end());
      }
      
      /*! \brief Get an iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns an Matrix::const_bidirectional
       * iterator that points to just after the last element in the
       * given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline 
      const_bidirectional_iterator end_bd () const
      {
        return (begin_bd().set_end());
      }

      /*! \brief Get a reverse iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns a
       * Matrix::reverse_bidirectional_iterator that points to the
       * last element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline reverse_bidirectional_iterator rbegin_bd()
      {
        return reverse_bidirectional_iterator (end_bd());
      }

      /*! \brief Get a reverse iterator pointing to the end of a Matrix.
       *
       * This is a factory that returns a
       * Matrix::const_reverse_bidirectional_iterator that points to
       * the last element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline const_reverse_bidirectional_iterator rbegin_bd () const
      {
        return const_reverse_bidirectional_iterator (end_bd());
      }

      /*! \brief Get a reverse iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a
       * Matrix::reverse_bidirectional_iterator that points to the
       * just before the first element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline reverse_bidirectional_iterator rend_bd ()
      {
        return reverse_bidirectional_iterator (begin_bd());
      }

      /*! \brief Get a reverse iterator pointing to the start of a Matrix.
       *
       * This is a factory that returns a
       * Matrix::const_reverse_bidirectional_iterator that points to
       * the just before the first element in the given Matrix.
       *
       * This is the default template of this function.  It allows the
       * user to generate iterators of a given Matrix without
       * explicitly stating the order of iteration.  The iterator
       * returned by this function always iterates in the same order
       * as the given Matrix' matrix_order.
       */
      inline const_reverse_iterator rend_bd () const
      {
        return const_reverse_bidirectiona_iterator (begin_bd());
      }

		protected:
			/**** INSTANCE VARIABLES ****/

			/* I know the point of C++ is to force you to write 20 times
			 * more code than should be necessary but "using" inherited ivs
       * is just stupid.
			 */
			using DBRef::data_;  // refer to inherited data pointer directly 
			using Base::rows_;   // " # of rows directly
			using Base::cols_;   // " # of cols directly

	}; // end class Matrix

  /**** EXTERNAL OPERATORS ****/

  /* External operators include a range of binary matrix operations
   * such as tests for equality, and arithmetic.  Style
   * (concrete/view) of the returned matrix is that of the left hand
   * side parameter by default
   *
   * There is also a question of the ordering of the returned matrix.
   * We adopt the convention of returning a matrix ordered like that
   * of the left hand side argument, by default.
   *
   * Whenever there is only one matrix argument (lhs is scalar) we use
   * its order and style as the default.
   *
   * A general template version of each operator also exists and users
   * can coerce the return type to whatever they prefer using some
   * ugly syntax; ex:
   *
   * Matrix<> A; ...  Matrix<double, Row> B = operator*<Row,Concrete>
   *                                          (A, A);
   *
   * In general, the matrix class copy constructor will quietly
   * convert whatever matrix template is returned to the type of the
   * matrix it is being copied into on return, but one might want to
   * specify the type for objects that only exist for a second (ex:
   * (operator*<Row,Concrete>(A, A)).begin()).  Also, note that the
   * fact that we return concrete matrices by default does not
   * preclude the user from taking advantage of fast view copies.  It
   * is the template type of the object being copy-constructed that
   * matters---in terms of underlying implementation all matrices are
   * views, concrete matrices just maintain a particular policy.
   *
   * TODO Consider the best type for scalar args to these functions.
   * For the most part, these will be primitives---doubles mostly.
   * Passing these by reference is probably less efficient than
   * passing by value.  But, for user-defined types pass-by-reference
   * might be the way to go and the cost in this case will be much
   * higher than the value-reference trade-off for primitives.  Right
   * now we use pass-by-reference but we might reconsider...
   */

  /**** ARITHMETIC OPERATORS ****/

  /* These macros provide templates for the basic definitions required
   * for all of the binary operators.  Each operator requires 6
   * definitions.  First, a general matrix definition must be
   * provided.  This definition can return a matrix of a different
   * style and order than its arguments but can only be called if its
   * template type is explicitly specified.  The actual logic of the
   * operator should be specified within this function.  The macros
   * provide definitions for the other 5 required templates, one
   * default matrix by matrix, general matrix by scalar, default
   * matrix by scalar, general scalar by matrix, default scalar by
   * matrix.  The default versions call the more general versions with
   * such that they will return concrete matrices with order equal to
   * the left-hand (or only) matrix passed to the default version.
   *
   */

#define SCYTHE_BINARY_OPERATOR_DMM(OP)                                \
  template <matrix_order ORDER, matrix_style L_STYLE,                 \
            matrix_order R_ORDER, matrix_style R_STYLE,               \
            typename T_type>                                          \
  inline Matrix<T_type, ORDER, Concrete>                              \
  OP (const Matrix<T_type, ORDER, L_STYLE>& lhs,                      \
      const Matrix<T_type, R_ORDER, R_STYLE>& rhs)                    \
  {                                                                   \
    return OP <T_type, ORDER, Concrete>(lhs, rhs);                    \
  }

#define SCYTHE_BINARY_OPERATOR_GMS(OP)                                \
  template <typename T_type, matrix_order ORDER, matrix_style STYLE,  \
            matrix_order L_ORDER, matrix_style L_STYLE>               \
  inline Matrix<T_type, ORDER, STYLE>                                 \
  OP (const Matrix<T_type, L_ORDER, L_STYLE>& lhs,                    \
      const typename Matrix<T_type>::ttype &rhs)                       \
  {                                                                   \
    return  (OP <T_type, ORDER, STYLE>                                \
        (lhs, Matrix<T_type, L_ORDER>(rhs)));                         \
  }

#define SCYTHE_BINARY_OPERATOR_DMS(OP)                                \
  template <matrix_order ORDER, matrix_style L_STYLE,                 \
            typename T_type>                                          \
  inline Matrix<T_type, ORDER, Concrete>                              \
  OP (const Matrix<T_type, ORDER, L_STYLE>& lhs,                      \
      const typename Matrix<T_type>::ttype &rhs)                      \
  {                                                                   \
    return (OP <T_type, ORDER, Concrete> (lhs, rhs));                 \
  }
  
#define SCYTHE_BINARY_OPERATOR_GSM(OP)                                \
  template <typename T_type, matrix_order ORDER, matrix_style STYLE,  \
            matrix_order R_ORDER, matrix_style R_STYLE>               \
  inline Matrix<T_type, ORDER, STYLE>                                 \
  OP (const typename Matrix<T_type>::ttype &lhs,                      \
      const Matrix<T_type, R_ORDER, R_STYLE>& rhs) \
  {                                                                   \
    return  (OP <T_type, ORDER, STYLE>                                \
        (Matrix<T_type, R_ORDER>(lhs), rhs));                         \
  }

#define SCYTHE_BINARY_OPERATOR_DSM(OP)                                \
  template <matrix_order ORDER, matrix_style R_STYLE,                 \
            typename T_type>                                          \
  inline Matrix<T_type, ORDER, Concrete>                              \
  OP (const typename Matrix<T_type>::ttype &lhs,                      \
      const Matrix<T_type, ORDER, R_STYLE>& rhs)                      \
  {                                                                   \
    return (OP <T_type, ORDER, Concrete> (lhs, rhs));                 \
  }

#define SCYTHE_BINARY_OPERATOR_DEFS(OP)                               \
  SCYTHE_BINARY_OPERATOR_DMM(OP)                                      \
  SCYTHE_BINARY_OPERATOR_GMS(OP)                                      \
  SCYTHE_BINARY_OPERATOR_DMS(OP)                                      \
  SCYTHE_BINARY_OPERATOR_GSM(OP)                                      \
  SCYTHE_BINARY_OPERATOR_DSM(OP)

  /* Matrix multiplication */
  
  /* General template version. Must be called with operator*<> syntax
   */
 
  /* We provide two symmetric algorithms for matrix multiplication,
   * one for col-major and the other for row-major matrices.  They are
   * designed to minimize cache misses.The decision is based on the
   * return type of the template so, when using matrices of multiple
   * orders, this can get ugly.  These optimizations only really start
   * paying dividends as matrices get big, because cache misses are
   * rare with smaller matrices.
   */

  /*! \brief Multiply two matrices.
   *
   * This operator multiplies the matrices \a lhs and \a rhs together,
   * returning the result in a new Matrix object.  This operator is
   * overloaded to provide both Matrix by Matrix multiplication and
   * Matrix by scalar multiplication.  In the latter case, the scalar
   * on the left- or right-hand side of the operator is promoted to a
   * 1x1 Matrix and then multiplied with the Matrix on the other side
   * of the operator.  In either case, the matrices must conform; that
   * is, the number of columns in the left-hand side argument must
   * equal the number of rows in the right-hand side argument.  The
   * one exception is when one matrix is a scalar.  In this case we
   * allow Matrix by scalar multiplication with the "*" operator that
   * is comparable to element-by-element multiplication of a Matrix by
   * a scalar value, for convenience.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * Scythe will use LAPACK/BLAS routines to multiply concrete
   * column-major matrices of double-precision floating point
   * numbers if LAPACK/BLAS is available and you compile your
   * program with the SCYTHE_LAPACK flag enabled.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \see operator*(const Matrix<T_type, L_ORDER, L_STYLE>& lhs, const Matrix<T_type, R_ORDER, R_STYLE>& rhs)
   * \see operator*(const Matrix<T_type, ORDER, L_STYLE>& lhs, const Matrix<T_type, R_ORDER, R_STYLE>& rhs)
   * \see operator*(const Matrix<T_type, L_ORDER, L_STYLE>& lhs, const T_type& rhs)
   * \see operator*(const Matrix<T_type, ORDER, L_STYLE>& lhs, const T_type& rhs)
   * \see operator*(const T_type& lhs, const Matrix<T_type, R_ORDER, R_STYLE>& rhs)
   * \see operator*(const T_type& lhs, const Matrix<T_type, ORDER, R_STYLE>& rhs)
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */
   
  template <typename T_type, matrix_order ORDER, matrix_style STYLE,
            matrix_order L_ORDER, matrix_style L_STYLE,
            matrix_order R_ORDER, matrix_style R_STYLE>
  inline Matrix<T_type, ORDER, STYLE>
  operator* (const Matrix<T_type, L_ORDER, L_STYLE>& lhs,
             const Matrix<T_type, R_ORDER, R_STYLE>& rhs)
  {
    if (lhs.size() == 1 || rhs.size() == 1)
      return (lhs % rhs);

    SCYTHE_CHECK_10 (lhs.cols() != rhs.rows(),
        scythe_conformation_error,
        "Matrices with dimensions (" << lhs.rows() 
        << ", " << lhs.cols()
        << ") and (" << rhs.rows() << ", " << rhs.cols()
        << ") are not multiplication-conformable");

    Matrix<T_type, ORDER, Concrete> result (lhs.rows(), rhs.cols(), false);

    T_type tmp;
    if (ORDER == Col) { // col-major optimized
     for (uint j = 0; j < rhs.cols(); ++j) {
       for (uint i = 0; i < lhs.rows(); ++i)
        result(i, j) = (T_type) 0;
       for (uint l = 0; l < lhs.cols(); ++l) {
         tmp = rhs(l, j);
         for (uint i = 0; i < lhs.rows(); ++i)
           result(i, j) += tmp * lhs(i, l);
       }
     }
    } else { // row-major optimized
     for (uint i = 0; i < lhs.rows(); ++i) {
       for (uint j = 0; j < rhs.cols(); ++j)
         result(i, j) = (T_type) 0;
       for (uint l = 0; l < rhs.rows(); ++l) {
         tmp = lhs(i, l);
         for (uint j = 0; j < rhs.cols(); ++j)
           result(i, j) += tmp * rhs(l,j);
       }
     }
    }

    SCYTHE_VIEW_RETURN(T_type, ORDER, STYLE, result)
  }

  SCYTHE_BINARY_OPERATOR_DEFS(operator*)

  /*! \brief Kronecker multiply two matrices.
   *
   * This functions computes the Kronecker product of two Matrix
   * objects. This function is overloaded to provide both Matrix by
   * Matrix addition and Matrix by scalar addition.  In the former
   * case, the dimensions of the two matrices must be the same.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */
  template <typename T_type, matrix_order ORDER, matrix_style STYLE,
            matrix_order L_ORDER, matrix_style L_STYLE,
            matrix_order R_ORDER, matrix_style R_STYLE>
  inline Matrix<T_type, ORDER, STYLE>
  kronecker (const Matrix<T_type, L_ORDER, L_STYLE>& lhs,
             const Matrix<T_type, R_ORDER, R_STYLE>& rhs)
  {
    Matrix<T_type,ORDER,Concrete> res = lhs;
    res.kronecker(rhs);
    return (res);
  }

  SCYTHE_BINARY_OPERATOR_DEFS(kronecker)

  /* Macro definition for general return type templates of standard
   * binary operators (this handles, +, -, %, /, but not *)
   */
    
#define SCYTHE_GENERAL_BINARY_OPERATOR(OP,FUNCTOR)                    \
  template <typename T_type, matrix_order ORDER, matrix_style STYLE,  \
            matrix_order L_ORDER, matrix_style L_STYLE,               \
            matrix_order R_ORDER, matrix_style R_STYLE>               \
  inline Matrix<T_type, ORDER, STYLE>                                 \
  OP (const Matrix<T_type, L_ORDER, L_STYLE>& lhs,                    \
      const Matrix<T_type, R_ORDER, R_STYLE>& rhs)                    \
  {                                                                   \
    SCYTHE_CHECK_10(lhs.size() != 1 && rhs.size() != 1 &&             \
        (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()),       \
        scythe_conformation_error,                                    \
        "Matrices with dimensions (" << lhs.rows()                    \
        << ", " << lhs.cols()                                         \
        << ") and (" << rhs.rows() << ", " << rhs.cols()              \
        << ") are not conformable");                                  \
                                                                      \
    if (lhs.size() == 1) {                                            \
      Matrix<T_type,ORDER,Concrete> res(rhs.rows(),rhs.cols(),false); \
      std::transform(rhs.begin_f(), rhs.end_f(),                      \
          res.template begin_f<R_ORDER>(),                            \
          std::bind1st(FUNCTOR <T_type>(), lhs(0)));                  \
      SCYTHE_VIEW_RETURN(T_type, ORDER, STYLE, res)                   \
    }                                                                 \
                                                                      \
    Matrix<T_type,ORDER,Concrete> res(lhs.rows(), lhs.cols(), false); \
                                                                      \
    if (rhs.size() == 1) {                                            \
      std::transform(lhs.begin_f(), lhs.end_f(),                      \
          res.template begin_f<L_ORDER> (),                           \
          std::bind2nd(FUNCTOR <T_type>(), rhs(0)));                  \
    } else {                                                          \
      std::transform(lhs.begin_f(), lhs.end_f(),                      \
          rhs.template begin_f<L_ORDER>(),                            \
          res.template begin_f<L_ORDER>(),                            \
          FUNCTOR <T_type> ());                                       \
    }                                                                 \
                                                                      \
    SCYTHE_VIEW_RETURN(T_type, ORDER, STYLE, res)                     \
  }

  /* Addition operators */

  /*! \fn operator+(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Add two matrices.
   *
   * This operator adds the matrices \a lhs and \a rhs together,
   * returning the result in a new Matrix object.  This operator is
   * overloaded to provide both Matrix by Matrix addition and
   * Matrix by scalar addition.  In the former case, the dimensions of
   * the two matrices must be the same.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_OPERATOR (operator+, std::plus)
  SCYTHE_BINARY_OPERATOR_DEFS (operator+)

  /* Subtraction operators */

  /*! \fn operator-(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Subtract two matrices.
   *
   * This operator subtracts the Matrix \a rhs from \a lhs, returning
   * the result in a new Matrix object.  This operator is overloaded
   * to provide both Matrix by Matrix subtraction and Matrix by scalar
   * subtraction.  In the former case, the dimensions of the two
   * matrices must be the same.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_OPERATOR (operator-, std::minus)
  SCYTHE_BINARY_OPERATOR_DEFS (operator-)

  /* Element-by-element multiplication operators */

  /*! \fn operator%(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Element multiply two matrices.
   *
   * This operator multiplies the elements of the matrices \a lhs and
   * \a rhs together, returning the result in a new Matrix object.
   * This operator is overloaded to provide both Matrix by Matrix
   * element-wise multiplication and Matrix by scalar element-wise
   * multiplication.  In the former case, the dimensions of the two
   * matrices must be the same.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_OPERATOR (operator%, std::multiplies)
  SCYTHE_BINARY_OPERATOR_DEFS(operator%)

  /* Element-by-element division */

  /*! \fn operator/(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Divide two matrices.
   *
   * This operator divides the Matrix \a lhs from \a rhs,
   * returning the result in a new Matrix object.  This operator is
   * overloaded to provide both Matrix by Matrix division and
   * Matrix by scalar division.  In the former case, the dimensions of
   * the two matrices must be the same.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_OPERATOR (operator/, std::divides)
  SCYTHE_BINARY_OPERATOR_DEFS (operator/)

  /* Element-by-element exponentiation */

  /*! \fn operator^(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Exponentiate one Matrix by another.
   *
   * This operator exponentiates the elements of Matrix \a lhs by
   * those in  \a rhs, returning the result in a new Matrix object.
   * This operator is overloaded to provide both Matrix by Matrix
   * exponentiation and Matrix by scalar exponentiation.  In the
   * former case, the dimensions of the two matrices must be the same.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_OPERATOR (operator^, exponentiate)
  SCYTHE_BINARY_OPERATOR_DEFS (operator^)

  /* Negation operators */

  // General return type version
  /*! \brief Negate a Matrix.
   *
   * This unary operator returns the negation of \a M.  This version
   * of the operator is a general template and can provide a Matrix
   * with any matrix_order or matrix_style as its return value.
   *
   * We also provide an overloaded default template that returns a
   * concrete matrix with the same matrix_order as \a M.
   *
   * \param M The Matrix to negate.
   *
   * \throw scythe_alloc_error (Level 1)
   */
  template <typename T_type, matrix_order R_ORDER, matrix_style R_STYLE,
            matrix_order ORDER, matrix_style STYLE>
  inline Matrix<T_type, R_ORDER, R_STYLE>
  operator- (const Matrix<T_type, ORDER, STYLE>& M)
  {
    Matrix<T_type, R_ORDER, Concrete> result(M.rows(), M.cols(), false);
    std::transform(M.template begin_f<ORDER>(), 
                   M.template end_f<ORDER>(), 
                   result.template begin_f<R_ORDER>(),
                   std::negate<T_type> ());
    SCYTHE_VIEW_RETURN(T_type, R_ORDER, R_STYLE, result)
  }
  
  // Default return type version
  template <matrix_order ORDER, matrix_style P_STYLE, typename T_type>
  inline Matrix<T_type, ORDER, Concrete>
  operator- (const Matrix<T_type, ORDER, P_STYLE>& M)
  {
    return operator-<T_type, ORDER, Concrete> (M);
  }

  /* Unary not operators */

  /*! \brief Logically NOT a Matrix.
   *
   * This unary operator returns NOT \a M.  This version of the
   * operator is a general template and can provide a boolean Matrix
   * with any matrix_order or matrix_style as its return value.
   *
   * We also provide a default template for this function that returns
   * a concrete boolean with the same matrix_order as \a M.
   *
   * \param M The Matrix to NOT.
   *
   * \see operator!(const Matrix<T_type, ORDER, P_STYLE>& M)
   *
   * \throw scythe_alloc_error (Level 1)
   */
  template <matrix_order R_ORDER, matrix_style R_STYLE, typename T_type,
            matrix_order ORDER, matrix_style STYLE>
  inline Matrix<bool, R_ORDER, R_STYLE>
  operator! (const Matrix<T_type, ORDER, STYLE>& M)
  {
    Matrix<bool, R_ORDER, Concrete> result(M.rows(), M.cols(), false);
    std::transform(M.template begin_f<ORDER>(), 
                   M.template end_f<ORDER>(), 
                   result.template begin_f<R_ORDER>(),
                   std::logical_not<T_type> ());
    SCYTHE_VIEW_RETURN(T_type, R_ORDER, R_STYLE, result)
  }
  
  // Default return type version
  template <typename T_type, matrix_order ORDER, matrix_style P_STYLE>
  inline Matrix<bool, ORDER, Concrete>
  operator! (const Matrix<T_type, ORDER, P_STYLE>& M)
  {
    return (operator!<ORDER, Concrete> (M));
  }
  /**** COMPARISON OPERATORS ****/

  /* These macros are analogous to those above, except they return
   * only boolean matrices and use slightly different template
   * parameter orderings.  Kind of redundant, but less confusing than
   * making omnibus macros that handle both cases.
   */
#define SCYTHE_GENERAL_BINARY_BOOL_OPERATOR(OP,FUNCTOR)               \
  template <matrix_order ORDER, matrix_style STYLE, typename T_type,  \
            matrix_order L_ORDER, matrix_style L_STYLE,               \
            matrix_order R_ORDER, matrix_style R_STYLE>               \
  inline Matrix<bool, ORDER, STYLE>                                   \
  OP (const Matrix<T_type, L_ORDER, L_STYLE>& lhs,                    \
      const Matrix<T_type, R_ORDER, R_STYLE>& rhs)                    \
  {                                                                   \
    SCYTHE_CHECK_10(lhs.size() != 1 && rhs.size() != 1 &&             \
        (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()),       \
        scythe_conformation_error,                                    \
        "Matrices with dimensions (" << lhs.rows()                    \
        << ", " << lhs.cols()                                         \
        << ") and (" << rhs.rows() << ", " << rhs.cols()              \
        << ") are not conformable");                                  \
                                                                      \
    if (lhs.size() == 1) {                                            \
      Matrix<bool,ORDER,Concrete> res(rhs.rows(),rhs.cols(),false);   \
      std::transform(rhs.begin_f(), rhs.end_f(),                      \
          res.template begin_f<R_ORDER>(),                            \
          std::bind1st(FUNCTOR <T_type>(), lhs(0)));                  \
      SCYTHE_VIEW_RETURN(T_type, ORDER, STYLE, res)                   \
    }                                                                 \
                                                                      \
    Matrix<bool,ORDER,Concrete> res(lhs.rows(), lhs.cols(), false);   \
                                                                      \
    if (rhs.size() == 1) {                                            \
      std::transform(lhs.begin_f(), lhs.end_f(),                      \
          res.template begin_f<L_ORDER> (),                           \
          std::bind2nd(FUNCTOR <T_type>(), rhs(0)));                  \
    } else {                                                          \
      std::transform(lhs.begin_f(), lhs.end_f(),                      \
          rhs.template begin_f<L_ORDER>(),                            \
          res.template begin_f<L_ORDER>(),                            \
          FUNCTOR <T_type> ());                                       \
    }                                                                 \
                                                                      \
    SCYTHE_VIEW_RETURN(T_type, ORDER, STYLE, res)                     \
  }

#define SCYTHE_BINARY_BOOL_OPERATOR_DMM(OP)                           \
  template <typename T_type, matrix_order ORDER, matrix_style L_STYLE,\
            matrix_order R_ORDER, matrix_style R_STYLE>               \
  inline Matrix<bool, ORDER, Concrete>                                \
  OP (const Matrix<T_type, ORDER, L_STYLE>& lhs,                      \
             const Matrix<T_type, R_ORDER, R_STYLE>& rhs)             \
  {                                                                   \
    return OP <ORDER, Concrete>(lhs, rhs);                            \
  }

#define SCYTHE_BINARY_BOOL_OPERATOR_GMS(OP)                           \
  template <matrix_order ORDER, matrix_style STYLE, typename T_type,  \
            matrix_order L_ORDER, matrix_style L_STYLE>               \
  inline Matrix<bool, ORDER, STYLE>                                   \
  OP (const Matrix<T_type, L_ORDER, L_STYLE>& lhs,                    \
      const typename Matrix<T_type>::ttype &rhs)                      \
  {                                                                   \
    return  (OP <ORDER, STYLE>                                        \
        (lhs, Matrix<T_type, L_ORDER>(rhs)));                         \
  }

#define SCYTHE_BINARY_BOOL_OPERATOR_DMS(OP)                           \
  template <typename T_type, matrix_order ORDER, matrix_style L_STYLE>\
  inline Matrix<bool, ORDER, Concrete>                                \
  OP (const Matrix<T_type, ORDER, L_STYLE>& lhs,                      \
      const typename Matrix<T_type>::ttype &rhs)                      \
  {                                                                   \
    return (OP <ORDER, Concrete> (lhs, rhs));                         \
  }
  
#define SCYTHE_BINARY_BOOL_OPERATOR_GSM(OP)                           \
  template <matrix_order ORDER, matrix_style STYLE, typename T_type,  \
            matrix_order R_ORDER, matrix_style R_STYLE>               \
  inline Matrix<bool, ORDER, STYLE>                                   \
  OP (const typename Matrix<T_type>::ttype &lhs,                      \
      const Matrix<T_type, R_ORDER, R_STYLE>& rhs)                    \
  {                                                                   \
    return  (OP <ORDER, STYLE>                                        \
        (Matrix<T_type, R_ORDER>(lhs), rhs));                         \
  }

#define SCYTHE_BINARY_BOOL_OPERATOR_DSM(OP)                           \
  template <typename T_type, matrix_order ORDER, matrix_style R_STYLE>\
  inline Matrix<bool, ORDER, Concrete>                                \
  OP (const typename Matrix<T_type>::ttype &lhs,                      \
      const Matrix<T_type, ORDER, R_STYLE>& rhs)                      \
  {                                                                   \
    return (OP <ORDER, Concrete> (lhs, rhs));                         \
  }

#define SCYTHE_BINARY_BOOL_OPERATOR_DEFS(OP)                          \
  SCYTHE_BINARY_BOOL_OPERATOR_DMM(OP)                                 \
  SCYTHE_BINARY_BOOL_OPERATOR_GMS(OP)                                 \
  SCYTHE_BINARY_BOOL_OPERATOR_DMS(OP)                                 \
  SCYTHE_BINARY_BOOL_OPERATOR_GSM(OP)                                 \
  SCYTHE_BINARY_BOOL_OPERATOR_DSM(OP)

  /* Element-wise Equality operator
   * See equals() method for straight equality checks
   */

  /*! \fn operator==(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Test Matrix equality.
   *
   * This operator compares the elements of \a lhs and \a rhs and
   * returns a boolean Matrix of true and false values, indicating
   * whether each pair of compared elements is equal.  This operator
   * is overloaded to provide both Matrix by Matrix equality testing
   * and Matrix by scalar equality testing.  In the former case, the
   * dimensions of the two matrices must be the same.  The boolean
   * Matrix returned has the same dimensions as \a lhs and \a rhs, or
   * matches the dimensionality of the larger Matrix object when one
   * of the two parameters is a scalar or a 1x1 Matrix.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_BOOL_OPERATOR (operator==, std::equal_to)
  SCYTHE_BINARY_BOOL_OPERATOR_DEFS (operator==)

  /*! \fn operator!=(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Test Matrix equality.
   *
   * This operator compares the elements of \a lhs and \a rhs and
   * returns a boolean Matrix of true and false values, indicating
   * whether each pair of compared elements is not equal.  This operator
   * is overloaded to provide both Matrix by Matrix inequality testing
   * and Matrix by scalar inequality testing.  In the former case, the
   * dimensions of the two matrices must be the same.  The boolean
   * Matrix returned has the same dimensions as \a lhs and \a rhs, or
   * matches the dimensionality of the larger Matrix object when one
   * of the two parameters is a scalar or a 1x1 Matrix.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_BOOL_OPERATOR (operator!=, std::not_equal_to)
  SCYTHE_BINARY_BOOL_OPERATOR_DEFS (operator!=)

  /*! \fn operator<(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Test Matrix inequality.
   *
   * This operator compares the elements of \a lhs and \a rhs and
   * returns a boolean Matrix of true and false values, indicating
   * whether each of the left-hand side elements is less than its
   * corresponding right hand side element.  This operator is
   * overloaded to provide both Matrix by Matrix inequality testing
   * and Matrix by scalar inequality testing.  In the former case, the
   * dimensions of the two matrices must be the same.  The boolean
   * Matrix returned has the same dimensions as \a lhs and \a rhs, or
   * matches the dimensionality of the larger Matrix object when one
   * of the two parameters is a scalar or a 1x1 Matrix.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_BOOL_OPERATOR (operator<, std::less)
  SCYTHE_BINARY_BOOL_OPERATOR_DEFS (operator<)

  /*! \fn operator<=(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Test Matrix inequality.
   *
   * This operator compares the elements of \a lhs and \a rhs and
   * returns a boolean Matrix of true and false values, indicating
   * whether each of the left-hand side elements is less than 
   * or equal to its
   * corresponding right hand side element.  This operator is
   * overloaded to provide both Matrix by Matrix inequality testing
   * and Matrix by scalar inequality testing.  In the former case, the
   * dimensions of the two matrices must be the same.  The boolean
   * Matrix returned has the same dimensions as \a lhs and \a rhs, or
   * matches the dimensionality of the larger Matrix object when one
   * of the two parameters is a scalar or a 1x1 Matrix.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_BOOL_OPERATOR (operator<=, std::less_equal)
  SCYTHE_BINARY_BOOL_OPERATOR_DEFS (operator<=)

  /*! \fn operator>(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Test Matrix inequality.
   *
   * This operator compares the elements of \a lhs and \a rhs and
   * returns a boolean Matrix of true and false values, indicating
   * whether each of the left-hand side elements is greater than its
   * corresponding right hand side element.  This operator is
   * overloaded to provide both Matrix by Matrix inequality testing
   * and Matrix by scalar inequality testing.  In the former case, the
   * dimensions of the two matrices must be the same.  The boolean
   * Matrix returned has the same dimensions as \a lhs and \a rhs, or
   * matches the dimensionality of the larger Matrix object when one
   * of the two parameters is a scalar or a 1x1 Matrix.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_BOOL_OPERATOR (operator>, std::greater)
  SCYTHE_BINARY_BOOL_OPERATOR_DEFS (operator>)

  /*! \fn operator>=(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Test Matrix inequality.
   *
   * This operator compares the elements of \a lhs and \a rhs and
   * returns a boolean Matrix of true and false values, indicating
   * whether each of the left-hand side elements is greater than 
   * or equal to its
   * corresponding right hand side element.  This operator is
   * overloaded to provide both Matrix by Matrix inequality testing
   * and Matrix by scalar inequality testing.  In the former case, the
   * dimensions of the two matrices must be the same.  The boolean
   * Matrix returned has the same dimensions as \a lhs and \a rhs, or
   * matches the dimensionality of the larger Matrix object when one
   * of the two parameters is a scalar or a 1x1 Matrix.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_BOOL_OPERATOR (operator>=, std::greater_equal)
  SCYTHE_BINARY_BOOL_OPERATOR_DEFS (operator>=)

  /*! \fn operator&(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Logically AND two matrices.
   *
   * This operator logically ANDs the elements of \a lhs and \a rhs
   * and returns a boolean Matrix of true and false values, with true
   * values in each position where both matrices' elements evaluate to
   * true (or the type specific analog to true, typically any non-zero
   * value).  This operator is overloaded to provide both Matrix by
   * Matrix AND and Matrix by scalar AND.  In the former case, the
   * dimensions of the two matrices must be the same.  The boolean
   * Matrix returned has the same dimensions as \a lhs and \a rhs, or
   * matches the dimensionality of the larger Matrix object when one
   * of the two parameters is a scalar or a 1x1 Matrix.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_BOOL_OPERATOR (operator&, std::logical_and)
  SCYTHE_BINARY_BOOL_OPERATOR_DEFS (operator&)


  /*! \fn operator|(const Matrix<T_type,L_ORDER,L_STYLE>&lhs,
   *                const Matrix<T_type,R_ORDER,R_STYLE>&rhs)
   *
   * \brief Logically OR two matrices.
   *
   * This operator logically ORs the elements of \a lhs and \a rhs
   * and returns a boolean Matrix of true and false values, with true
   * values in each position where either Matrix's elements evaluate to
   * true (or the type specific analog to true, typically any non-zero
   * value).  This operator is overloaded to provide both Matrix by
   * Matrix OR and Matrix by scalar OR.  In the former case, the
   * dimensions of the two matrices must be the same.  The boolean
   * Matrix returned has the same dimensions as \a lhs and \a rhs, or
   * matches the dimensionality of the larger Matrix object when one
   * of the two parameters is a scalar or a 1x1 Matrix.
   *
   * In addition, we define multiple templates of the overloaded
   * operator to provide maximal flexibility when working with
   * matrices with differing matrix_order and/or matrix_style.  Each
   * version of the overloaded operator (Matrix by Matrix, scalar by
   * Matrix, and Matrix by scalar) provides both a default and
   * general behavior, using templates.  By default, the returned
   * Matrix is concrete and has the same matrix_order as the
   * left-hand (or only) Matrix argument.  Alternatively, one may
   * coerce the matrix_order and matrix_style of the returned Matrix
   * to preferred values by using the full template declaration of
   * the operator.
   *
   * \param lhs The left-hand-side Matrix or scalar.
   * \param rhs The right-hand-side Matrix or scalar.
   *
   * \throw scythe_conformation_error (Level 1)
   * \throw scythe_alloc_error (Level 1)
   */

  SCYTHE_GENERAL_BINARY_BOOL_OPERATOR (operator|, std::logical_or)
  SCYTHE_BINARY_BOOL_OPERATOR_DEFS (operator|)

  /**** INPUT-OUTPUT ****/


  /* This function simply copies values from an input stream into a
   * matrix.  It relies on the iterators for bounds checking.
   */

  /*! \brief Populate a Matrix from a stream.
   *
   * This operator reads values from a stream and enters them into an
   * existing Matrix in order.
   *
   * \param is The istream to read from.
   * \param M  The Matrix to populate.
   *
   * \see operator<<(std::ostream& os, const Matrix<T,O,S>& M)
   * \see Matrix::Matrix(const std::string& file)
   *
   * \throw scythe_bounds_error (Level 3)
   */
  template <typename T, matrix_order O, matrix_style S>
  std::istream& operator>> (std::istream& is, Matrix<T,O,S>& M)
  {
    std::copy(std::istream_iterator<T> (is), std::istream_iterator<T>(),
         M.begin_f());

    return is;
  }

  /* Writes a matrix to an ostream in readable format.  This is
   * intended to be used to pretty-print to the terminal.
   */

  /*!\brief Write a Matrix to a stream.
   *
   * Writes a matrix to an ostream in a column-aligned format.  This
   * operator is primarily intended for pretty-printing to the
   * terminal and uses two passes in order to correctly align the
   * output.  If you wish to write a Matrix to disk, Matrix::save() is
   * probably a better option.
   *
   * \param os The ostream to write to.
   * \param M  The Matrix to write out.
   *
   * \see operator>>(std::istream& is, Matrix<T,O,S>& M)
   * \see Matrix::save()
   */
  template <typename T, matrix_order O, matrix_style S>
  std::ostream& operator<< (std::ostream& os, const Matrix<T,O,S>& M)
  {
    /* This function take two passes to figure out appropriate field
     * widths.  Speed isn't really the point here.
     */

    // Store previous io settings
    std::ios_base::fmtflags preop = os.flags();

    uint mlen = os.width();
    std::ostringstream oss;
    oss.precision(os.precision());
    oss << std::setiosflags(std::ios::fixed);
    
    typename Matrix<T,O,S>::const_forward_iterator last = M.end_f();
    for (typename Matrix<T,O,S>::const_forward_iterator i = M.begin_f();
        i != last; ++i) {
      oss.str("");
      oss << (*i);
      if (oss.str().length() > mlen)
        mlen = oss.str().length();
    }


    /* Write the stream */
    // Change to a fixed with format.  Users should control precision
    os << std::setiosflags(std::ios::fixed);

    
    for (uint i = 0; i < M.rows(); ++i) {
      Matrix<T, O, View> row = M(i, _);
      //for (uint i = 0; i < row.size(); ++i)
      //  os << std::setw(mlen) << row[i] << " ";
      typename Matrix<T,O,View>::const_forward_iterator row_last 
        = row.end_f();
      for (typename 
          Matrix<T,O,View>::forward_iterator el = row.begin_f();
          el != row_last; ++el) {
        os << std::setw(mlen) << *el << " ";
      }
      os << std::endl;
    }
    
    
    // Restore pre-op flags
    os.flags(preop);

    return os;
  }

#ifdef SCYTHE_LAPACK
  /* A template specialization of operator* for col-major, concrete
   * matrices of doubles that is only visible when a LAPACK library is
   * available.  This function is an analog of the above function and
   * the above doxygen documentation serves for both.
   *
   * This needs to go below % so it can see the template definition
   * (since it isn't actually in the template itself.
   */

  template<>
  Matrix<>
  operator*<double,Col,Concrete,Col,Concrete>
  (const Matrix<>& lhs, const Matrix<>& rhs)
  {
    if (lhs.size() == 1 || rhs.size() == 1)
      return (lhs % rhs);

    SCYTHE_DEBUG_MSG("Using lapack/blas for matrix multiplication");
    SCYTHE_CHECK_10 (lhs.cols() != rhs.rows(),
        scythe_conformation_error,
        "Matrices with dimensions (" << lhs.rows() 
        << ", " << lhs.cols()
        << ") and (" << rhs.rows() << ", " << rhs.cols()
        << ") are not multiplication-conformable");

    Matrix<> result (lhs.rows(), rhs.cols(), false);

    // Get pointers to the internal arrays and set up some vars
    double* lhspnt = lhs.getArray();
    double* rhspnt = rhs.getArray();
    double* resultpnt = result.getArray();
    const double one(1.0);
    const double zero(0.0);
    int rows = (int) lhs.rows();
    int cols = (int) rhs.cols();
    int innerDim = (int) rhs.rows();

    // Call the lapack routine.
    lapack::dgemm_("N", "N", &rows, &cols, &innerDim, &one, lhspnt,
                   &rows, rhspnt, &innerDim, &zero, resultpnt, &rows);

    return result;
  }
#endif

} // end namespace scythe

#endif /* SCYTHE_MATRIX_H */
