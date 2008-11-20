
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
 *  scythestat/datablock.h
 */

/*! \file datablock.h
 * \brief Definitions of internal Matrix data management classes
 *
 * DataBlock and DataBlockReference objects provide the data half of
 * the data/view container model used by Scythe's matrices.  A
 * DataBlock contains a data array of a given type, some information
 * about the DataBlock, and a reference count.  Matrix objects
 * provide views to the DataBlock, thus allowing us to provide
 * Matrix objects that reference subsets of other matrices.  When no
 * matrices remain that reference the DataBlock the reference count
 * falls to zero and the block is automatically deallocated.
 *
 * DataBlock uses a simple doubling/halving memory allocation scheme
 * but this may change in later releases.
 *
 * The DataBlock classes are used exclusively within the library and
 * do not constitute a part of Scythe's public interface.
 *
 * Based on code in Blitz++ (http://www.oonumerics.org/blitz/) by
 * Todd Veldhuizen <tveldhui@oonumerics.org>.  Blitz++ is
 * distributed under the terms of the GNU GPL.
 */

#ifndef SCYTHE_DATABLOCK_H
#define SCYTHE_DATABLOCK_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "error.h"
#else
#include "scythestat/error.h"
#endif

namespace scythe {
	/* Convenience typedefs */
  namespace { // local to this file
	  typedef unsigned int uint;
  }

  /*!  \brief Handles Matrix data internals.
	 * 
   * Handles data allocation, reallocation, and deletion of blocks of
   * elements; the actual data Matrix objects point to.  Keeps a
   * reference count.
	 */
  template <typename T_type>
  class DataBlock { 
    public:
      /**** CONSTRUCTORS ****/
			
      /*
       * Create an empty data block.
       */
			
			DataBlock ()
				:	data_ (0),
					size_ (0),
					refs_ (0)
			{}

      /* 
       * Create a block of a given size.
       */
			explicit
			DataBlock (uint size)
				:	data_ (0),
					size_ (0),
					refs_ (0)
			{
				resize(size);
				SCYTHE_DEBUG_MSG("Constructed new " << size << "(" << size_
						<< ") DataBlock at address " << data_);
			}

      /*
       * Create an exact copy of another data block.
       */
			DataBlock (const DataBlock<T_type>& b)
				:	data_ (b.data_),
					size_ (b.size_),
					refs_ (b.refs_)
			{}

			/**** DESTRUCTOR ****/

			~DataBlock ()
			{
				SCYTHE_DEBUG_MSG("Destructing block at " << data_);
				deallocate();
			}

			/**** REFERENCE COUNTING ****/

			inline uint addReference ()
			{
        SCYTHE_DEBUG_MSG("Added reference to DataBlock at address "
            << data_);
				return ++refs_;
			}

			inline uint removeReference ()
			{
        SCYTHE_DEBUG_MSG("Removed reference to DataBlock at address "
           << data_);
				return --refs_ ;
			}

			inline uint references ()
			{
				return refs_;
			}

			/**** ACCESSORS ****/

			inline T_type* data()
			{
				return data_;
			}

			inline const T_type* data() const
			{
				return data_;
			}

			inline uint size () const
			{
				return size_;
			}

		protected:
			/**** (DE)ALLOCATION AND RESIZING ****/
			
			/* Allocate data given the current block size. */
			inline void allocate (uint size)
			{
				/* TODO Think about cache boundary allocations for big blocks
				 * see blitz++ */

				if (data_ != 0) // Get rid of previous allocation if it exists
					deallocate();

				data_ = new (std::nothrow) T_type[size];
				
				SCYTHE_CHECK_10(data_ == 0, scythe_alloc_error,
						"Failure allocating DataBlock of size " << size);
			}

			/* Deallocate a block's data */
			inline void deallocate ()
			{
				SCYTHE_DEBUG_MSG("  Deallocating DataBlock of size " << size_
						<< " at address " << data_);
				delete[] data_;
				data_ = 0;
			}

		public:
			/* TODO At the moment, references call this method directly.  Not
			 * sure if this is the best interface choice. */
			/* Resize a block. */
			void resize (uint newsize)
			{
				if (newsize > size_)
					grow(newsize);
				else if (newsize < size_ / 4)
					shrink();
			}

		protected:
			/* Make a block larger. Expects to be called by resize and does
			 * not reset the size_ variable. */
			inline void grow (uint newsize)
			{
				size_ = size_ ? size_ : 1; // make sure not zero

				/* TODO Can we speed this up?  In 20 iters we're at
				 * 1048576 elems doing the math might be more costly...
				 */
				while (size_ < newsize)
					size_ <<= 1;

				allocate(size_);
			}

			/* Make a block smaller. Expects to be called by resize */
			inline void shrink ()
			{
				size_ >>= 1;
				allocate(size_);
			}

		private:
			/**** INSTANCE VARIABLES ****/
			T_type *data_;   // The data array
			uint size_;  // The number of elements in the block
			uint refs_;  // The number of views looking at this block
	}; // end class DataBlock

	/*! \brief Null data block object.
   *
   * A nice little way to represent empty data blocks.
   */
	template <class T_type>
	class NullDataBlock : public DataBlock<T_type>
	{
		typedef DataBlock<T_type> T_base;
		public:
			
			NullDataBlock ()
				: DataBlock<T_type> ()
			{
        // never want to deallocate (or resize) this one
				T_base::addReference(); 
				SCYTHE_DEBUG_MSG("Constructed NULL datablock");
			}

			~NullDataBlock ()
			{}

	}; // end class NullDataBlock


  /*! 
   * \brief Handle to DataBlock objects.
   *
	 * Matrices inherit from this object.  It provides a handle into
	 * DataBlock objects and automates cleanup when the referenced
	 * object runs out of referants.
	 */
	template <class T_type>
	class DataBlockReference {
		public:
			/**** CONSTRUCTORS ****/

			/* Default constructor: points the object at a static null block
			 */
			DataBlockReference ()
				:	data_ (0),
					block_ (&nullBlock_)
			{
				block_->addReference();
			}

			/* New block constructor: creates a new underlying block of a
			 * given size and points at it. */
			explicit
			DataBlockReference (uint size)
				:	data_ (0),
					block_ (0)
			{
				block_ = new (std::nothrow) DataBlock<T_type> (size);
				SCYTHE_CHECK_10 (block_ == 0, scythe_alloc_error,
						"Could not allocate DataBlock object");
				
				data_ = block_->data();
				block_->addReference();
			}

			/* Refrence to an existing block constructor: points to an
			 * offset within an existing block. */
			DataBlockReference (const DataBlockReference<T_type>& reference,
					uint offset = 0)
				:	data_ (reference.data_ + offset),
					block_ (reference.block_)
			{
				block_->addReference();
			}
			
			/**** DESTRUCTOR ****/
			/* Automates removal of underlying block objects when refcount
			 * hits nil.
			 */
			virtual ~DataBlockReference ()
			{
				withdrawReference();
			}

		protected:

			/**** MEMBERS CALLED BY DERIVED CLASS ****/
			void referenceOther (const DataBlockReference<T_type>& ref,
					uint offset = 0)
			{
				withdrawReference ();
				block_ = ref.block_;
				block_->addReference();
				data_ = ref.data_ + offset;
			}

			void referenceNew (uint size)
			{
				/* If we are the only referent to this data block, resize it.
				 * Otherwise, shift the reference to point to a newly
				 * constructed block.
				 */
				if (block_->references() == 1) {
					block_->resize(size);
          data_ = block_->data(); // This is a pretty good indication
          // that the interface and implementation are too tightly
          // coupled for resizing.
				} else {
					withdrawReference();
					block_ = 0;
					block_ = new (std::nothrow) DataBlock<T_type> (size);
					SCYTHE_CHECK_10(block_ == 0, scythe_alloc_error,
							"Could not allocate new data block");
					data_ = block_->data();
					block_->addReference();
				}
			}

		private:
			/**** INTERNAL MEMBERS ****/
			void withdrawReference ()
			{
				if (block_->removeReference() == 0
						&& block_ != &nullBlock_)
					delete block_;
			}

			void referenceNull ()
			{
				withdrawReference();
				block_ = &nullBlock_;
				block_->addReference();
				data_ = 0;
			}


		/**** INSTANCE VARIABLES ****/
		protected:
			T_type* data_;  // Pointer to the underlying data (offset)
		
		private:
			DataBlock<T_type>* block_;
			static NullDataBlock<T_type> nullBlock_;

	}; // end class DataBlockReference

	/* Instantiation of the static null memory block */
	template <typename T>
	NullDataBlock<T> DataBlockReference<T>::nullBlock_;

} // end namespace scythe

#endif /* SCYTHE_DATABLOCK_H */
