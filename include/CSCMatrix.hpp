// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file CSCMatrix.hpp
 * Defines a compressed sparse column matrix class.
 */
#pragma once
#include "SparseTriplet.hpp"
#include "Types.hpp"
#include <suitesparse/cs.h>
#include <ostream>
#include <memory>

namespace blitzdg {
    /**
     * Implements a compressed sparse column (CSC) matrix class
     * that wraps a CXSparse cs_di object.
     */
	class CSCMat {
        /**
         * Custom deleter used by std::unique_ptr for
         * objects of type cs_di created by cs_di_spalloc.
         */
		struct deleter {
			void operator()(cs_di* ptr) const {
				cs_di_spfree(ptr);
				ptr = nullptr;
			}
		};
	public:
		using cs_di_smart_ptr = std::unique_ptr<cs_di, deleter>; /**< smart pointer for cs_di objects */
		
        /**
         * Default constructor that creates an empty 0 by 0 CSC matrix.
         */
		CSCMat()
			: CSCMat(0, 0, 0)
		{}
		
        /**
         * Constructor that creates a rows x cols CSC matrix with nnz nonzero elements.
         * 
         * On creation, all the elements of the matrix are set to zero. Each
         * element should be assigned a nonzero value prior to using this matrix, 
         * unless storing explicit zeroes is acceptable. Note that when storing 
         * explicit zeroes, nnz() will overestimate the number of nonzero elements.
         * @param[in] rows The number of rows.
         * @param[in] cols The number of columns.
         * @param[in] nnz The number of nonzero elements.
         */
		CSCMat(index_type rows, index_type cols, index_type nnz);

        /**
         * Constructor that creates a CSC matrix from a sparse triplet.
         * @param[in] triplet The sparse triplet matrix.
         */
        explicit CSCMat(const SparseTriplet& triplet);

        /**
         * Constructor that creates a CSC matrix from a dense matrix.
         * @param[in] mat The dense matrix.
         * @param[in] dropTol Drop tolerance for identifying nonzero elements in mat. Defaults to 0.
         * @note An element \f$a_{ij}\f$ of mat is considered nonzero if
         * \f$|a_{ij}| > \mathrm{dropTol}\f$.
         */
        explicit CSCMat(const matrix_type& mat, real_type dropTol = real_type(0));
		
        /**
         * Constructor that creates a matrix from a cs_di object
         * that is owned by a cs_di_smart_ptr. 
         * 
         * This constructor implies a transfer of ownership from
         * the smart pointer mat to the CSCMat object. Note that
         * the cs_di object must have its values attribute set to
         * true, otherwise an std::runtime_error exception is thrown.
         * @param[in] mat Smart pointer to cs_di object.
         */
		explicit CSCMat(cs_di_smart_ptr mat);
		
        /**
         * Copy constructor. Does a deep copy.
         */
		CSCMat(const CSCMat& other);

        /**
         * Assignment operator. Does a deep copy.
         */
		CSCMat& operator=(CSCMat other);

        /**
         * Move constructor.
         */
		CSCMat(CSCMat&& other) = default;

        /**
         * Move assignment operator.
         */
		CSCMat& operator=(CSCMat&& other) = default;
		
        /**
         * Returns true if the matrix is empty, i.e.,
         * has a zero dimension.
         */
		bool empty() const {
			return (mat_->m == 0 || mat_->n == 0);
		}
		
        /**
         * Returns the number of rows.
         */
		index_type rows() const {
			return mat_->m;
		}
		
        /**
         * Returns the number of columns.
         */
		index_type cols() const {
			return mat_->n;
		}
		
        /**
         * Returns the number of nonzero elements.
         */
		index_type nnz() const {
			return mat_->nzmax;
		}
		
        /**
         * Provides read-write access to the kth
         * row index.
         */
		index_type& rowInds(index_type k) {
			return mat_->i[k];
		}
		
        /**
         * Provides read-write access to the kth
         * column pointer index.
         */
		index_type& colPtrs(index_type k) {
			return mat_->p[k];
		}
		
        /**
         * Provides read-write access to the kth
         * nonzero element.
         */
		real_type& elems(index_type k) {
			return mat_->x[k];
		}
		
        /**
         * Provides write access to the kth
         * row index.
         */
		index_type rowInds(index_type k) const {
			return mat_->i[k];
		}
		
        /**
         * Provides write access to the kth
         * column pointer index.
         */
		index_type colPtrs(index_type k) const {
			return mat_->p[k];
		}
		
        /**
         * Provides write access to the kth
         * nonzero value.
         */
		real_type elems(index_type k) const {
			return mat_->x[k];
		}
		
        /**
         * Provides access to the underlying raw pointer
         * that points to the managed cs_di object.
         * This pointer is used to interface with the 
         * functions in the CXSparse library.
         */
		const cs_di* matPtr() const {
			return mat_.get();
		}
		
        /**
         * Returns a const pointer to the  
         * array of row indices.
         */
		const index_type* rowInds() const {
			return mat_->i;
		}
		
        /**
         * Returns a const pointer to the 
         * array that stores the locations of 
         * the start of each column.
         */
		const index_type* colPtrs() const {
			return mat_->p;
		}
		
        /**
         * Returns a const pointer to the array
         * of values.
         */
		const real_type* elems() const {
			return mat_->x;
		}

        /**
         * Removes any elements with duplicate row and
         * column indices by summing their values.
         */
        void removeDuplicates();

        /**
         * Removes any elements that are <= dropTol
         * in absolute value.
         * @param[in] dropTol Drop tolerance. Defaults to 0.
         */
        void prune(real_type dropTol = real_type(0));

        /**
         * Overwrites this matrix with its transpose.
         */
        void transpose();
		
        /**
         * Swaps the contents of two CSCMat objects.
         */
		friend void swap(CSCMat& lhs, CSCMat& rhs);	
	private:
        // The cs_di object is wrapped in a unique_ptr
        // that manages the lifetime of the object.
		cs_di_smart_ptr mat_;
	};

    /**
     * Writes the CSC matrix to the output stream 
     * in coordinate format (i, j, elem).
     */
    std::ostream& operator<<(std::ostream& strm, const CSCMat& mat);
	
    /**
     * Returns the transpose of a CSC matrix.
     * @param[in] mat The matrix to be transposed.
     */
	CSCMat transpose(const CSCMat& mat);
	
    /**
     * Returns the matrix-matrix product of two 
     * CSC matrices: lhs * rhs.
     * @param[in] lhs The left-hand side matrix.
     * @param[in] rhs The right-hand side matrix.
     */
	CSCMat multiply(const CSCMat& lhs, const CSCMat& rhs);
} // namespace blitzdg
