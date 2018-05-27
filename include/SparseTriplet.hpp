// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file SparseTriplet.hpp
 * @brief Defines a sparse triplet matrix class.
 */
#pragma once
#include "Types.hpp"
#include <algorithm>
#include <vector>

namespace blitzdg {
	/**
	 * Implements a sparse triplet matrix class that stores a sparse matrix
	 * in terms of three arrays that hold the row indices, column indices,
	 * and values of the nonzero elements.
	 */
	class SparseTriplet {
	public:
		/**
		 * Default constructor that creates an empty 0 by 0 sparse triplet matrix.
		 */
		SparseTriplet()
			: SparseTriplet(0, 0, 0)
		{}
		
		/**
		 * Constructor that creates a rows by cols sparse triplet matrix with 
		 * space allocated for nzmax nonzero elements.
		 * @param[in] rows The number of rows.
		 * @param[in] cols The number of columns.
		 * @param[in] nzmax The amount of storage allocated for nonzero elements.
		 */
		SparseTriplet(index_type rows, index_type cols, index_type nzmax)
			: rows_{ rows }, cols_{ cols }, nnz_{ 0 }, nzmax_{ nzmax },
			row_(nzmax, 0), col_(nzmax, 0), elems_(nzmax, real_type(0))
		{
			if (rows == 0 || cols == 0) {
				rows_ = cols_ = nzmax_ = 0;
				row_.resize(0);
				col_.resize(0);
				elems_.resize(0);
			}
		}
		
		/**
		 * Constructs a rows by cols sparse triplet matrix with nnz nonzero elements
		 * from the input arrays.
		 * @param[in] rows The number of rows.
		 * @param[in] cols The number of columns.
		 * @param[in] nnz The number of nonzero elements.
		 * @param[in] rowPtr Pointer to the array of row indices.
		 * @param[in] colPtr Pointer to the array of columned indices.
		 * @param[in] elemsPtr Pointer to the array of values.
		 * @note The input arrays must have at least nnz elements allocated.
		 */
		SparseTriplet(index_type rows, index_type cols, index_type nnz,
			const index_type* rowPtr, const index_type* colPtr, const real_type* elemsPtr)
			: rows_{ rows }, cols_{ cols }, nnz_{ nnz }, nzmax_{ nnz },
			row_(rowPtr, rowPtr + nnz), col_(colPtr, colPtr + nnz),
			elems_(elemsPtr, elemsPtr + nnz)
		{}
		
		/**
		 * Constructs a sparse triplet matrix with nnz nonzero elements
		 * from the input arrays. The dimensions are determined from the
		 * maximum row and column indices.
		 * @param[in] nnz The number of nonzero elements.
		 * @param[in] rowPtr Pointer to the array of row indices.
		 * @param[in] colPtr Pointer to the array of columned indices.
		 * @param[in] elemsPtr Pointer to the array of values.
		 * @note The input arrays must have at least nnz elements allocated.
		 */
		SparseTriplet(index_type nnz, const index_type* rowPtr, 
			const index_type* colPtr, const real_type* elemsPtr)
			: SparseTriplet(*std::max_element(rowPtr, rowPtr + nnz) + 1,
			*std::max_element(colPtr, colPtr + nnz) + 1, 
			nnz, rowPtr, colPtr, elemsPtr)
		{}
		
		/**
		 * Constructor that creates a sparse triplet matrix from a dense matrix.
		 * @param[in] mat The dense matrix.
		 * @param[in] dropTol Drop tolerance for identifying nonzero elements in mat. Defaults to 0.
		 * @note An element \f$a_{ij}\f$ of mat is considered nonzero if
         * \f$|a_{ij}| > \mathrm{dropTol}\f$.
		 */
		explicit SparseTriplet(const matrix_type& mat, real_type dropTol = real_type(0));
		
		/**
		 * Copy constructor.
		 */
		SparseTriplet(const SparseTriplet& other) = default;
		
		/**
		 * Copy assignment operator.
		 */
		SparseTriplet& operator=(const SparseTriplet& other) = default;
		
		/**
		 * Move constructor.
		 */
		SparseTriplet(SparseTriplet&& other) = default;
		
		/**
		 * Move assignment operator.
		 */
		SparseTriplet& operator=(SparseTriplet&& other) = default;
		
		/**
         * Returns true if the matrix is empty, i.e.,
         * has a zero dimension.
         */
		bool empty() const {
			return (rows_ == 0 || cols_ == 0);
		}
		
		/**
		 * Returns the number of rows.
		 */
		index_type rows() const {
			return rows_;
		}
		
		/**
		 * Returns the number of columns.
		 */
		index_type cols() const {
			return cols_;
		}
		
		/**
		 * Returns the number of nonzero elements.
		 */
		index_type nnz() const {
			return nnz_;
		}
		
		/**
		 * Returns the maximum number of nonzero elements
		 * currently allocated.
		 */
		index_type nzmax() const {
			return nzmax_;
		}
		
		/**
		 * Inserts an element into the sparse triplet matrix.
		 * If the matrix is full, then additional space is allocated.
		 * @param[in] row The row index.
		 * @param[in] col The column index.
		 * @param[in] elem The value.
		 */
		void insert(index_type row, index_type col, real_type elem) {
			if (nnz_ >= nzmax_)
				grow(newSize());
			row_[nnz_] = row;
			col_[nnz_] = col;
			elems_[nnz_++] = elem;
		}
		
		/**
		 * Increases the allocated storage to nzmax.
		 * Does nothing if nzmax is less than the current
		 * amount of allocated storage.
		 * @param[in] nzmax The new amount of allocated storage.
		 */
		void expand(index_type nzmax) {
			grow(nzmax);
		}
		
		/**
		 * Removes all the current elements from the matrix.
		 * Does not reallocate any memory.
		 */
		void clear() {
			std::fill(row_.begin(), row_.end(), 0);
			std::fill(col_.begin(), col_.end(), 0);
			std::fill(elems_.begin(), elems_.end(), real_type(0));
			nnz_ = 0;
		}
		
		/**
		 * Provides read-write access to the
		 * kth row index.
		 */
		index_type& row(index_type k) {
			return row_[k];
		}
		
		/**
		 * Provides read access to the kth row index.
		 */
		index_type row(index_type k) const {
			return row_[k];
		}
		
		/**
		 * Provides read-write access to the
		 * kth column index.
		 */
		index_type& col(index_type k) {
			return col_[k];
		}
		
		/**
		 * Provides read access to the kth column index.
		 */
		index_type col(index_type k) const {
			return col_[k];
		}
		
		/**
		 * Provides read-write access to the
		 * kth element.
		 */ 
		real_type& elem(index_type k) {
			return elems_[k];
		}
		
		/**
		 * Provides read access to the kth element.
		 */
		real_type elem(index_type k) const {
			return elems_[k];
		}
		
		/**
		 * Returns a pointer to the array storing the row indices.
		 */
		const index_type* rowPtr() const {
			return row_.data();
		}
		
		/**
		 * Returns a pointer to the array storing the column indices.
		 */
		const index_type* colPtr() const {
			return col_.data();
		}
		
		/**
		 * Returns a pointer to the array storing the element values.
		 */
		const real_type* elemsPtr() const {
			return elems_.data();
		}

		/**
		 * Swaps the contents of two sparse triplet matrices.
		 */
		friend void swap(SparseTriplet& lhs, SparseTriplet& rhs);
	private:
		index_type rows_;
		index_type cols_;
		index_type nnz_;
		index_type nzmax_;
		std::vector<index_type> row_;
		std::vector<index_type> col_;
		std::vector<real_type> elems_;
		
		/**
		 * Returns the new value of nzmax when automatically
		 * expanding the amount of storage.
		 */
		index_type newSize() const;
		
		/**
		 * Expands the size of the internal arrays to nznew
		 * if nznew > nzmax. Otherwise does nothing.
		 */
		void grow(index_type nznew);
	};

	/**
	 * Writes a sparse triplet matrix to the output stream
	 * in coordinate format (i, j, elem).
	 */
	std::ostream& operator<<(std::ostream& strm, const SparseTriplet& mat);
} // namespace blitzdg