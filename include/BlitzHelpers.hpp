// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file BlitzHelpers.hpp
 * @brief Defines a set of utility functions for working with blitz arrays
 * and their matrix and vector aliases.
 */
#pragma once
#include "Traits.hpp"
#include "Types.hpp"
#include <blitz/array.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace blitzdg {
    /**
     * Used for constructing a matrix with column-major ordering.
     * 
     * The following example shows how to build a 5 x 5
     * real matrix with column-major ordering.
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * real_matrix_type mat(5, 5, ColumnMajorOrder());
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    using ColumnMajorOrder = blitz::ColumnMajorArray<2>;

    /**
     * Used for constructing a matrix with row-major ordering.
     * 
     * The following example shows how to build a 5 x 5
     * real matrix with row-major ordering.
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * real_matrix_type mat(5, 5, RowMajorOrder());
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * @note: Row-major ordering is the default ordering
     * and does not need to be explicitly specified.
     */
    using RowMajorOrder = blitz::GeneralArrayStorage<2>;

    /**
     * Returns true if the matrix use row-major storage.
     */
    template <typename T>
    bool isRowMajor(const matrix_type<T>& mat) {
        return mat.isMajorRank(0);
    }

    /**
     * Returns true if the matrix uses column-major storage.
     */
    template <typename T>
    bool isColumnMajor(const matrix_type<T>& mat) {
        return mat.isMajorRank(1);
    }

    /**
     * Returns the number of nonzero elements in the dense matrix.
     * @param[in] mat The dense matrix.
     * @param[in] dropTol The drop tolerance. Defaults to 0.
     * @note An element \f$a_{ij}\f$ of the input matrix is considered
     * nonzero if \f$|a_{ij}| > \mathrm{dropTol}\f$. 
     */
    template <typename T>
    index_type countNonzeros(const matrix_type<T>& mat, real_type dropTol = real_type(0)) {
        std::size_t nnz = 0;
        for (typename matrix_type<T>::const_iterator itr = mat.begin(); itr != mat.end(); ++itr) {
            if (std::abs(*itr) > dropTol)
                ++nnz;
        }
        if (static_cast<index_type>(nnz) > std::numeric_limits<index_type>::max())
            throw std::runtime_error("countNonzeros: number of nonzero elements exceeds maximum allowable");
        return static_cast<index_type>(nnz);
    }

    /**
     * Vectorizes a dense matrix by writing it to an array.
     * @param[in] mat The Dense matrix.
     * @param[in] byRows If true, the matrix is copied to the array by rows. Defaults to true.
     * @param[out] arrItr An output iterator to an array.
     * @note The array must have space allocated for at 
     * least m*n elements, where m is the the number of rows of mat and n is 
     * the number of columns. 
     */
    template <typename T, typename OutputItr>
    void reshapeMatTo1D(const matrix_type<T>& mat, OutputItr arrItr, bool byRows = true) {
        // Fail at compile time if the type T is not
        // the same as the value type of OutputItr.
        static_assert(iteratorValueTypeSameAs<OutputItr, T>(),
        "Matrix value type differs from array value type");
        if (byRows) {
            for (index_type i = 0; i < mat.rows(); ++i) {
                for (index_type j = 0; j < mat.cols(); ++j)
                    *arrItr++ = mat(i, j);
            }
        }
        else {
            for (index_type j = 0; j < mat.cols(); ++j) {
                for (index_type i = 0; i < mat.rows(); ++i)
                    *arrItr++ = mat(i, j);
            }
        }
    }

    /**
     * Reshapes an array to a dense matrix.
     * @param[in] arrItr An input iterator to the array.
     * @param[in] byRows If true, the array is stored in mat rowwise. Defaults to true.
     * @param[out] mat The dense matrix.
     * @note The array must have space allocated for at 
     * least m*n elements, where m is the the number of rows of mat and n is 
     * the number of columns. 
     */
    template <typename T, typename InputItr>
    void reshape1DToMat(InputItr arrItr, matrix_type<T>& mat, bool byRows = true) {
        // Fail at compile time if the type T is not
        // the same as the value type of InputItr.
        static_assert(iteratorValueTypeSameAs<InputItr, T>(),
        "Matrix value type differs from array value type");
        if (byRows) {
            for (index_type i = 0; i < mat.rows(); ++i) {
                for (index_type j = 0; j < mat.cols(); ++j)
                    mat(i, j) = *arrItr++;
            }
        }
        else {
            for (index_type j = 0; j < mat.cols(); ++j) {
                for (index_type i = 0; i < mat.rows(); ++i)
                   mat(i, j) = *arrItr++;
            }
        }
    }

	/**
	 * Evaluate vector at an array of indices.
	 * @param[in] vec The input vector to be evaluated at a list of indices (or map).
	 * @param[in] map The map or list of indices.
	 * @param[out] out The output vector containing the result of applying the map.
	 */
	template <typename T, typename U>
	void applyIndexMap(const vector_type<T>& vec, const vector_type<U>& map, vector_type<T>& out) {
		static_assert(isIntegral<U>(), "map value_type is not integral");
        for (index_type k = 0; k < map.length(0); ++k)
            out(k) = vec(map(k));
    }

	/**
	 * Convert a vector to a blitz matrix.
	 * @param[in] vec The vector.
	 * @param[out] mat The output dense matrix.
	 * @param[in] byRows whether to convert row-wise (default) or column-wise.
	 */
	template <typename T>
	void vectorToFull(const vector_type<T>& vec, matrix_type<T>& mat, bool byRows = true) {
		reshape1DToMat(vec.begin(), mat, byRows);
	}

	/**
     * Convert a blitz matrix to a vector.
	 * @param[in] mat The dense matrix.
	 * @param[out] vec The output vector.
	 * @param[in] byRows whether to convert row-wise (default) or column-wise.
     */
	template <typename T>
	void fullToVector(const matrix_type<T>& mat, vector_type<T>& vec, bool byRows = true) {
		reshapeMatTo1D(mat, vec.begin(), byRows);
	}

    /**
     *  Computes Kronecker product of matrices A and B.
     *  @param[in] A The left matrix.
     *  @param[in] B The right matrix.
     *  @returns C = kron(A,B)
     */
    template <typename T>
    matrix_type<T> kron(const matrix_type<T>& A, const matrix_type<T>&B) {
        using range = blitz::Range;
        const index_type Br = B.rows(), Bc = B.cols();
        const index_type Ar = A.rows(), Ac = A.cols();
        matrix_type<T> C(Ar*Br, Ac*Bc);

        for (index_type i = 0; i < Ar; ++i) {
            for (index_type j = 0; j < Ac; ++j) {
                C(range(j * Br, (j + 1) * Br - 1), 
                range(i * Bc, (i + 1) * Bc - 1)) = A(i,j) * B;
            }
        }

        return C;
    }
} // namespace blitzdg