// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file DenseMatrixHelpers.hpp
 * @brief Defines a set of utility functions for working with dense matrices.
 */
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
        if (nnz > std::numeric_limits<index_type>::max())
            throw std::runtime_error("countNonzeros: number of nonzero elements exceeds maximum allowable");
        return static_cast<index_type>(nnz);
    }

    /**
     * Vectorizes a dense matrix. 
     * @param[in] mat The Dense matrix.
     * @param[in] byRows If true, the rows of mat are stored contiguously. Defaults to true.
     * @param[out] arr Pointer to array that stores the vectorized matrix.
     * @note The array pointed to by arr must have space allocated for at 
     * least m*n elements, where m is the the number of rows of mat and n is 
     * the number of columns. 
     */
    template <typename T>
    void fullToPodArray(const matrix_type<T>& mat, T* arr, bool byRows = true) {
        std::size_t nnz = 0;
        if (byRows) {
            for (index_type i = 0; i < mat.rows(); ++i) {
                for (index_type j = 0; j < mat.cols(); ++j)
                    arr[nnz++] = mat(i, j);
            }
        }
        else {
            for (index_type j = 0; j < mat.cols(); ++j) {
                for (index_type i = 0; i < mat.rows(); ++i)
                    arr[nnz++] = mat(i, j);
            }
        }
    }

    /**
     * Reshapes an array to a dense matrix.
     * @param[in] arr Pointer to the array.
     * @param[in] byRows If true, arr is stored in mat rowwise. Defaults to true.
     * @param[out] mat The dense matrix.
     * @note The array pointed to by arr must have space allocated for at 
     * least m*n elements, where m is the the number of rows of mat and n is 
     * the number of columns. 
     */
    template <typename T>
    void podArrayToFull(const T* arr, matrix_type<T>& mat, bool byRows = true) {
        std::size_t nnz = 0;
        if (byRows) {
            for (index_type i = 0; i < mat.rows(); ++i) {
                for (index_type j = 0; j < mat.cols(); ++j)
                    mat(i, j) = arr[nnz++];
            }
        }
        else {
            for (index_type j = 0; j < mat.cols(); ++j) {
                for (index_type i = 0; i < mat.rows(); ++i)
                   mat(i, j) = arr[nnz++];
            }
        }
    }
} // namespace blitzdg