// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file DenseMatrixHelpers.hpp
 * @brief Defines a set of utility functions for working with dense matrices.
 */
#include "Types.hpp"

namespace blitzdg {
    /**
     * Returns true if the matrix use row-major storage.
     */
    inline bool isRowMajor(const matrix_type& mat) {
        return mat.isMajorRank(1);
    }

    /**
     * Returns true if the matrix uses column-major storage.
     */
    inline bool isColumnMajor(const matrix_type& mat) {
        return mat.isMajorRank(2);
    }

    /**
     * Returns the number of nonzero elements in the dense matrix.
     * @param[in] mat The dense matrix.
     * @param[in] dropTol The drop tolerance. Defaults to 0.
     * @note An element \f$a_{ij}\f$ of the input matrix is considered
     * nonzero if \f$|a_{ij}| > \mathrm{dropTol}\f$. 
     */
    index_type countNonzeros(const matrix_type& mat, real_type dropTol = real_type(0));

    /**
     * Vectorizes a dense matrix. 
     * @param[in] mat The Dense matrix.
     * @param[in] byRows If true, the rows of mat are stored contiguously. Defaults to true.
     * @param[out] arr Pointer to array that stores the vectorized matrix.
     * @note The array pointed to by arr must have space allocated for at 
     * least m*n elements, where m is the the number of rows of mat and n is 
     * the number of columns. 
     */
    void fullToPodArray(const matrix_type& mat, real_type* arr, bool byRows = true);

    /**
     * Reshapes an array to a dense matrix.
     * @param[in] arr Pointer to the array.
     * @param[in] byRows If true, arr is stored in mat rowwise. Defaults to true.
     * @param[out] mat The dense matrix.
     * @note The array pointed to by arr must have space allocated for at 
     * least m*n elements, where m is the the number of rows of mat and n is 
     * the number of columns. 
     */
    void podArrayToFull(const real_type* arr, matrix_type& mat, bool byRows = true);
} // namespace blitzdg