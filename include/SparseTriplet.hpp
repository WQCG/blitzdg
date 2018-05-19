// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file SparseTriplet.hpp
 * @brief Defines the SparseTriplet structure for storing sparse matrices in a conceptually simple way.
 */
#pragma once
#include "Types.hpp"

namespace blitzdg {
    struct SparseTriplet {
        index_type nz;   /**< Number of nonzero entries. */
        index_type* row; /**< Pointer to array of row indices. */
        index_type* col; /**< Pointer to array of column indices. */
        real_type* val;  /**< Pointer to array of elements. */

        /**
         * Constructor.
         * Creates a sparse triplet matrix with space for nz_ elements.
         * @param nz_ Number of nonzero elements.
         * @note If nz_ <= 0, then the empty matrix is created. The number of 
         * nonzeros is set to zero and the pointers are set to null.
         */
        explicit SparseTriplet(index_type nz_ = 0)
            : nz{ nz_ > 0 ? nz_ : 0 }, 
            row{ nz_ > 0 ? new index_type[nz_] : nullptr }, 
            col{ nz_ > 0 ? new index_type[nz_] : nullptr }, 
            val{ nz_ > 0 ? new real_type[nz_]{0} : nullptr }
        {}

        /**
         * Destructor.
         */
        ~SparseTriplet() {
            delete[] row; row = nullptr;
            delete[] col; col = nullptr;
            delete[] val; val = nullptr;
        }
    };
} // namespace blitzdg