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
    };
} // namespace blitzdg