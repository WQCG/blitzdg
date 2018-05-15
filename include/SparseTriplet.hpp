// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file SparseTriplet.hpp
 * @brief Defines the SparseTriplet struct for storing sparse matrices in a conceptually simple way.
 */

#pragma once
#include "Types.hpp"

namespace blitzdg {
    struct SparseTriplet {
        index_type nz;
        index_type* row;
        index_type* col;
        real_type* val;
    };
} // namespace blitzdg