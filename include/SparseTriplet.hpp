// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

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