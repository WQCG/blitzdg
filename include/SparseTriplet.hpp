// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
namespace blitzdg {
    struct SparseTriplet {
        int nz;
        int *row;
        int *col;
        double *val;
    };
} // namespace blitzdg