// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file SparseTriplet.hpp
 * @brief Defines the SparseTriplet struct for storing sparse matrices in a conceptually simple way.
 */

#pragma once
struct SparseTriplet {
    int nz;
    int *row;
    int *col;
    double *val;
};