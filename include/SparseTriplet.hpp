// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
struct SparseTriplet {
    int nz;
    int *row;
    int *col;
    double *val;
};