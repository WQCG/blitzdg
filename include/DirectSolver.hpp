// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file DirectSolver.hpp
 * @brief Defines the DirectSolver class that implements LAPACK direct solution
 * routine DSGESV. Documentation at http://www.netlib.org/lapack/lapack-3.1.1/html/dsgesv.f.html.
 */

#pragma once
#include <blitz/array.h>
#include <SparseMatrixConverter.hpp>

using namespace blitz;

class DirectSolver {
    int N;
    SparseMatrixConverter MatrixConverter;

  public:
    DirectSolver(SparseMatrixConverter const &);

    void solve(const Array<double,2> & A, const Array<double, 2> & B, Array<double, 2> & X);
    
    ~DirectSolver();
};

