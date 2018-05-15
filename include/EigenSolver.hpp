// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file EigenSolver.hpp
 * @brief Defines the EigenSolver class that implements LAPACK solution
 * routine DSYEVD. For symmetric matrices only. Documentation at http://www.netlib.org/lapack/explore-3.1.1-html/dsyevd.f.html.
 */

#pragma once
#include <blitz/array.h>
#include <SparseMatrixConverter.hpp>

using namespace blitz;

class EigenSolver {
    int N;
    SparseMatrixConverter MatrixConverter;

  public:
    EigenSolver(SparseMatrixConverter const &);

    void solve(const Array<double,2> & A, Array<double,1> & eigenvalues, Array<double, 2> & eigenvectors);
    
    ~EigenSolver();
};

