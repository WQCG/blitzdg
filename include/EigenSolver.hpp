// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

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

