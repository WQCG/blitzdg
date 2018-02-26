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

