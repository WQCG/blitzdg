#pragma once
#include <blitz/array.h>
#include <SparseMatrixConverter.hpp>

using namespace blitz;

class EigenSolver {
    int N;
    Array<double, 2> * A;
    SparseMatrixConverter MatrixConverter;

  public:
    EigenSolver(Array<double, 2> * const &, SparseMatrixConverter const &);

    void solve(Array<double,1> & eigenvalues, Array<double, 2> & eigenvectors);
    
    ~EigenSolver();
};

