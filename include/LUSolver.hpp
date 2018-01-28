#pragma once
#include <blitz/array.h>
#include <SparseMatrixConverter.hpp>

using namespace blitz;

class LUSolver {
    int N;
    Array<double, 2> * A;

    // Umfpack-specific fields
    int * Ap;
    int * Ai;
    double * Ax;
    int * Map;
    double * null = (double *) NULL ;

    void * Symbolic;
    void * Numeric;

    SparseTriplet Triplet;

    SparseMatrixConverter MatrixConverter;
  
  public:
    LUSolver(Array<double, 2> * const &, SparseMatrixConverter const &);
    
    Array<double, 2> & get_A();

    void factorize();

    void solve(Array<double, 1> const &, Array<double, 1> &);

    ~LUSolver();
};

