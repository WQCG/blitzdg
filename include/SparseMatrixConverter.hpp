#pragma once
#include <blitz/array.h>
#include <SparseTriplet.hpp>
#include <suitesparse/umfpack.h>

using namespace blitz;

class SparseMatrixConverter {

  public:
    SparseMatrixConverter();

    void fullToSparseTriplet(const Array<double, 2> & A, SparseTriplet & triplet);

    void sparseTripletToCompressedColumn(const int numRows, const int numCols, const SparseTriplet & triplet, int * Aptr, int * Aind, double * Avalues);

    void fullToCompressedColumn(const Array<double, 2> & A, int * Aptr, int * Aind, double * Avalues);

    void fullToPodArray(const Array<double, 2> & A, double * Apod);

    int getNumNonZeros(const Array<double, 2> & A);

    ~SparseMatrixConverter();
};
