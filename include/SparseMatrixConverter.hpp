// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file SparseMatrixConverter.hpp
 * @brief Defines the SparseMatrixConverter class converts between blitz++
 * 2D arrays, sparse triplets, sparse compressed column formats, and contiguous C-style (POD) arrays.
 */

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

    void podArrayToFull(const double * Apod, Array<double, 2> & A);

    int getNumNonZeros(const Array<double, 2> & A);

    ~SparseMatrixConverter();
};
