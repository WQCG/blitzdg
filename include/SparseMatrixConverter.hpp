// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include "SparseTriplet.hpp"
#include <blitz/array.h>
#include <suitesparse/umfpack.h>

using namespace blitz;

namespace blitzdg {
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
} // namespace blitzdg
