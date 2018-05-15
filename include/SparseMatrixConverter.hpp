// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include "SparseTriplet.hpp"
#include "Types.hpp"

namespace blitzdg {
  class SparseMatrixConverter {

    public:
      SparseMatrixConverter();

      void fullToSparseTriplet(const matrix_type& A, SparseTriplet & triplet);

      void sparseTripletToCompressedColumn(const index_type numRows, const index_type numCols, const SparseTriplet & triplet, index_type * Aptr, index_type * Aind, real_type * Avalues);

      void fullToCompressedColumn(const matrix_type& A, index_type * Aptr, index_type * Aind, real_type * Avalues);

      void fullToPodArray(const matrix_type& A, real_type * Apod);

      void podArrayToFull(const real_type * Apod, matrix_type& A);

      index_type getNumNonZeros(const matrix_type& A);

      ~SparseMatrixConverter();
  };
} // namespace blitzdg
