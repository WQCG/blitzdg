// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include "SparseMatrixConverter.hpp"
#include "Types.hpp"

namespace blitzdg {
  class EigenSolver {
      index_type N;
      SparseMatrixConverter MatrixConverter;

    public:
      EigenSolver(SparseMatrixConverter const &);

      void solve(const matrix_type& A, vector_type& eigenvalues, matrix_type& eigenvectors);
      
      ~EigenSolver();
  };
} // namespace blitzdg

