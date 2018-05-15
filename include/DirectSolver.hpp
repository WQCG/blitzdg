// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include "SparseMatrixConverter.hpp"
#include "Types.hpp"

namespace blitzdg {
  class DirectSolver {
      index_type N;
      SparseMatrixConverter MatrixConverter;

    public:
      DirectSolver(SparseMatrixConverter const &);

      void solve(const matrix_type& A, const matrix_type& B, matrix_type& X);
      
      ~DirectSolver();
  };
} // namespace blitzdg

