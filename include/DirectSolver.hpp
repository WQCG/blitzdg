// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file DirectSolver.hpp
 * @brief Defines the DirectSolver class that implements LAPACK direct solution
 * routine DSGESV. Documentation at http://www.netlib.org/lapack/lapack-3.1.1/html/dsgesv.f.html.
 */

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

