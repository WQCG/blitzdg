// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file EigenSolver.hpp
 * @brief Defines the EigenSolver class that implements LAPACK solution
 * routine DSYEVD. For symmetric matrices only. Documentation at http://www.netlib.org/lapack/explore-3.1.1-html/dsyevd.f.html.
 */

#pragma once
#include "SparseMatrixConverter.hpp"
#include "Types.hpp"

namespace blitzdg {
  class EigenSolver {
      SparseMatrixConverter MatrixConverter;

  public:
      EigenSolver()
        : MatrixConverter()
      {}

      void solve(const matrix_type& A, vector_type& eigenvalues, matrix_type& eigenvectors) const;      
  };
} // namespace blitzdg

