// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file EigenSolver.hpp
 * @brief Defines the EigenSolver class that implements LAPACK solution
 * routine DSYEVD. For symmetric matrices only. 
 * Documentation at http://www.netlib.org/lapack/explore-3.1.1-html/dsyevd.f.html.
 */

#pragma once
#include "SparseMatrixConverter.hpp"
#include "Types.hpp"

namespace blitzdg {
  class EigenSolver {
      SparseMatrixConverter MatrixConverter;

  public:
      /**
       * Default constructor.
       */
      EigenSolver()
        : MatrixConverter()
      {}

      /**
       * Computes the eigenvalues and eigenvectors of the symmetric matrix \f$A\f$.
       * @param[in] A The \f$n\times n\f$ symmetric coefficient matrix.
       * @param[out] eigenvalues The length \f$n\f$ array of eigenvalues.
       * @param[out] eigenvectors The \f$n\times n\f$ matrix of eigenvectors (one per column).
       * @note We assume that the matrix \f$A\f$ uses the default rowwise storage order
       * of blitz++ 2D arrays.
       */
      void solve(const matrix_type& A, vector_type& eigenvalues, matrix_type& eigenvectors) const;      
  };
} // namespace blitzdg

