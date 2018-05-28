// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file EigenSolver.hpp
 * @brief Defines the EigenSolver class that implements LAPACK solution
 * routine DSYEVD. For symmetric matrices only. 
 * Documentation at http://www.netlib.org/lapack/explore-3.1.1-html/dsyevd.f.html.
 */

#pragma once
#include "Types.hpp"

namespace blitzdg {
  class EigenSolver {
  public:
      /**
       * Computes the eigenvalues and eigenvectors of the symmetric matrix \f$A\f$.
       * @param[in] A The \f$n\times n\f$ symmetric coefficient matrix.
       * @param[out] eigenvalues The length \f$n\f$ array of eigenvalues.
       * @param[out] eigenvectors The \f$n\times n\f$ matrix of eigenvectors (one per column).
       * @note We assume that the matrix \f$A\f$ uses the default rowwise storage order
       * of blitz++ 2D arrays.
       */
      void solve(const real_matrix_type& A, real_vector_type& eigenvalues, real_matrix_type& eigenvectors) const;      
  };
} // namespace blitzdg

