// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file DirectSolver.hpp
 * @brief Defines the DirectSolver class that implements LAPACK direct solution
 * routine DSGESV. Documentation at http://www.netlib.org/lapack/lapack-3.1.1/html/dsgesv.f.html.
 */
#pragma once
#include "Types.hpp"

namespace blitzdg {
    /**
     * Implements the direct solution of a linear system via the LAPACK routine DSGESV,
     * which employs LU factorization with row pivoting and iterative refinement.
     */
    class DirectSolver {  
    public:
        /**
         * Solves the linear system \f$AX = B\f$.
         * @param[in] A The \f$n\times n\f$ coefficient matrix.
         * @param[in] B The \f$n\times k\f$ right-hand side matrix (one per column).
         * @param[out] X The \f$n\times k\f$ solution matrix.
         * @note We assume that the matrix \f$A\f$ uses the default rowwise storage order
         * of blitz++ 2D arrays.
         */
        void solve(const real_matrix_type& A, const real_matrix_type& B, real_matrix_type& X) const;
  };
} // namespace blitzdg

