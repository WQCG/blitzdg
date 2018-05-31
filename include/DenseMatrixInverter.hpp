// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file DirectSolver.hpp
 * @brief Defines the DenseMatrixInverter class that ... does something.
 */
#pragma once
#include "Types.hpp"

namespace blitzdg {
    /**
     * Implements the dense matrix LAPACK routine DSGESV,
     */
    class DenseMatrixInverter {  
    public:
        /**
         * Inverts the dense matrix \f$A\f$.
         * @param[in] A The \f$n\times n\f$ coefficient matrix.
         * @param[out] Ainv The \f$n\times n\f$ inverse of matrix A.
         * @note We assume that the matrix \f$A\f$ uses the default rowwise storage order
         * of blitz++ 2D arrays.
         */
        void computeInverse(const real_matrix_type& A, real_matrix_type& Ainv) const;
  };
} // namespace blitzdg