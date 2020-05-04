// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file DenseCholeskyFactorizer.hpp
 * @brief Defines the DenseCholeskyFactorizer class that computes the Cholesky factorization 
 * of a real symmetric positive semidefinite matrix A.
 * Documentation at: http://www.netlib.org/lapack/explore-html/d1/d7a/group__double_p_ocomputational_ga2f55f604a6003d03b5cd4a0adcfb74d6.html
 */
#pragma once
#include "Types.hpp"

namespace blitzdg {
    /**
     * Implements the explicit Cholesky factorization of a dense matrix via the LAPACK routine
     * DPOTRF.
     */
    class DenseCholeskyFactorizer {  
    public:
        /**
         * Inverts the dense matrix \f$A\f$.s
         * @param[in] A The \f$n\times n\f$ coefficient matrix.
         * @param[out] R The \f$n\times n\f$ upper triangular matrix such that A = R^T * R
         * @note We assume that the matrix \f$A\f$ uses the default rowwise storage order
         * of blitz++ 2D arrays.
         */
        void computeCholesky(const real_matrix_type& A, real_matrix_type& R) const;
  };
} // namespace blitzdg