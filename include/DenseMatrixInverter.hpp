// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file DenseMatrixInverter.hpp
 * @brief Defines the DenseMatrixInverter class that computes the explicit
 * inverse of a dense matrix via the LAPACK routines DGETRF and DGETRI.
 * Documentation at: http://www.netlib.org/lapack/lapack-3.1.1/html/dgetrf.f.html
 * and http://www.netlib.org/lapack/lapack-3.1.1/html/dgetri.f.html.
 */
#pragma once
#include "Types.hpp"

namespace blitzdg {
    /**
     * Implements the explicit inverse of a dense matrix via the LAPACK routines
     * DGETRF and DGETRI. The routine DGETRF computes the LU factorization of the
     * dense matrix and DGETRI computes the inverse from the LU factorization.
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