// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file VandermondeBuilders.hpp
 * @brief Defines a class for computating the generalized Vandermonde, its gradient,
 * and its inverse for 1D orthonormal polynomials.
 */
#pragma once
#include "DirectSolver.hpp"
#include "JacobiBuilders.hpp"
#include "DenseMatrixInverter.hpp"
#include "Types.hpp"
#include <memory>

namespace blitzdg {
    /**
     * Provides facilities for the construction of one-dimensional
     * nodes, operators, and geometric factors.
     * @note This class is moveable but not copyable.
     */ 
    class VandermondeBuilders {
        JacobiBuilders Jacobi;
        DenseMatrixInverter Inverter;
public:
        /**
         * Constructor.
         */
        VandermondeBuilders() 
            : Jacobi{}, Inverter{}
        {}

       /**
        * Computes the Vandermonde matrix which maps modal coefficients to nodal values.
        * Also computes the inverse and stores in an output parameter.
        * @param[in] r The nodal points at which to evaluate the modal basis functions.
        * @param[out] V The generalized Vandermonde matrix (dense).
        * @param[out] Vinv The inverse of the generalized Vandermonde matrix.
        */
        void computeVandermondeMatrix(const real_vector_type & r, real_matrix_type& V, real_matrix_type& Vinv) const {
            const index_type Np = r.length(0);
            real_vector_type p(Np);
            for (index_type j=0; j < Np; j++) {
                Jacobi.computeJacobiPolynomial(r, 0.0, 0.0, j, p);
                V(blitz::Range::all(), j) = p;
            }

            Inverter.computeInverse(V, Vinv);
        }

       /**
        * Computes the elementwise derivative of the Vandermonde matrix.
        * @param[in] r The nodal points at which to evaluate the modal basis functions.
        * @param[out] DVr Elementwise derivative of Vandermonde matrix.
        */
        void computeGradVandermonde(const real_vector_type & r, real_matrix_type& DVr) const {
            const index_type Np = r.length(0);
            blitz::firstIndex ii;
            real_vector_type dp(Np);
            for (index_type i=0; i< Np; i++) {
                dp = 0.*ii;
                Jacobi.computeGradJacobi(r, 0.0, 0.0, i, dp);
                DVr(blitz::Range::all(), i) = dp;
            }    
        }
    };
}