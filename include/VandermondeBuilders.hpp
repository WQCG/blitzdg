// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
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
#include "boost/python.hpp"
#include "boost/python/numpy.hpp"
#include <memory>

namespace blitzdg {
    using pyarray = boost::python::numpy::ndarray;
    using boost::python::numpy::dtype;
    using boost::python::numpy::zeros;
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
        void computeVandermondeMatrix(const real_vector_type & r, real_matrix_type& V, real_matrix_type& Vinv, bool includeInverse=true) const {
            index_type Np = V.cols();
            index_type Nr = V.rows();
            real_vector_type p(Nr);
            for (index_type j=0; j < Np; j++) {
                Jacobi.computeJacobiPolynomial(r, 0.0, 0.0, j, p);
                V(blitz::Range::all(), j) = p;
            }
            if (!includeInverse) {
                return;
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

        boost::python::tuple buildVandermondeMatrix_numpy(pyarray r, bool includeInverse=true, index_type order=-1) {
            index_type M = r.shape(0);
            index_type N = M;
            if (order > -1) {
                N = order + 1;
            }
            real_matrix_type V(M, N), Vinv(N,N);
            real_vector_type rblitz(M);

            rblitz = 0.0; V = 0.0; Vinv = 0.0;

            char * raw = r.get_data();
            std::copy(&raw[0], &raw[M*sizeof(real_type)], reinterpret_cast<char*>(rblitz.data()));

            computeVandermondeMatrix(rblitz, V, Vinv, includeInverse);

            Py_intptr_t shape[2] = { M, N };
            pyarray Vpy = zeros(2, shape, dtype::get_builtin<real_type>());
            pyarray Vinvpy = zeros(2, shape, dtype::get_builtin<real_type>());

            raw = reinterpret_cast<char*>(V.data());
            std::copy(&raw[0], &raw[M*N*sizeof(real_type)], Vpy.get_data());

            if (includeInverse) {
                raw = reinterpret_cast<char*>(Vinv.data());
                std::copy(&raw[0], &raw[M*N*sizeof(real_type)], Vinvpy.get_data());
                return boost::python::make_tuple<pyarray, pyarray>(Vpy, Vinvpy);
            }
            return boost::python::make_tuple<pyarray>(Vpy);
        }
    };
}