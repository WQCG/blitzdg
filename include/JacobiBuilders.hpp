// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file JacobiBuilders.hpp
 * @brief Defines the JacobiBuilders class that constructs Jacobi polynomials, their
 * derivatives, and associated quadrature points. This class is also responsible for
 * constructing the Jacobi(Legendre)-Gauss-Lobotto nodes.
 */
#pragma once
#include "EigenSolver.hpp"
#include "Types.hpp"

namespace blitzdg {
    class JacobiBuilders {
        EigenSolver EigSolver;

    public:

        /**
         * Default constructor.
         */
        JacobiBuilders()
            : EigSolver()
        {}

        /**
         * Computes the value of the Nth order Jacobi polynomial with parameters 
         * \f$\alpha\f$ and \f$\beta\f$ at each point in the input vector x.
         * @param[in] x Points at which to evaluate the Jacobi polynomial.
         * @param[in] alpha Jacobi polynomial parameter.
         * @param[in] beta Jacobi polynomial parameter.
         * @param[in] N Order of the Jacobi polynomial.
         * @param[out] p Jacobi polynomial values at the points in x.
         * @note Jacobi parameters must satisfy \f$\alpha,\beta > -1\f$ and
         * \f$\alpha,\beta \neq -1/2\f$.
         */
        void computeJacobiPolynomial(const real_vector_type& x, real_type alpha, real_type beta, index_type N, real_vector_type& p) const;
        
        /**
         * Computes the nodes and weights of the Nth order Gauss-Jacobi-Lobotto quadrature rule
         * with parameters \f$\alpha\f$ and \f$\beta\f$.
         * @param[in] alpha Jacobi polynomial parameter.
         * @param[in] beta Jacobi polynomial parameter.
         * @param[in] N Order of the quadrature rule.
         * @param[out] x Length N array of quadrature rule nodes.
         * @param[out] w Length N array of quadrature weights. 
         * @note Jacobi parameters must satisfy \f$\alpha,\beta > -1\f$ and
         * \f$\alpha,\beta \neq -1/2\f$.
         */
        void computeJacobiQuadWeights(real_type alpha, real_type beta, index_type N, real_vector_type & x, real_vector_type & w) const;
        
        /**
         * Computes the nodes of the Nth order Gauss-Jacobi-Lobotto quadrature rule
         * with parameters \f$\alpha\f$ and \f$\beta\f$.
         * @param[in] alpha Jacobi polynomial parameter.
         * @param[in] beta Jacobi polynomial parameter.
         * @param[in] N Order of the quadrature rule.
         * @param[out] x Length N array of quadrature rule nodes.
         * @note Jacobi parameters must satisfy \f$\alpha,\beta > -1\f$ and
         * \f$\alpha,\beta \neq -1/2\f$.
         */
        void computeGaussLobottoPoints(real_type alpha, real_type beta, index_type N, real_vector_type & x) const;
        
        /**
         * Computes the first derivative of the Nth order Jacobi polynomial with parameters
         * \f$\alpha\f$ and \f$\beta\f$ at each point in the input vector x.
         * @param[in] x Points at which to evaluate the Jacobi polynomial derivative.
         * @param[in] alpha Jacobi polynomial parameter.
         * @param[in] beta Jacobi polynomial parameter.
         * @param[in] N Order of the Jacobi polynomial.
         * @param[out] dp Jacobi polynomial first derivative values at the points in x.
         * @note Jacobi parameters must satisfy \f$\alpha,\beta > -1\f$ and
         * \f$\alpha,\beta \neq -1/2\f$.
         */
        void computeGradJacobi(const real_vector_type& x, real_type alpha, real_type beta, index_type N, real_vector_type& dp) const;
    }; // class JacobiBuilders
} // namespace blitzdg