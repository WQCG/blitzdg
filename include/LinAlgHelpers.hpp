// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file LinAlgHelpers.hpp
 * @brief Implements a set functions for basic linear algebra operations.
 */
#pragma once
#include "Types.hpp"
#include <cmath>

namespace blitzdg {
    /**
     * Given Cartesian coordinates \f$(a,b)\f$, computes the
     * parameters \f$c,\,s\f$ of the associated Givens rotation.
     * Calls the BLAS routine DROTG.
     * 
     * The Givens rotation \f$G(c,s)\f$ is defined so that
     * \f{eqnarray*}{
     *     ca + sb &=& r\\
     *     cb - sa &=& 0
     * \f}
     * @param[in,out] a On input the \f$x\f$-coordinate.
     * @param[in,out] b On input the \f$y\f$-coordinate.
     * @param[out] c A real scalar.
     * @param[out] s A real scalar.
     * 
     * On output: \f$c,\,s\f$ are the Givens rotation parameters,
     * \f$a = r\f$, and \f$b = 0\f$.
     */
    void drotg(real_type& a, real_type& b, real_type& c, real_type& s);

    /**
     * Applies a Givens rotation \f$G(c,s)\f$ to the tuple \f$(x,y)\f$.
     * @param[in] c Givens rotation parameter.
     * @param[in] s Givens rotation parameter.
     * @param[in,out] x On input the \f$x\f$-coordinate.
     * @param[in,out] y On input the \f$y\f$-coordinate.
     * 
     * On output:
     * \f{eqnarray*}{
     *     x &=& cx + sy\\
     *     y &=& cy - sx
     * \f}
     */ 
    inline void applyGivens(real_type c, real_type s, real_type& x, real_type& y) {
        real_type tmp = c * x + s * y;
        y = c * y - s * x;
        x = tmp;
    }

    /**
     * Returns the one-norm of a vector x.
     */
    template <typename T>
    inline T norm1(const vector_type<T>& x) {
        return blitz::sum(blitz::abs(x));
    }

    /**
     * Returns the two-norm of a vector x.
     */
    template <typename T>
    inline real_type norm2(const vector_type<T>& x) {
        return std::sqrt(blitz::sum(x * x));
    }

    /**
     * Returns the Frobenius norm (two-norm) of a vector x. 
     */
    template <typename T>
    inline real_type normFro(const vector_type<T>& x) {
        return std::sqrt(blitz::sum(x * x));
    }

    /**
     * Returns the infinity-norm of a vector x.
     */
    template <typename T>
    inline T normInf(const vector_type<T>& x) {
        return blitz::max(blitz::abs(x));
    }

	/**
	 * Returns the maximum absolute value of all elements of a matrix mat.
	 */
    template <typename T>
	inline T normMax(const matrix_type<T>& mat) {
		return blitz::max(blitz::abs(mat));
	}

    /**
     * Returns the Frobenius norm of a matrix mat.
     */
    template <typename T>
    inline real_type normFro(const matrix_type<T>& mat) {
        return std::sqrt(blitz::sum(mat * mat));
    }

    /**
     * Helper method for checking that two points (x1,y1) and (x2,y2) are less than a distance eps apart.
     */
    inline bool distanceLessThanEps(real_type x1, real_type y1, real_type x2, real_type y2, real_type eps) {
        return std::hypot(x2 - x1, y2 - y1) < eps;
    }
} // namespace blitzdg
