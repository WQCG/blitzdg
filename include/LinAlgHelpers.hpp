/**
 * @file LinAlgHelpers.hpp
 * @brief Implements a set functions for basic linear algebra operations.
 */
#include "Types.hpp"
#include <cmath>

namespace blitzdg {
    /**
     * Given Cartesian coordinates \f$(a,b)\f$, computes the
     * parameters \f$c,\,s\f$ of the associated Givens rotation.
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
    void DROTG(real_type& a, real_type& b, real_type& c, real_type& s);

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
     * Returns the one-norm of the input vector x.
     */
    inline real_type norm1(const vector_type& x) {
        return blitz::sum(blitz::abs(x));
    }

    /**
     * Returns the two-norm of the input vector x.
     */
    inline real_type norm2(const vector_type& x) {
        return std::sqrt(blitz::sum(x * x));
    }

    /**
     * Returns the infinity-norm of the input vector x.
     */
    inline real_type normInf(const vector_type& x) {
        return blitz::max(blitz::abs(x));
    }
} // namespace blitzdg