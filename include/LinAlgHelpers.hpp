#include "Types.hpp"
#include <cmath>

namespace blitzdg {
    /**
     * Given Cartesian coordinates (a,b), computes the
     * parameters c, s of the associated Givens rotation.
     * 
     * The Givens rotation G(c,s) is defined so that
     *     c*a + s*b = r
     *     c*b - s*a = 0
     * @param[in,out] a real scalar
     * @param[in,out] b real scalar
     * @param[in,out] c real scalar
     * @param[in,out] s real scalar
     * On output:
     *     a = r
     *     b = 0
     */
    void DROTG(real_type& a, real_type& b, real_type& c, real_type& s);

    /**
     * Applies a Givens rotation G(c,s) to the tuple (x,y).
     * @param[in] c real scalar
     * @param[in] s real scalar
     * @param[in,out] x real scalar
     * @param[in,out] y real scalar
     * On output:
     *     x = c*x + s*y
     *     y = c*y - s*x
     */ 
    inline void applyGivens(real_type c, real_type s, real_type& x, real_type& y) {
        real_type tmp = c * x + s * y;
        y = c * y - s * x;
        x = tmp;
    }

    /**
     * Vector 1-norm.
     * @param[in] x real-valued vector
     */
    inline real_type norm1(const vector_type& x) {
        return blitz::sum(blitz::abs(x));
    }

    /**
     * Vector 2-norm.
     * @param[in] x real-valued vector
     */
    inline real_type norm2(const vector_type& x) {
        return std::sqrt(blitz::sum(x * x));
    }

    /**
     * Vector infinity-norm.
     * @param[in] x real-valued vector
     */
    inline real_type normInf(const vector_type& x) {
        return blitz::max(blitz::abs(x));
    }
} // namespace blitzdg