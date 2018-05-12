#include "Types.hpp"
#include <cmath>

namespace blitzdg {
// Computes the Givens rotation values c and s so that
// [ c s ] [a] = [r]
// [-s c ] [b]   [0]
// Overwrites a and b with r and 0, respectively.
void DROTG(real_type& a, real_type& b, real_type& c, real_type& s);

// Applies a Givens rotation
// x <-- c*x + s*y
// y <-- c*y - s*x
inline void applyGivens(real_type c, real_type s, real_type& x, real_type& y) {
    real_type tmp = c * x + s * y;
    y = c * y - s * x;
    x = tmp;
}

// Vector 1-norm
inline real_type norm1(const vector_type& x) {
   return blitz::sum(blitz::abs(x));
}

// Vector 2-norm
inline real_type norm2(const vector_type& x) {
    return std::sqrt(blitz::sum(x * x));
}

// vector infinity-norm
inline real_type normInf(const vector_type& x) {
    return blitz::max(blitz::abs(x));
}
} // namespace blitzdg