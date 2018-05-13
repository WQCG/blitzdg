#include "LinAlgHelpers.hpp"

namespace blitzdg {
    extern "C" {
        void drotg_(double* a, double* b, double* c, double* s);
    }

    void DROTG(real_type& a, real_type& b, real_type& c, real_type& s) {
        drotg_(&a, &b, &c, &s);
        b = real_type(0);
    }
} // namespace blitzdg
