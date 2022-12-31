// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "LinAlgHelpers.hpp"

namespace blitzdg {
    extern "C" {
        void drotg_(double* a, double* b, double* c, double* s);
    }

    void drotg(real_type& a, real_type& b, real_type& c, real_type& s) {
        drotg_(&a, &b, &c, &s);
        b = real_type(0);
    }
} // namespace blitzdg
