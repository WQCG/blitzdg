// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file BCtypes.hpp
 * @brief Defines enum of boundary condition tags.
 */
namespace blitzdg {
    enum BCTag {
        In = 1,
        Out = 2,
        Wall = 3,
        Far = 4,
        Cyl = 5,
        Dirichlet = 6,
        Neuman = 7,
        Slip = 8
    };
}