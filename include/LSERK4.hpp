// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file LSERK4.hpp
 * @brief Header file containing coefficients for Low-Storage Explicit Runge-Kutta time-stepper (4th order).
 */
#pragma once
#include "Types.hpp"

namespace blitzdg {
    /**
     * Coefficients for Low-Storage Explicit Runge-Kutta 4th Order Time-Stepper.
     */
    namespace LSERK4 {
        const index_type numStages = 5;
        const real_type rk4a[numStages] = { 0.0,
                                            -567301805773.0/1357537059087.0,
                                            -2404267990393.0/2016746695238.0,
                                            -3550918686646.0/2091501179385.0,
                                            -1275806237668.0/842570457699.0 };

        const real_type rk4b[numStages] = { 1432997174477.0/9575080441755.0,
                                            5161836677717.0/13612068292357.0,
                                            1720146321549.0/2090206949498.0,
                                            3134564353537.0/4481467310338.0,
                                            2277821191437.0/14882151754819.0 };

    } // namespace LSERK4
} // namespace blitzdg