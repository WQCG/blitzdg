// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include "Types.hpp"
#include <boost/python/numpy.hpp>

namespace blitzdg {
    class GaussFaceContext2D {
    private:
        index_type NGauss_;
        real_mat_smart_ptr nx_;
        real_mat_smart_ptr ny_;
        real_mat_smart_ptr sJ_;
        real_mat_smart_ptr Jac_;
        real_mat_smart_ptr rx_;
        real_mat_smart_ptr ry_;
        real_mat_smart_ptr sx_;
        real_mat_smart_ptr sy_;
        std::unique_ptr<index_hashmap> bcMap_;
        real_mat_smart_ptr x_;
        real_mat_smart_ptr y_;
        real_mat_smart_ptr W_;

    public:
        GaussFaceContext2D() = default;

        GaussFaceContext2D(
            index_type NGauss,
            const real_matrix_type& nx,
            const real_matrix_type& ny,
            const real_matrix_type& sJ,
            const real_matrix_type& Jac,
            const real_matrix_type& rx,
            const real_matrix_type& ry,
            const real_matrix_type& sx,
            const real_matrix_type& sy,
            const index_hashmap& bcMap,
            const real_matrix_type& x,
            const real_matrix_type& y,
            const real_matrix_type& W
        ) : 
            NGauss_{ NGauss },
            nx_{ std::make_unique<real_matrix_type>(nx) },
            ny_{ std::make_unique<real_matrix_type>(ny) },
            sJ_{ std::make_unique<real_matrix_type>(sJ) },
            Jac_{ std::make_unique<real_matrix_type>(Jac) },
            rx_{ std::make_unique<real_matrix_type>(rx) },
            ry_{ std::make_unique<real_matrix_type>(ry) },
            sx_{ std::make_unique<real_matrix_type>(sx) },
            sy_{ std::make_unique<real_matrix_type>(sy) },
            bcMap_{ std::make_unique<index_hashmap>(bcMap) },
            x_{ std::make_unique<real_matrix_type>(x) },
            y_{ std::make_unique<real_matrix_type>(y) },
            W_{ std::make_unique<real_matrix_type>(W) }

        {}
        
        const index_type NGauss() const { return NGauss_; }
        const real_matrix_type& nx() const { return *nx_; }
        const real_matrix_type& ny() const { return *ny_; }
        const real_matrix_type& sJ() const { return *sJ_; }
        const real_matrix_type& Jac() const { return *Jac_; }
        const real_matrix_type& rx() const { return *rx_; }
        const real_matrix_type& ry() const { return *ry_; }
        const real_matrix_type& sx() const { return *sx_; }
        const real_matrix_type& sy() const { return *sy_; }
        const index_hashmap& bcMap() const { return *bcMap_; }
        const real_matrix_type& W() const { return *W_; }
    };
}