// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include "Types.hpp"
#include "BlitzHelpers.hpp"
#include <vector>
#include <memory>
#include <boost/python/numpy.hpp>

namespace blitzdg {
	class CubatureContext2D {
	private:
		index_type NCubature_;
        index_type NumCubaturePoints_;
        real_vec_smart_ptr r_;
        real_vec_smart_ptr s_;
        real_vec_smart_ptr w_;
        real_mat_smart_ptr V_;
        real_mat_smart_ptr rx_;
        real_mat_smart_ptr sx_;
        real_mat_smart_ptr ry_;
        real_mat_smart_ptr sy_;
        real_mat_smart_ptr J_;
        real_mat_smart_ptr Dr_;
        real_mat_smart_ptr Ds_;
        std::unique_ptr<real_tensor3_type> MM_;
        std::unique_ptr<real_tensor3_type> MMChol_;
        real_mat_smart_ptr x_;
        real_mat_smart_ptr y_;
        real_mat_smart_ptr W_;



    public:
        CubatureContext2D() = default;
        CubatureContext2D(
            index_type NCubature,
            index_type NumCubaturePoints,
            const real_vector_type& r,
            const real_vector_type& s,
            const real_vector_type& w,
            const real_matrix_type& V,
            const real_matrix_type& rx,
            const real_matrix_type& sx,
            const real_matrix_type& ry,
            const real_matrix_type& sy,
            const real_matrix_type& J,
            const real_matrix_type& Dr,
            const real_matrix_type& Ds,
            const real_tensor3_type& MM,
            const real_tensor3_type& MMChol,
            const real_matrix_type& x,
            const real_matrix_type& y,
            const real_matrix_type& W    
        ) : 
            NCubature_{ NCubature },
            NumCubaturePoints_{ NumCubaturePoints },
            r_{ std::make_unique<real_vector_type>(r) },
            s_{ std::make_unique<real_vector_type>(s) },
            w_{ std::make_unique<real_vector_type>(w) },
            V_{ std::make_unique<real_matrix_type>(V) },
            rx_{ std::make_unique<real_matrix_type>(rx) },
            sx_{ std::make_unique<real_matrix_type>(sx) },
            ry_{ std::make_unique<real_matrix_type>(ry) },
            sy_{ std::make_unique<real_matrix_type>(sy) },
            J_{ std::make_unique<real_matrix_type>(J) },
            Dr_{ std::make_unique<real_matrix_type>(Dr) },
            Ds_{ std::make_unique<real_matrix_type>(Ds) },
            MM_{ std::make_unique<real_tensor3_type>(MM) },
            MMChol_{ std::make_unique<real_tensor3_type>(MMChol) },
            x_{ std::make_unique<real_matrix_type>(x) },
            y_{ std::make_unique<real_matrix_type>(y) },
            W_{ std::make_unique<real_matrix_type>(W) }
        {}

        const index_type NCubature() const { return NCubature_; }
        const index_type NumCubaturePoints() const { return NumCubaturePoints_; }

        const real_vector_type& r() const { return *r_; }
        const real_vector_type& s() const { return *s_; }
        const real_vector_type& w() const { return *w_; }
        
        const real_matrix_type& V() const { return *V_; }
        const real_matrix_type& rx() const { return *rx_; }
        const real_matrix_type& sx() const { return *sx_; }
        const real_matrix_type& ry() const { return *ry_; }
        const real_matrix_type& sy() const { return *sy_; }
        const real_matrix_type& Jac() const { return *J_; }
        const real_matrix_type& Dr() const { return *Dr_; }
        const real_matrix_type& Ds() const { return *Ds_; }

        const real_tensor3_type& MM() const { return *MM_; }
        const real_tensor3_type& MMChol() const { return *MMChol_; }

        const real_matrix_type& x() const { return *x_; }
        const real_matrix_type& y() const { return *y_; }
        const real_matrix_type& W() const { return *W_; }

        using numpyarray = boost::python::numpy::ndarray;

        // for py bindings
        numpyarray r_numpy() const;
        numpyarray s_numpy() const;
        numpyarray w_numpy() const;
        numpyarray V_numpy() const;
        numpyarray rx_numpy() const;
        numpyarray ry_numpy() const;
        numpyarray sx_numpy() const;
        numpyarray sy_numpy() const;
        numpyarray Jac_numpy() const;
        numpyarray Dr_numpy() const;
        numpyarray Ds_numpy() const;

        numpyarray MM_numpy() const;
        numpyarray MMChol_numpy() const;
        
        numpyarray x_numpy() const;
        numpyarray y_numpy() const;
        numpyarray W_numpy() const;

    };
}