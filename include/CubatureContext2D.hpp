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
        std::shared_ptr<real_vector_type> r_;
        std::shared_ptr<real_vector_type> s_;
        std::shared_ptr<real_vector_type> w_;
        std::shared_ptr<real_matrix_type> V_;
        std::shared_ptr<real_matrix_type> rx_;
        std::shared_ptr<real_matrix_type> sx_;
        std::shared_ptr<real_matrix_type> ry_;
        std::shared_ptr<real_matrix_type> sy_;
        std::shared_ptr<real_matrix_type> J_;
        std::shared_ptr<real_matrix_type> Dr_;
        std::shared_ptr<real_matrix_type> Ds_;
        std::shared_ptr<real_tensor3_type> MM_;
        std::shared_ptr<real_tensor3_type> MMChol_;
        std::shared_ptr<real_matrix_type> x_;
        std::shared_ptr<real_matrix_type> y_;
        std::shared_ptr<real_matrix_type> W_;

    public:
        CubatureContext2D() = default;
        CubatureContext2D(
            index_type NCubature,
            index_type NumCubaturePoints,
            std::shared_ptr<real_vector_type> r,
            std::shared_ptr<real_vector_type> s,
            std::shared_ptr<real_vector_type> w,
            std::shared_ptr<real_matrix_type> V,
            std::shared_ptr<real_matrix_type> rx,
            std::shared_ptr<real_matrix_type> sx,
            std::shared_ptr<real_matrix_type> ry,
            std::shared_ptr<real_matrix_type> sy,
            std::shared_ptr<real_matrix_type> J,
            std::shared_ptr<real_matrix_type> Dr,
            std::shared_ptr<real_matrix_type> Ds,
            std::shared_ptr<real_tensor3_type> MM,
            std::shared_ptr<real_tensor3_type> MMChol,
            std::shared_ptr<real_matrix_type> x,
            std::shared_ptr<real_matrix_type> y,
            std::shared_ptr<real_matrix_type> W    
        ) : 
            NCubature_{ NCubature },
            NumCubaturePoints_{ NumCubaturePoints },
            r_{ r },
            s_{ s },
            w_{ w },
            V_{ V },
            rx_{ rx },
            sx_{ sx },
            ry_{ ry },
            sy_{ sy },
            J_{ J },
            Dr_{ Dr },
            Ds_{ Ds },
            MM_{ MM },
            MMChol_{ MMChol },
            x_{ x },
            y_{ y },
            W_{ W }
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