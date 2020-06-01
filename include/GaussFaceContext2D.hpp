// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include "Types.hpp"
#include <boost/python/numpy.hpp>

namespace blitzdg {
    
    class GaussFaceContext2D {
    private:
        index_type NGauss_;
        std::shared_ptr<real_matrix_type> nx_;
        std::shared_ptr<real_matrix_type> ny_;
        std::shared_ptr<real_matrix_type> sJ_;
        std::shared_ptr<real_matrix_type> Jac_;
        std::shared_ptr<real_matrix_type> rx_;
        std::shared_ptr<real_matrix_type> ry_;
        std::shared_ptr<real_matrix_type> sx_;
        std::shared_ptr<real_matrix_type> sy_;
        std::shared_ptr<index_hashmap> bcMap_;
        std::shared_ptr<real_matrix_type> x_;
        std::shared_ptr<real_matrix_type> y_;
        std::shared_ptr<real_matrix_type> W_;
        std::shared_ptr<real_matrix_type> Interp_;
        std::shared_ptr<index_vector_type> mapM_;
        std::shared_ptr<index_vector_type> mapP_;

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
            const real_matrix_type& W,
            const real_matrix_type& Interp,
            const index_vector_type& mapM,
            const index_vector_type& mapP
        ) : 
            NGauss_{ NGauss },
            nx_{ std::make_shared<real_matrix_type>(nx) },
            ny_{ std::make_shared<real_matrix_type>(ny) },
            sJ_{ std::make_shared<real_matrix_type>(sJ) },
            Jac_{ std::make_shared<real_matrix_type>(Jac) },
            rx_{ std::make_shared<real_matrix_type>(rx) },
            ry_{ std::make_shared<real_matrix_type>(ry) },
            sx_{ std::make_shared<real_matrix_type>(sx) },
            sy_{ std::make_shared<real_matrix_type>(sy) },
            bcMap_{ std::make_shared<index_hashmap>(bcMap) },
            x_{ std::make_shared<real_matrix_type>(x) },
            y_{ std::make_shared<real_matrix_type>(y) },
            W_{ std::make_shared<real_matrix_type>(W) },
            Interp_{ std::make_shared<real_matrix_type>(Interp) },
            mapM_{ std::make_shared<index_vector_type>(mapM) },
            mapP_{ std::make_shared<index_vector_type>(mapP) }
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
        const real_matrix_type& x() const { return *x_; }
        const real_matrix_type& y() const { return *y_; }
        const real_matrix_type& W() const { return *W_; }
        const real_matrix_type& Interp() const { return *Interp_; }

        const index_vector_type& mapM() const { return *mapM_; };
        const index_vector_type& mapP() const { return *mapP_; }; 

        using numpyarray = boost::python::numpy::ndarray;

        // for py bindings
        numpyarray sJ_numpy() const;
        numpyarray Jac_numpy() const;
        numpyarray rx_numpy() const;
        numpyarray ry_numpy() const;
        numpyarray sx_numpy() const;
        numpyarray sy_numpy() const;
        numpyarray nx_numpy() const;
        numpyarray ny_numpy() const;
        boost::python::dict bcMap_numpy() const;
        numpyarray x_numpy() const;
        numpyarray y_numpy() const;
        numpyarray W_numpy() const;
        numpyarray Interp_numpy() const;

        numpyarray mapM_numpy() const;
        numpyarray mapP_numpy() const;
    };
}