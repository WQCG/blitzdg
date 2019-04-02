// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include "Types.hpp"

namespace blitzdg {
    class DGContext2D {
    private:
        index_type Np_;
        index_type Nfp_;
        index_type K_;
        index_type NumFaces_;
        real_matrix_type* Filt_;
        real_matrix_type* x_;
        real_matrix_type* y_;
        real_matrix_type* xGrid_;
        real_matrix_type* yGrid_;
        real_matrix_type* Fscale_;
        real_matrix_type* V_;
        real_matrix_type* Vinv_;
        real_matrix_type* J_;
        real_matrix_type* rx_;
        real_matrix_type* ry_;
        real_matrix_type* sx_;
        real_matrix_type* sy_;
        real_matrix_type* nx_;
        real_matrix_type* ny_;
        real_matrix_type* Dr_;
        real_matrix_type* Ds_;
        real_matrix_type* Lift_;
        index_vector_type* vmapM_;
        index_vector_type* vmapP_;
        index_hashmap* bcHash_;
    public:
        DGContext2D() = default;
        
        DGContext2D(
            index_type numLocalPoints,
            index_type numFacePoints,
            index_type numElems,
            index_type numFaces,
            real_matrix_type* filter,
            real_matrix_type* xgrid,
            real_matrix_type* ygrid,
            real_matrix_type* fscale,
            real_matrix_type* vandermonde2d,
            real_matrix_type* vandermonde2dinv,
            real_matrix_type* jacobian,
            real_matrix_type* rx,
            real_matrix_type* ry,
            real_matrix_type* sx,
            real_matrix_type* sy,
            real_matrix_type* nx,
            real_matrix_type* ny,
            real_matrix_type* Dr,
            real_matrix_type* Ds,
            real_matrix_type* lift,
            index_vector_type* vmapM,
            index_vector_type* vmapP,
            index_hashmap* bcmap)
            : Np_{ numLocalPoints }, Nfp_{ numFacePoints },
            K_{ numElems }, NumFaces_{ numFaces },
            Filt_{ filter }, xGrid_{ xgrid }, yGrid_{ ygrid },
            Fscale_{ fscale }, V_ { vandermonde2d }, Vinv_ {vandermonde2dinv}, J_{ jacobian }, rx_{ rx },
            ry_{ ry }, sx_{ sx }, sy_{ sy }, nx_{ nx }, ny_{ ny },
            Dr_{ Dr }, Ds_{ Ds }, Lift_{ lift }, vmapM_{ vmapM },
            vmapP_{ vmapP }, bcHash_{ bcmap }
        {}
        
        index_type numLocalPoints() const {
            return Np_;
        }
        
        index_type numFacePoints() const {
            return Nfp_;
        }
        
        index_type numElements() const {
            return K_;
        }
        
        index_type numFaces() const {
            return NumFaces_;
        }
        
        const real_matrix_type& filter() const {
            return *Filt_;
        }
        
        const real_matrix_type& x() const {
            return *xGrid_;
        }
        
        const real_matrix_type& y() const {
            return *yGrid_;
        }
        
        const real_matrix_type& fscale() const {
            return *Fscale_;
        }

        const real_matrix_type& V() const {
            return *V_;
        }

        const real_matrix_type& Vinv() const {
            return *Vinv_;
        }
        
        const real_matrix_type& jacobian() const {
            return *J_;
        }
        
        const real_matrix_type& rx() const {
            return *rx_;
        }
        
        const real_matrix_type& sx() const {
            return *sx_;
        }
        
        const real_matrix_type& ry() const {
            return *ry_;
        }
        
        const real_matrix_type& sy() const {
            return *sy_;
        }
        
        const real_matrix_type& nx() const {
            return *nx_;
        }
        
        const real_matrix_type& ny() const {
            return *ny_;
        }
        
        const real_matrix_type& Dr() const {
            return *Dr_;
        }
        
        const real_matrix_type& Ds() const {
            return *Ds_;
        }
        
        const real_matrix_type& lift() const {
            return *Lift_;
        }
        
        const index_vector_type& vmapM() const {
            return *vmapM_;
        }
        
        const index_vector_type& vmapP() const {
            return *vmapP_;
        }
        
        const index_hashmap& bcmap() const {
            return *bcHash_;
        }
    };
}