// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include "Types.hpp"
#include <boost/python/numpy.hpp>

namespace blitzdg {
    class DGContext2D {
    private:
        index_type N_;
        index_type Np_;
        index_type Nfp_;
        index_type K_;
        index_type NumFaces_;
        real_matrix_type* Filt_;
        real_vector_type* r_;
        real_vector_type* s_;
        real_matrix_type* xGrid_;
        real_matrix_type* yGrid_;
        real_matrix_type* Fscale_;
        index_matrix_type* Fmask_;
        std::vector<index_type>* Gather_;
        std::vector<index_type>* Scatter_;
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
            index_type order,
            index_type numLocalPoints,
            index_type numFacePoints,
            index_type numElems,
            index_type numFaces,
            real_matrix_type* filter,
            real_vector_type* rgrid,
            real_vector_type* sgrid,
            real_matrix_type* xgrid,
            real_matrix_type* ygrid,
            real_matrix_type* fscale,
            index_matrix_type* fmask,
            std::vector<index_type>* gather,
            std::vector<index_type>* scatter,
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
            : N_{ order }, Np_{ numLocalPoints }, Nfp_{ numFacePoints },
            K_{ numElems }, NumFaces_{ numFaces },
            Filt_{ filter }, r_{ rgrid }, s_{ sgrid },  xGrid_{ xgrid }, yGrid_{ ygrid },
            Fscale_{ fscale }, Fmask_ { fmask }, Gather_{ gather }, Scatter_{ scatter }, V_ { vandermonde2d }, Vinv_ {vandermonde2dinv}, J_{ jacobian }, rx_{ rx },
            ry_{ ry }, sx_{ sx }, sy_{ sy }, nx_{ nx }, ny_{ ny },
            Dr_{ Dr }, Ds_{ Ds }, Lift_{ lift }, vmapM_{ vmapM },
            vmapP_{ vmapP }, bcHash_{ bcmap }
        {}

        index_type order() const {
            return N_;
        }
        
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

        const real_vector_type& r() const {
            return *r_;
        }

        const real_vector_type& s() const {
            return *s_;
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

        const index_matrix_type& fmask() const {
            return *Fmask_;
        }

        const std::vector<index_type>& gather() const {
            return *Gather_;
        }

        const std::vector<index_type>& scatter() const {
            return *Scatter_;
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

        boost::python::numpy::ndarray filter_numpy() const;
        boost::python::numpy::ndarray r_numpy() const;
        boost::python::numpy::ndarray s_numpy() const;
        boost::python::numpy::ndarray x_numpy() const;
        boost::python::numpy::ndarray y_numpy() const;
        boost::python::numpy::ndarray fscale_numpy() const;
        boost::python::numpy::ndarray fmask_numpy() const;
        boost::python::numpy::ndarray gather_numpy() const;
        boost::python::numpy::ndarray scatter_numpy() const;
        boost::python::numpy::ndarray jacobian_numpy() const;
        boost::python::numpy::ndarray rx_numpy() const;
        boost::python::numpy::ndarray ry_numpy() const;
        boost::python::numpy::ndarray sx_numpy() const;
        boost::python::numpy::ndarray sy_numpy() const;
        boost::python::numpy::ndarray nx_numpy() const;
        boost::python::numpy::ndarray ny_numpy() const;
        boost::python::numpy::ndarray Dr_numpy() const;
        boost::python::numpy::ndarray Ds_numpy() const;
        boost::python::numpy::ndarray lift_numpy() const;
        boost::python::numpy::ndarray vmapM_numpy() const;
        boost::python::numpy::ndarray vmapP_numpy() const;
        boost::python::dict bcmap_numpy() const;
        boost::python::numpy::ndarray V_numpy() const;

        void computeDifferentiationMatrices(const real_vector_type& x, const real_vector_type& y, const real_matrix_type& V, real_matrix_type& Dx, real_matrix_type& Dy) {
            const real_matrix_type& Dr = *Dr_;
            const real_matrix_type& Ds = *Ds_;

            index_type Nout = V.rows();

            blitz::firstIndex ii;
            blitz::secondIndex jj;
            blitz::thirdIndex kk;

            real_vector_type xr(Np_), xs(Np_), yr(Np_), ys(Np_), J(Np_),
                rx(Np_), sx(Np_), ry(Np_), sy(Np_);

            xr = blitz::sum(Dr(ii, kk)*x(kk, jj), kk);
            xs = blitz::sum(Ds(ii, kk)*x(kk, jj), kk);

            yr = blitz::sum(Dr(ii, kk)*y(kk, jj), kk);
            ys = blitz::sum(Ds(ii, kk)*y(kk, jj), kk);

            J = -xs*yr + xr*ys;

            rx = ys/J;
            sx =-yr/J;
            ry =-xs/J;
            sy = xr/J;

            Dx = 0.0; Dy = 0.0;


            for (index_type i=0; i < Nout; ++i) {
                for (index_type j=0; j < Np_; ++j) {
                    Dx(i, j) = rx(i) * Dr(i, j) + sx(i) * Ds(i, j);
                    Dy(i, j) = ry(i) * Dr(i, j) + sy(i) * Ds(i, j);
                }
            }
        }
    };
}