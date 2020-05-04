// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "Nodes1DProvisioner.hpp"
#include "CSCMatrix.hpp"
#include "BlitzHelpers.hpp"
#include "VandermondeBuilders.hpp"
#include "Types.hpp"
#include "DGContext2D.hpp"
#include <blitz/array.h>
#include <cmath>
#include <limits>
#include <memory>
#include <boost/python/numpy.hpp>

using boost::python::numpy::ndarray;
using boost::python::numpy::zeros;
using boost::python::numpy::dtype;
using boost::python::dict;
using boost::python::list;


namespace blitzdg {
    ndarray DGContext2D::filter_numpy() const {
        Py_intptr_t shape[2] = { Np_, Np_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Filt_->begin(), Filt_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::r_numpy() const {
        Py_intptr_t shape[1] = { Np_ };
        ndarray result = zeros(1, shape, dtype::get_builtin<real_type>());
        std::copy(r_->begin(), r_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::s_numpy() const {
        Py_intptr_t shape[1] = { Np_ };
        ndarray result = zeros(1, shape, dtype::get_builtin<real_type>());
        std::copy(s_->begin(), s_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::x_numpy() const {
        Py_intptr_t shape[2] = { Np_, K_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(xGrid_->begin(), xGrid_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::y_numpy() const {
        Py_intptr_t shape[2] = { Np_, K_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(yGrid_->begin(), yGrid_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::fscale_numpy() const {
        Py_intptr_t shape[2] = { NumFaces_*Nfp_, K_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Fscale_->begin(), Fscale_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::jacobian_numpy() const {
        Py_intptr_t shape[2] = { Np_, K_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(J_->begin(), J_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::rx_numpy() const {
        Py_intptr_t shape[2] = { Np_, K_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(rx_->begin(), rx_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::ry_numpy() const {
        Py_intptr_t shape[2] = { Np_, K_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(ry_->begin(), ry_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::sx_numpy() const {
        Py_intptr_t shape[2] = { Np_, K_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(sx_->begin(), sx_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::sy_numpy() const {
        Py_intptr_t shape[2] = { Np_, K_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(sy_->begin(), sy_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::nx_numpy() const {
        Py_intptr_t shape[2] = { NumFaces_*Nfp_, K_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(nx_->begin(), nx_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::ny_numpy() const {
        Py_intptr_t shape[2] = { NumFaces_*Nfp_, K_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(ny_->begin(), ny_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::Dr_numpy() const {
        Py_intptr_t shape[2] = { Np_, Np_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Dr_->begin(), Dr_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::Ds_numpy() const {
        Py_intptr_t shape[2] = { Np_, Np_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Ds_->begin(), Ds_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::lift_numpy() const {
        Py_intptr_t shape[2] = { Np_, NumFaces_*Nfp_ };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Lift_->begin(), Lift_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::vmapM_numpy() const {
        Py_intptr_t shape[1] = {  NumFaces_*Nfp_*K_ };
        ndarray result = zeros(1, shape, dtype::get_builtin<index_type>());
        std::copy(vmapM_->begin(), vmapM_->end(), reinterpret_cast<index_type*>(result.get_data()));
        return result;
    }

    ndarray DGContext2D::vmapP_numpy() const {
        Py_intptr_t shape[1] = {  NumFaces_*Nfp_*K_ };
        ndarray result = zeros(1, shape, dtype::get_builtin<index_type>());
        std::copy(vmapP_->begin(), vmapP_->end(), reinterpret_cast<index_type*>(result.get_data()));
        return result;
    }

    dict DGContext2D::bcmap_numpy() const {
        boost::python::dict bcDict;

        for (const auto& t : *bcHash_) {
            auto key = t.first;
            const auto& vec = t.second;

            boost::python::list bcList;
            for (auto n : vec) {
                bcList.append( n );
            }
            bcDict[key] = bcList;
        }

        return bcDict;
    }
}