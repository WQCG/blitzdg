// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "Nodes1DProvisioner.hpp"
#include "CSCMatrix.hpp"
#include "BlitzHelpers.hpp"
#include "VandermondeBuilders.hpp"
#include "Types.hpp"
#include "GaussFaceContext2D.hpp"
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

    ndarray GaussFaceContext2D::x_numpy() const {
        Py_intptr_t shape[2] = { x_->rows(), x_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(x_->begin(), x_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::y_numpy() const {
        Py_intptr_t shape[2] = { y_->rows(), y_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(y_->begin(), y_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::sJ_numpy() const {
        Py_intptr_t shape[2] = { sJ_->rows(), sJ_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(sJ_->begin(), sJ_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::Jac_numpy() const {
        Py_intptr_t shape[2] = { Jac_->rows(), Jac_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Jac_->begin(), Jac_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::rx_numpy() const {
        Py_intptr_t shape[2] = { rx_->rows(), rx_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(rx_->begin(), rx_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::ry_numpy() const {
        Py_intptr_t shape[2] = { ry_->rows(), ry_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(ry_->begin(), ry_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::sx_numpy() const {
        Py_intptr_t shape[2] = { sx_->rows(), sx_->rows() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(sx_->begin(), sx_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::sy_numpy() const {
        Py_intptr_t shape[2] = { sy_->rows(), sy_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(sy_->begin(), sy_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::nx_numpy() const {
        Py_intptr_t shape[2] = { nx_->rows(), nx_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(nx_->begin(), nx_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::ny_numpy() const {
        Py_intptr_t shape[2] = { ny_->rows(), ny_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(ny_->begin(), ny_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    dict GaussFaceContext2D::bcMap_numpy() const {
        boost::python::dict bcDict;

        for (const auto& t : *bcMap_) {
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

    ndarray GaussFaceContext2D::W_numpy() const {
        Py_intptr_t shape[2] = { W_->rows(), W_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(W_->begin(), W_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::Interp_numpy() const {
        Py_intptr_t shape[2] = { Interp_->rows(), Interp_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Interp_->begin(), Interp_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    const real_matrix_type& GaussFaceContext2D::Interp() const {
        return *Interp_;
    }

    ndarray GaussFaceContext2D::mapM_numpy() const {
        Py_intptr_t shape[1] = { mapM_->length(0) };
        ndarray result = zeros(1, shape, dtype::get_builtin<index_type>());
        std::copy(mapM_->begin(), mapM_->end(), reinterpret_cast<index_type*>(result.get_data()));
        return result;
    }

    ndarray GaussFaceContext2D::mapP_numpy() const {
        Py_intptr_t shape[1] = { mapP_->length(0) };
        ndarray result = zeros(1, shape, dtype::get_builtin<index_type>());
        std::copy(mapP_->begin(), mapP_->end(), reinterpret_cast<index_type*>(result.get_data()));
        return result;
    }
}