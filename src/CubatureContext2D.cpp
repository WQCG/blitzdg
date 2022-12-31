// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "Nodes1DProvisioner.hpp"
#include "CSCMatrix.hpp"
#include "BlitzHelpers.hpp"
#include "VandermondeBuilders.hpp"
#include "Types.hpp"
#include "CubatureContext2D.hpp"
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

    ndarray CubatureContext2D::x_numpy() const {
        Py_intptr_t shape[2] = { x_->rows(), x_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(x_->begin(), x_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::y_numpy() const {
        Py_intptr_t shape[2] = { y_->rows(), y_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(y_->begin(), y_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::W_numpy() const {
        Py_intptr_t shape[2] = { W_->rows(), W_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(W_->begin(), W_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::r_numpy() const {
        Py_intptr_t shape[1] = { r_->length(0) };
        ndarray result = zeros(1, shape, dtype::get_builtin<real_type>());
        std::copy(r_->begin(), r_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::s_numpy() const {
        Py_intptr_t shape[1] = { s_->length(0) };
        ndarray result = zeros(1, shape, dtype::get_builtin<real_type>());
        std::copy(s_->begin(), s_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::w_numpy() const {
        Py_intptr_t shape[1] = { w_->length(0) };
        ndarray result = zeros(1, shape, dtype::get_builtin<real_type>());
        std::copy(w_->begin(), w_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::Jac_numpy() const {
        Py_intptr_t shape[2] = { J_->rows(), J_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(J_->begin(), J_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::rx_numpy() const {
        Py_intptr_t shape[2] = { rx_->rows(), rx_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(rx_->begin(), rx_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::ry_numpy() const {
        Py_intptr_t shape[2] = { ry_->rows(), ry_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(ry_->begin(), ry_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::sx_numpy() const {
        Py_intptr_t shape[2] = { sx_->rows(), sx_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(sx_->begin(), sx_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::sy_numpy() const {
        Py_intptr_t shape[2] = { sy_->rows(), sy_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(sy_->begin(), sy_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::Dr_numpy() const {
        Py_intptr_t shape[2] = { Dr_->rows(), Dr_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Dr_->begin(), Dr_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::Ds_numpy() const {
        Py_intptr_t shape[2] = { Ds_->rows(), Ds_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Ds_->begin(), Ds_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::MM_numpy() const {
        Py_intptr_t shape[3] = { MM_->length(0), MM_->length(1),MM_->length(2) };
        ndarray result = zeros(3, shape, dtype::get_builtin<real_type>());
        std::copy(MM_->begin(), MM_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::MMChol_numpy() const {
        Py_intptr_t shape[3] = { MMChol_->length(0), MMChol_->length(1),MMChol_->length(2) };
        ndarray result = zeros(3, shape, dtype::get_builtin<real_type>());
        std::copy(MMChol_->begin(), MMChol_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray CubatureContext2D::V_numpy() const {
        Py_intptr_t shape[2] = { V_->rows(), V_->cols() };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(V_->begin(), V_->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }
}