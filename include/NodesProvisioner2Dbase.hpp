#pragma once
#include "Types.hpp"

namespace blitzdg {
    class NodesProvisioner2DBase {
public:
    virtual int get_NOrder() const = 0;
    virtual int get_NumLocalPoints() const = 0;
    virtual void splitElements(const real_matrix_type& x, const real_matrix_type& y, const real_matrix_type& field,
        real_matrix_type& xnew, real_matrix_type& ynew, real_matrix_type& fieldnew) const = 0;
    virtual const real_matrix_type& get_xGrid() const = 0;
    virtual const real_matrix_type& get_yGrid() const = 0;
    virtual ~NodesProvisioner2DBase() = default;
    };
}