#pragma once

#include <suitesparse/umfpack.h>
#include "CSCMatrix.hpp"
#include "SparseTriplet.hpp"
#include "DGContext2D.hpp"
#include "MeshManager.hpp"

namespace blitzdg {
    class Poisson2DSparseMatrix {
    public:
        Poisson2DSparseMatrix(DGContext2D& dg, MeshManager& mshManager);
    };
}