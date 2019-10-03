#pragma once

#include <memory>
#include <suitesparse/umfpack.h>
#include <boost/python/numpy.hpp>
#include "CSCMatrix.hpp"
#include "SparseTriplet.hpp"
#include "DGContext2D.hpp"
#include "MeshManager.hpp"

namespace blitzdg {
    class Poisson2DSparseMatrix {
    private:
    std::unique_ptr<CSCMat> OP_, MM_;

    public:
        Poisson2DSparseMatrix(DGContext2D& dg, MeshManager& mshManager);

        const CSCMat& getMM() const { return *OP_; };
        const CSCMat& getOP() const { return *MM_; };

        const boost::python::numpy::ndarray getOP_numpy() const;
        const boost::python::numpy::ndarray getMM_numpy() const;
    };
}
