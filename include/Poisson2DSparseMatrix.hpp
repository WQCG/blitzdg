#pragma once

#include <memory>
#include <suitesparse/umfpack.h>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "CSCMatrix.hpp"
#include "SparseTriplet.hpp"
#include "DGContext2D.hpp"
#include "MeshManager.hpp"
#include "Types.hpp"

namespace blitzdg {
    class Poisson2DSparseMatrix {
    private:
    std::unique_ptr<CSCMat> OP_, MM_;
    std::unique_ptr<real_matrix_type> BcRhs_;

    void buildPoissonOperator(DGContext2D& dg, MeshManager& mshManager, const index_vector_type& bcType);

    public:
        Poisson2DSparseMatrix(DGContext2D& dg, MeshManager& mshManager);

        void buildBcRhs(DGContext2D& dg, const MeshManager& mshManager, const real_matrix_type& ubc, const real_matrix_type& qbc, const index_vector_type& bcType);
        const boost::python::numpy::ndarray buildBcRhs_numpy(DGContext2D& dg, const MeshManager& mshManager, const boost::python::numpy::ndarray& ubc, const boost::python::numpy::ndarray& qbc);

        const CSCMat& getMM() const { return *OP_; };
        const CSCMat& getOP() const { return *MM_; };
        const real_matrix_type& getBcRhs() const { return *BcRhs_;}

        const boost::python::numpy::ndarray getOP_numpy() const;
        const boost::python::numpy::ndarray getMM_numpy() const;
        const boost::python::numpy::ndarray getBcRhs_numpy() const;

    };
}
