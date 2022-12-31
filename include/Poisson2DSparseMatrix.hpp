#pragma once

#include <memory>
#include <suitesparse/umfpack.h>
#include <boost/python/numpy.hpp>
#include "CSCMatrix.hpp"
#include "CubatureContext2D.hpp"
#include "GaussFaceContext2D.hpp"
#include "SparseTriplet.hpp"
#include "DGContext2D.hpp"
#include "MeshManager.hpp"
#include "Types.hpp"

namespace blitzdg {
    using python_array = boost::python::numpy::ndarray;

    class Poisson2DSparseMatrix {
    private:
    std::unique_ptr<CSCMat> OP_, MM_;
    std::unique_ptr<real_matrix_type> BcRhs_;

    void buildPoissonOperator(DGContext2D& dg, MeshManager& mshManager, const index_vector_type& bcType, index_type bordered, index_type skipDG);
    void buildSEMPoissonOperator(DGContext2D& dg, MeshManager& mshManager, index_type bordered);
    void buildPoissonOperator(DGContext2D& dg, MeshManager& mshManager, const index_vector_type& bcType, GaussFaceContext2D& gCtx, CubatureContext2D& cubCtx);

    public:
        Poisson2DSparseMatrix(DGContext2D& dg, MeshManager& mshManager, index_type bordered, index_type skipDG);
        Poisson2DSparseMatrix(DGContext2D& dg, MeshManager& mshManager, GaussFaceContext2D& gCtx, CubatureContext2D& cubCtx);
        

        void buildBcRhs(DGContext2D& dg, const MeshManager& mshManager, const real_matrix_type& ubc, const real_matrix_type& qbc, const index_vector_type& bcType);
        void buildBcRhs(DGContext2D& dg, const MeshManager& mshManager, const GaussFaceContext2D& gctx, const CubatureContext2D& cubCtx, const real_matrix_type& ubc, const real_matrix_type& qbc, const index_vector_type& bcType);
        const python_array buildBcRhs_numpy(DGContext2D& dg, const MeshManager& mshManager, const python_array& ubc, const python_array& qbc);
        const python_array buildCubatureBcRhs_numpy(DGContext2D& dg, const MeshManager& mshManager, const GaussFaceContext2D& gctx, const CubatureContext2D& cubCtx, const python_array& ubc, const python_array& qbc);

        const CSCMat& getMM() const { return *OP_; };
        const CSCMat& getOP() const { return *MM_; };
        const real_matrix_type& getBcRhs() const { return *BcRhs_;}

        const python_array getOP_numpy() const;
        const python_array getMM_numpy() const;
        const python_array getBcRhs_numpy() const;

    };
}
