// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "TriangleNodesProvisioner.hpp"
#include "CSCMatrix.hpp"
#include "BlitzHelpers.hpp"
#include "DenseMatrixInverter.hpp"
#include "MeshManager.hpp"
#include "Types.hpp"
#include <blitz/array.h>
#include <cmath>
#include <limits>
#include <memory>

using blitz::firstIndex;
using blitz::Range;
using blitz::secondIndex;
using blitz::sum;
using blitz::thirdIndex;
using std::numeric_limits;
using std::unique_ptr;

namespace blitzdg {
    const index_type TriangleNodesProvisioner::NumFaces = 3;
    const real_type TriangleNodesProvisioner::NodeTol = 1.e-5;

    TriangleNodesProvisioner::TriangleNodesProvisioner(index_type _NOrder, index_type _NumElements, const MeshManager * _MeshManager) 
        : NumElements{ _NumElements }, NOrder{ _NOrder },
        NumLocalPoints{ (int)0.5*(_NOrder + 2)*(_NOrder+1) }, 
        xGrid{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), NumElements) },
        rGrid{ new real_vector_type((int)0.5*(_NOrder + 2)*(_NOrder+1)) },
        V{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), (int)0.5*(_NOrder + 2)*(_NOrder+1)) }, 
        Dr{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), (int)0.5*(_NOrder + 2)*(_NOrder+1)) },
        Lift{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), (_NOrder+1)*NumFaces) },
        J{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), NumElements) },
        rx{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), NumElements) },
        nx{ new real_matrix_type((_NOrder+1)*NumFaces, NumElements) },
        Vinv{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), (int)0.5*(_NOrder + 2)*(_NOrder+1)) },
        Fmask{ new index_vector_type((_NOrder+1)*NumFaces) },
        Fx{ new real_matrix_type((_NOrder+1)*NumFaces, NumElements) },
        Fscale{ new real_matrix_type((_NOrder+1)*NumFaces, NumElements) },
        EToV{ new index_matrix_type(NumElements, NumFaces) },
        EToE{ new index_matrix_type(NumElements, NumFaces) },
        EToF{ new index_matrix_type(NumElements, NumFaces) },
        vmapM{ new index_vector_type((_NOrder+1)*NumFaces*NumElements) },
        vmapP{ new index_vector_type((_NOrder+1)*NumFaces*NumElements) },
        Mesh2D { _MeshManager },
        Nodes1D{ new Nodes1DProvisioner(_NOrder, NumElements, -1.0, 1.0) }
    {}

    void TriangleNodesProvisioner::buildVandermondeMatrix() {
    }
}