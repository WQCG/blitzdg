// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file pyblitzdg.cpp
 * @brief Python bindings for blitzdg's public API.
 */

#include <boost/python.hpp>
#include "Nodes1DProvisioner.hpp"
#include "VtkOutputter.hpp"
#include "Types.hpp"
#include "LSERK4.hpp"
#include "MeshManager.hpp"
#include "NodesProvisioner2Dbase.hpp"
#include "TriangleNodesProvisioner.hpp"
#include "QuadNodesProvisioner.hpp"
#include "Poisson2DSparseMatrix.hpp"
#include "GaussFaceContext2D.hpp"
#include "CubatureContext2D.hpp"
#include "VandermondeBuilders.hpp"

using namespace boost::python;
using namespace boost::python::numpy;

namespace blitzdg {
#ifdef _WIN32
    class __declspec(dllexport) lserk4wrapper {
#else
    class lserk4wrapper {
#endif
        public:

        static boost::python::list rk4a() {

            boost::python::list a;
            for (int i = 0; i < 5; ++i) {
                a.append(LSERK4::rk4a[i]);
            }
            return a;
        }

        static boost::python::list rk4b() {
            boost::python::list a;
            for (int i = 0; i < 5; ++i) {
                a.append(LSERK4::rk4b[i]);
            }
            return a;
        }
    };
}

BOOST_PYTHON_MODULE(pyblitzdg)
{
    using namespace boost::python;
    using namespace blitzdg;

    boost::python::numpy::initialize();

    class_<Nodes1DProvisioner, boost::noncopyable>("Nodes1DProvisioner", "Class for building a one-dimensional DG mesh from line segments discretized using Legendre-Gauss-Lobatto polynomial interpolation nodes.", init<index_type, index_type, real_type, real_type>(args("self", "NOrder", "K", "xLeft", "xRight")))
        .def("buildNodes", &Nodes1DProvisioner::buildNodes, "Builds the LGL nodes specified by the parameters passed into the constructor."
            , args("self"))
        .def("computeJacobian", &Nodes1DProvisioner::computeJacobian, "Computes the geometric Jacobian of the mapping to the standard element."
            , args("self"))
        .add_property("numLocalPoints", &Nodes1DProvisioner::get_NumLocalPoints, "Property containing the number of points per element.")
        .add_property("xGrid", &Nodes1DProvisioner::get_xGrid_numpy, "Property containing the physical or 'x'-grid, a subset of the real line.")
        .add_property("Dr", &Nodes1DProvisioner::get_Dr_numpy, "Property containing the differentiation matrix for the standard element, a discretization of coordinate r in [-1, 1].")
        .add_property("rx", &Nodes1DProvisioner::get_rx_numpy, "Property containing the derivative of the coordinate r with respect to x.")
        .add_property("Fscale", &Nodes1DProvisioner::get_Fscale_numpy, "Property containing the geometric scaling factor for boundary (i.e., end-point) integrals, i.e., F = 1/Jacobian evaluated along the faces of all elements.")
        .add_property("Lift", &Nodes1DProvisioner::get_Lift_numpy, "Property containing the linear lifting operator to transform from face nodes to elemental nodes left-multiplied by the standard element's inverse mass matrix.")
        .add_property("vmapM", &Nodes1DProvisioner::get_vmapM_numpy, "Property containing the volume-to-surface maps for interior (-) traces.")
        .add_property("vmapP", &Nodes1DProvisioner::get_vmapP_numpy, "Property containing the volume-to-surface maps for exterior (+) traces.")
        .add_property("mapI", &Nodes1DProvisioner::get_mapI, "Property containing the surface index of the Inflow boundary.")
        .add_property("mapO", &Nodes1DProvisioner::get_mapO, "Property containing the surface index of the Outflow boundary.")
        .add_property("nx", &Nodes1DProvisioner::get_nx_numpy, "Property containing the unit outward-pointing normal along elemental surface boundaries.");

    class_<VandermondeBuilders>("VandermondeBuilder")
        .def("buildVandermondeMatrix", &VandermondeBuilders::buildVandermondeMatrix_numpy, "Build Generalized Vandermonde matrix for a numpy 1D array of input points, r.", args("self", "r"));
        //.def("buildGradVandermonde", &VandermondeBuilders::buildGradVandermonde_numpy, "Build Gradient of Generalized Vandermonde matrix for a numpy 1D array of input points, r.", args("self", "r"))

    class_<lserk4wrapper>("LSERK4")
        .def_readonly("numStages", new int(5), "The number of stages in the LSERK4 time-stepper (5).")
        .add_static_property("rk4a", &lserk4wrapper::rk4a, "The a_n coefficients in the LSERK4 time-stepper, where n=0,...,4.")
        .add_static_property("rk4b", &lserk4wrapper::rk4b, "The b_n coefficients in the LSERK4 time-stepper, where n=0,...,4.");

    class_<MeshManager, boost::noncopyable>("MeshManager", init<>())
        .def("readMesh", &MeshManager::readMesh_python, "Method to read in a gmsh input file. Currently only .msh files of version 2.x are supported.")
        .def("partitionMesh", &MeshManager::partitionMesh, "Method to partition the mesh into 'numPartitions' subdomains.")
        .def("buildMesh", &MeshManager::buildMesh, "Method to build a mesh from numpy arrays, 'EToV' (List of elements as an Nv x K 2D array), and 'Vert' (List of vertices as an Nvx3 array), where Nv is the number of vertices and K is the number of elements.")
        .def("setBCType", &MeshManager::set_BCType_numpy, "Method to set the bcType 2D numpy array (NumFaces x NumElements). Useful when rewriting the bcType array for different types of boundary conditions.")
        .add_property("vertexPartitionMap", &MeshManager::get_VertexPartitionMap_numpy, "Property containing the vertex partion map. Can be called after the partitionMesh method has been called on a mesh.")
        .add_property("elementPartitionMap", &MeshManager::get_ElementPartitionMap_numpy, "Property ontaining the element partition map. Can be called after the partitionMesh method has been called on a mesh.")
        .add_property("numElements", &MeshManager::get_NumElements, "Property containing the number of elements in the mesh.")
        .add_property("vertices", &MeshManager::get_Vertices_numpy, "Property containing the Nv x 3 numpy array of vertices in the mesh.")
        .add_property("elements", &MeshManager::get_Elements_numpy, "Property containing the NumElements x NumFaces array of elements.")
        .add_property("bcType", &MeshManager::get_BCType_numpy, "Property containing the NumElements x NumFaces array of boundary condition types. Values of 0 corresponds to interior (non-boundary) faces.");

    class_<TriangleNodesProvisioner, boost::noncopyable>("TriangleNodesProvisioner", "Class for building a two-dimensional triangular DG nodal mesh from triangular mesh discretized using the Hesthaven & Warburton explicit construction for triangular polynomial interpolation nodes.", init<index_type, MeshManager&>(args("self", "NOrder", "MeshManager")))
        .def("buildFilter", &TriangleNodesProvisioner::buildFilter, "Method to construct a polynomial spectral filter that can be retrieved from the dgContext.")
        .def("buildGaussFaceNodes", &TriangleNodesProvisioner::buildGaussFaceNodes, "Method to construct a Gaussian quadrature mesh along faces of order NGauss", args("self", "NGauss (integer)"))
        .def("buildCubatureVolumeMesh", &TriangleNodesProvisioner::buildCubatureVolumeMesh, "Method to construct a 2D cubature mesh on each element.", args("self", "NCubature (integer)"))
        .def("dgContext", &TriangleNodesProvisioner::get_DGContext, "Read-only property containing the 2D triangular DG Context containing all the data necessary for a spatial discretization of the domain of interest.");

    class_<QuadNodesProvisioner, boost::noncopyable>("QuadNodesProvisioner", init<index_type, MeshManager&>())
        .def("buildFilter", &QuadNodesProvisioner::buildFilter, "Method to construct a polynomial spectral filter that can be retrieved from the dgContext.")
        .def("dgContext", &QuadNodesProvisioner::get_DGContext, "Read-only property containing the 2D quadrilateral DG Context containing all the data necessary for a spatial discretization of the domain of interest.");
    
    class_<GaussFaceContext2D, boost::noncopyable>("GaussFaceContext2D", init<>())
        .add_property("NGauss", &GaussFaceContext2D::NGauss)
        .add_property("nx", &GaussFaceContext2D::nx_numpy)
        .add_property("ny", &GaussFaceContext2D::ny_numpy)
        .add_property("sJ", &GaussFaceContext2D::sJ_numpy)
        .add_property("J",  &GaussFaceContext2D::Jac_numpy)
        .add_property("rx", &GaussFaceContext2D::rx_numpy)
        .add_property("ry", &GaussFaceContext2D::ry_numpy)
        .add_property("sx", &GaussFaceContext2D::sx_numpy)
        .add_property("sy", &GaussFaceContext2D::sy_numpy)
        .add_property("BCmap", &GaussFaceContext2D::bcMap_numpy)
        .add_property("x", &GaussFaceContext2D::x_numpy)
        .add_property("y", &GaussFaceContext2D::y_numpy)
        .add_property("W", &GaussFaceContext2D::W_numpy);

    class_<CubatureContext2D, boost::noncopyable>("CubatureContext2D", init<>())
        .add_property("NCubature", &CubatureContext2D::NCubature)
        .add_property("NumCubaturePoints", &CubatureContext2D::NumCubaturePoints)
        .add_property("r", &CubatureContext2D::r_numpy)
        .add_property("s", &CubatureContext2D::s_numpy)
        .add_property("w", &CubatureContext2D::w_numpy)
        .add_property("V", &CubatureContext2D::V_numpy)
        .add_property("rx", &CubatureContext2D::rx_numpy)
        .add_property("ry", &CubatureContext2D::ry_numpy)
        .add_property("sx", &CubatureContext2D::sx_numpy)
        .add_property("sy", &CubatureContext2D::sy_numpy)
        .add_property("J", &CubatureContext2D::Jac_numpy)
        .add_property("Dr", &CubatureContext2D::Dr_numpy)
        .add_property("Ds", &CubatureContext2D::Ds_numpy)
        .add_property("MM", &CubatureContext2D::MM_numpy)
        .add_property("MMChol", &CubatureContext2D::MMChol_numpy);

    class_<DGContext2D>("DGContext2D", init<>())
        .add_property("numLocalPoints", &DGContext2D::numLocalPoints, "The number of nodal points in an element.")
        .add_property("numFacePoints", &DGContext2D::numFacePoints, "The number of nodal points along an element's face.")
        .add_property("numElements", &DGContext2D::numElements, "The total number of elements in the mesh.")
        .add_property("numFaces", &DGContext2D::numFaces, "The number of faces in one element.")
        .add_property("filter", &DGContext2D::filter_numpy, "The spectral element-wise filter matrix with dimension NumLocalPoints x NumLocalPoints.")
        .add_property("r", &DGContext2D::r_numpy, "The r-coordinate on the standard element.")
        .add_property("s", &DGContext2D::s_numpy, "The s-coordinate on the standard element.")
        .add_property("x", &DGContext2D::x_numpy, "The physical x-coordinates in a NumLocalPoints x NumElements array.")
        .add_property("y", &DGContext2D::y_numpy, "The physical y-coordinates in a NumLocalPoints x NumElements array.")
        .add_property("Fscale", &DGContext2D::fscale_numpy, "Scaling factor for surface integrals along the faces of each element. Fscale = |normal|/J.")
        .add_property("Fmask", &DGContext2D::fmask_numpy, "Mask that when applied to an element volume, returns the values along faces of the element, e.g., Fx1 = x[Fmask, 0].")
        .add_property("gather", &DGContext2D::gather_numpy, "Index mapping from duplicated-interfaces global node set to unique global node numbering.")
        .add_property("scatter", &DGContext2D::scatter_numpy, "Index mapping from unique global node numbering to standard DG global node numbering. The inverse map of gather.")
        .add_property("J", &DGContext2D::jacobian_numpy, "The Jacobian of the mapping to the standard element: J = xr*ys - xs*yr.")
        .add_property("rx", &DGContext2D::rx_numpy, "The x-derivative of the r-coordinate on each element.")
        .add_property("ry", &DGContext2D::ry_numpy, "The y-derivative of the r-coordinate on each element.")
        .add_property("sx", &DGContext2D::sx_numpy, "The x-derivative of the s-coordinate on each element.")
        .add_property("sy", &DGContext2D::sy_numpy, "The y-derivative of the s-coordinate on each element.")
        .add_property("nx", &DGContext2D::nx_numpy, "x-component of the element-wise unit outward normal.")
        .add_property("ny", &DGContext2D::ny_numpy, "y-component of the element-wise unit outward normal.")
        .add_property("Dr", &DGContext2D::Dr_numpy, "Differentiation matrix with respect to the 'r'-coordinate on the standard 2D element.")
        .add_property("Ds", &DGContext2D::Ds_numpy, "Differentiation matrix with respect to the 's'-coordinate on the standard 2D element.")
        .add_property("Lift", &DGContext2D::lift_numpy, "Matrix operator transforming surface integral contributions to volume contributions. Lift = MassMatrix^{-1} * E, where E contains the edge-mass matrices of the standard element.")
        .add_property("vmapM", &DGContext2D::vmapM_numpy, "Volume-to-Surface index map. Interior (-) traces.")
        .add_property("vmapP", &DGContext2D::vmapP_numpy, "Volume-to-surface index map. Exterior (+) traces.")
        .add_property("BCmap", &DGContext2D::bcmap_numpy, "Dictionary containing a list of boundary condition nodes for each boundary type (key).");
    
    class_<VtkOutputter>("VtkOutputter", "Class for writting output to unstructured visualization toolkit (vtk) files.", init<TriangleNodesProvisioner&>(args("TriangleNodesProvisioner")))
        .def(init<QuadNodesProvisioner&>(args("QuadNodesProvisioner")))
        .def("writeFieldToFile", &VtkOutputter::writeFieldToFile_numpy, "Write a field to a .vtu file for a given time-index number.")
        .def("writeFieldsToFiles", &VtkOutputter::writeFieldsToFiles_numpy, "Write a dictionary of fields to a set of .vtu files for a given time-index number.");

    class_<Poisson2DSparseMatrix, boost::noncopyable>("Poisson2DSparseMatrix", init<DGContext2D&, MeshManager&>(args("DGContext2D", "MeshManager")))
        .def("buildBcRhs", &Poisson2DSparseMatrix::buildBcRhs_numpy, "Build boundary conditions contribution to right-hand side of the linear system for the Poisson problem.")
        .def("getOP", &Poisson2DSparseMatrix::getOP_numpy, "Read-only property to retrieve the DG-discretized sparse 2D Poisson operator")
        .def("getMM", &Poisson2DSparseMatrix::getMM_numpy, "Read-only property containing the DG-discretize sparse 2D Mass Matrix.");
}


