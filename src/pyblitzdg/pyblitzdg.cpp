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

    class_<lserk4wrapper>("LSERK4")
        .def_readonly("numStages", new int(5))
        .add_static_property("rk4a", &lserk4wrapper::rk4a)
        .add_static_property("rk4b", &lserk4wrapper::rk4b);

    class_<MeshManager, boost::noncopyable>("MeshManager", init<>())
        .def("readMesh", &MeshManager::readMesh)
        .def("partitionMesh", &MeshManager::partitionMesh)
        .def("buildMesh", &MeshManager::buildMesh)
        .def("setBCType", &MeshManager::set_BCType_numpy)
        .add_property("vertexPartitionMap", &MeshManager::get_VertexPartitionMap_numpy)
        .add_property("elementPartitionMap", &MeshManager::get_ElementPartitionMap_numpy)
        .add_property("numElements", &MeshManager::get_NumElements)
        .add_property("vertices", &MeshManager::get_Vertices_numpy)
        .add_property("elements", &MeshManager::get_Elements_numpy)
        .add_property("bcType", &MeshManager::get_BCType_numpy);


    class_<TriangleNodesProvisioner, boost::noncopyable>("TriangleNodesProvisioner", init<index_type, MeshManager&>())
        .def("buildFilter", &TriangleNodesProvisioner::buildFilter)
        .def("dgContext", &TriangleNodesProvisioner::get_DGContext);

    class_<QuadNodesProvisioner, boost::noncopyable>("QuadNodesProvisioner", init<index_type, MeshManager&>())
        .def("buildFilter", &QuadNodesProvisioner::buildFilter)
        .def("dgContext", &QuadNodesProvisioner::get_DGContext);

    class_<DGContext2D>("DGContext2D", init<>())
        .add_property("numLocalPoints", &DGContext2D::numLocalPoints)
        .add_property("numFacePoints", &DGContext2D::numFacePoints)
        .add_property("numElements", &DGContext2D::numElements)
        .add_property("numFaces", &DGContext2D::numFaces)
        .add_property("filter", &DGContext2D::filter_numpy)
        .add_property("x", &DGContext2D::x_numpy)
        .add_property("y", &DGContext2D::y_numpy)
        .add_property("Fscale", &DGContext2D::fscale_numpy)
        .add_property("J", &DGContext2D::jacobian_numpy)
        .add_property("rx", &DGContext2D::rx_numpy)
        .add_property("ry", &DGContext2D::ry_numpy)
        .add_property("sx", &DGContext2D::sx_numpy)
        .add_property("sy", &DGContext2D::sy_numpy)
        .add_property("nx", &DGContext2D::nx_numpy)
        .add_property("ny", &DGContext2D::ny_numpy)
        .add_property("Dr", &DGContext2D::Dr_numpy)
        .add_property("Ds", &DGContext2D::Ds_numpy)
        .add_property("Lift", &DGContext2D::lift_numpy)
        .add_property("vmapM", &DGContext2D::vmapM_numpy)
        .add_property("vmapP", &DGContext2D::vmapP_numpy)
        .add_property("BCmap", &DGContext2D::bcmap_numpy);
    
    class_<VtkOutputter>("VtkOutputter", init<TriangleNodesProvisioner&>())
        .def(init<QuadNodesProvisioner&>())
        .def("writeFieldToFile", &VtkOutputter::writeFieldToFile_numpy)
        .def("writeFieldsToFiles", &VtkOutputter::writeFieldsToFiles_numpy);

    class_<Poisson2DSparseMatrix, boost::noncopyable>("Poisson2DSparseMatrix", init<DGContext2D&, MeshManager&>())
        .def("buildBcRhs", &Poisson2DSparseMatrix::buildBcRhs_numpy)
        .def("getOP", &Poisson2DSparseMatrix::getOP_numpy)
        .def("getMM", &Poisson2DSparseMatrix::getMM_numpy);
}


