// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file pyblitzdg.cpp
 * @brief Python bindings for blitzdg's public API.
 */

#include <boost/python.hpp>
#include "Nodes1DProvisioner.hpp"
#include "Types.hpp"
#include "LSERK4.hpp"
#include "MeshManager.hpp"
#include "TriangleNodesProvisioner.hpp"

using namespace boost::python;
using namespace boost::python::numpy;

namespace blitzdg {
    class lserk4wrapper{
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

    class_<Nodes1DProvisioner, boost::noncopyable>("Nodes1DProvisioner", init<index_type, index_type, real_type, real_type>())
        .def("buildNodes", &Nodes1DProvisioner::buildNodes)
        .def("computeJacobian", &Nodes1DProvisioner::computeJacobian)
        .add_property("numLocalPoints", &Nodes1DProvisioner::get_NumLocalPoints)
        .add_property("xGrid", &Nodes1DProvisioner::get_xGrid_numpy)
        .add_property("Dr", &Nodes1DProvisioner::get_Dr_numpy)
        .add_property("rx", &Nodes1DProvisioner::get_rx_numpy)
        .add_property("Fscale", &Nodes1DProvisioner::get_Fscale_numpy)
        .add_property("Lift", &Nodes1DProvisioner::get_Lift_numpy)
        .add_property("vmapM", &Nodes1DProvisioner::get_vmapM_numpy)
        .add_property("vmapP", &Nodes1DProvisioner::get_vmapP_numpy)
        .add_property("mapI", &Nodes1DProvisioner::get_mapI)
        .add_property("mapO", &Nodes1DProvisioner::get_mapO)
        .add_property("nx", &Nodes1DProvisioner::get_nx_numpy);

    class_<lserk4wrapper>("LSERK4")
        .def_readonly("numStages", new int(5))
        .add_static_property("rk4a", &lserk4wrapper::rk4a)
        .add_static_property("rk4b", &lserk4wrapper::rk4b);

    class_<MeshManager, boost::noncopyable>("MeshManager", init<>())
        .def("readMesh", &MeshManager::readMesh)
        .add_property("numElements", &MeshManager::get_NumElements);

    class_<TriangleNodesProvisioner, boost::noncopyable>("TriangleNodesProvisioner", init<index_type, MeshManager&>())
        .def("buildFilter", &TriangleNodesProvisioner::buildFilter)
        .def("dgContext", &TriangleNodesProvisioner::get_DGContext);

    class_<DGContext2D>("DGContext2D", init<>())
        .add_property("numLocalPoints", &DGContext2D::numLocalPoints)
        .add_property("numFacePoints", &DGContext2D::numFacePoints)
        .add_property("numElements", &DGContext2D::numElements)
        .add_property("numFaces", &DGContext2D::numFaces)
        .add_property("filter", &DGContext2D::filter_numpy);
}


