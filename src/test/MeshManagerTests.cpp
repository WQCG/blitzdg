// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#include <igloo/igloo_alt.h>
#include <MeshManager.hpp>

using namespace igloo;
using namespace std;

MeshManager * meshManager=nullptr;

namespace MeshManagerTests {
    Describe(MeshManager_Object) {
        void SetUp() {
            meshManager = new MeshManager();
        }

        It(Reads_Vertex_Files) {
            cout << "MeshManager" << endl;
            MeshManager & mgr = *meshManager;

            mgr.readVertices("input/2box.V");

            Assert::That(mgr.get_NumVerts(), Equals(6));
            Assert::That(mgr.get_Dim(), Equals(2));

            double * verts = mgr.get_Vertices();

            Assert::That(verts[0], Equals(0.0));
            Assert::That(verts[1], Equals(0.0));

            Assert::That(verts[2], Equals(0.5));
            Assert::That(verts[3], Equals(0.0));

            Assert::That(verts[4], Equals(1.0));
            Assert::That(verts[5], Equals(0.0));

            Assert::That(verts[6], Equals(1.0));
            Assert::That(verts[7], Equals(1.0));

            Assert::That(verts[8], Equals(0.5));
            Assert::That(verts[9], Equals(1.0));

            Assert::That(verts[10], Equals(0.0));
            Assert::That(verts[11], Equals(1.0));
        }

        It(Reads_Element_Files) {
            cout << "Reads Elements Files" << endl;
            MeshManager & mgr = *meshManager;

            mgr.readElements("input/2box.E2V");

            Assert::That(mgr.get_NumElements(), Equals(2));
            Assert::That(mgr.get_ElementType(), Equals(4));

            int * elements = mgr.get_Elements();

            Assert::That(elements[0], Equals(0));
            Assert::That(elements[1], Equals(1));
            Assert::That(elements[2], Equals(4));
            Assert::That(elements[3], Equals(5));

            Assert::That(elements[4], Equals(1));
            Assert::That(elements[5], Equals(2));
            Assert::That(elements[6], Equals(3));
            Assert::That(elements[7], Equals(4));
        }

        It(Can_Print_Vertices_And_DoesNotThrow) {
            cout << "Can_Print_Vertices_And_DoesNotThrow" << endl;
            MeshManager & mgr = *meshManager;
            mgr.readVertices("input/2box.V");
            cout << endl << "Vertices:" << endl;
            mgr.printVertices();
        }

        It(Can_Print_Elements_And_DoesNotThrow) {
            cout << "Can_Print_Elements_And_DoesNotThrow" << endl;
            MeshManager & mgr = *meshManager;
            mgr.readElements("input/2box.E2V");
            cout << endl << "Elements" << endl;
            mgr.printElements();
        }

       It(Can_Partition_A_Mesh) {
            cout << "Can_Partition_A_Mesh" << endl;
            MeshManager & mgr = *meshManager;
            mgr.readVertices("input/2box.V");
            mgr.readElements("input/2box.E2V");

            mgr.partitionMesh(2);

            int * & epMap = mgr.get_ElementPartitionMap();
            Assert::That(epMap[0], Equals(1));
            Assert::That(epMap[1], Equals(0));

			int * & vpMap = mgr.get_VertexPartitionMap();
            Assert::That(vpMap[0], Equals(1));
            Assert::That(vpMap[1], Equals(1));
            Assert::That(vpMap[2], Equals(0));
            Assert::That(vpMap[3], Equals(0));
            Assert::That(vpMap[4], Equals(0));
            Assert::That(vpMap[5], Equals(1));
        } 
    };
}