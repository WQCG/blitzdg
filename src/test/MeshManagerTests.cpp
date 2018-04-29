// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#include <igloo/igloo_alt.h>
#include <MeshManager.hpp>
#include <whereami.h>

using namespace igloo;
using namespace std;
using namespace boost;

namespace MeshManagerTests {

    Describe(MeshManager_Object) {
        typedef vector< iterator_range<string::iterator> > find_vector_type;

        MeshManager * meshManager=nullptr;
        string ExePath;
        string PathDelimeter;

        void SetUp() {
            meshManager = new MeshManager();

            // Deal with paths to the test input files.
            int cap = 1024;
            char * pathBuffer = new char[cap];
            wai_getExecutablePath(pathBuffer, cap, NULL);
            ExePath = string(pathBuffer);

            find_vector_type FindVec;

            PathDelimeter = "/";
            replace_last(ExePath, ".exe", "");
            replace_last(ExePath, "/bin/test", "");
            find_all( FindVec, ExePath, "\\" );
            if (FindVec.size() > 0) {
                PathDelimeter = "\\";
                replace_last(ExePath, "\\bin\\test", "");
            }
        }

        string get_VertexFilePath() {
            std::vector<std::string> pathVec;
            pathVec.push_back(ExePath);
            pathVec.push_back("input");
            pathVec.push_back("2box.V");
            return join(pathVec, PathDelimeter);
        }

        string get_EToVFilePath() {
            std::vector<std::string> pathVec;
            pathVec.push_back(ExePath);
            pathVec.push_back("input");
            pathVec.push_back("2box.E2V");
            return join(pathVec, PathDelimeter);
        }

        void TearDown() {
            delete meshManager;
        }

        It(Reads_Vertex_Files) {

            string vertexFile = get_VertexFilePath();
            cout << "MeshManager Reads Vertex File: " << vertexFile << endl;
            MeshManager & mgr = *meshManager;

            mgr.readVertices(vertexFile);

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
            string eToVFile = get_EToVFilePath();

            cout << "MeshManager Reads Elements Files: " << eToVFile << endl;
            MeshManager & mgr = *meshManager;

            mgr.readElements(eToVFile);

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
            string vertexFile = get_VertexFilePath();
            cout << "Can_Print_Vertices_And_DoesNotThrow" << endl;
            MeshManager & mgr = *meshManager;
            mgr.readVertices(vertexFile);
            cout << "Vertices:" << endl;
            mgr.printVertices();
        }

        It(Can_Print_Elements_And_DoesNotThrow) {
            string eToVFile = get_EToVFilePath();
            cout << "Can_Print_Elements_And_DoesNotThrow" << endl;
            MeshManager & mgr = *meshManager;
            mgr.readElements(eToVFile);
            cout << endl << "Elements" << endl;
            mgr.printElements();
        }

       It(Can_Partition_A_Mesh) {
            cout << "Can_Partition_A_Mesh" << endl;
            string eToVFile = get_EToVFilePath();
            string vertexFile = get_VertexFilePath();

            MeshManager & mgr = *meshManager;
            mgr.readVertices(vertexFile);
            mgr.readElements(eToVFile);

            cout << "K: " << mgr.get_NumElements() << endl;
            cout << "Nv: " << mgr.get_NumVerts() << endl;
            mgr.partitionMesh(2);

            int * & epMap = mgr.get_ElementPartitionMap();
            Assert::That(epMap[0], Equals(1));
            Assert::That(epMap[1], Equals(0));

			int * & vpMap = mgr.get_VertexPartitionMap();
            Assert::That(vpMap[0], IsGreaterThan(-1));
            Assert::That(vpMap[1], IsGreaterThan(-1));
            Assert::That(vpMap[2], IsGreaterThan(-1));
            Assert::That(vpMap[3], IsGreaterThan(-1));
            Assert::That(vpMap[4], IsGreaterThan(-1));
            Assert::That(vpMap[5], IsGreaterThan(-1));
        } 
    };
}