// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "MeshManager.hpp"
#include "Types.hpp"
#include "LinAlgHelpers.hpp"
#include <igloo/igloo_alt.h>
#include <boost/algorithm/string.hpp>
#include <whereami.h>
#include <string>
#include <vector>

using boost::algorithm::find_all;
using boost::algorithm::join;
using boost::algorithm::replace_last;
using boost::algorithm::trim_right;
using boost::iterator_range;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::unique_ptr;

namespace blitzdg {
    namespace MeshManagerTests {
        using namespace igloo;
        Describe(MeshManager_Object) {
            using find_vector_type = vector<iterator_range<string::iterator>>;

            string PathDelimeter = "/";
            MeshManager* meshManager=nullptr;
            string* ExePath = nullptr;
            unique_ptr<vector<string>> InputPathVec = nullptr;

            void SetUp() {
                meshManager = new MeshManager();
                if (ExePath == nullptr) {
                    // Deal with paths to the test input files.
                    index_type cap = 1024;
                    char * pathBuffer = new char[cap];
                    index_type length = wai_getExecutablePath(pathBuffer, cap, NULL);
                    ExePath = new string();

                    for(index_type i=0; i < length; i++) {
                        *ExePath += pathBuffer[i]; 
                    }
                    trim_right(*ExePath);
                    delete [] pathBuffer;
                }
                InputPathVec =  unique_ptr<vector<string>>(new vector<string>());
                resolveInputPath();
            }
            void TearDown() {
                delete ExePath;
                delete meshManager;
            }

            void resolveInputPath() {
                string path(*ExePath);
                find_vector_type FindVec;

                replace_last(path, ".exe", "");
                replace_last(path, "/bin/test", "");
                find_all( FindVec, path, "\\" );
                if (FindVec.size() > 0) {
                    PathDelimeter = "\\";
                    replace_last(path, "\\bin\\test", "");
                }

                vector<string> pathVec;
                pathVec.push_back(path);
                pathVec.push_back("input");

                *InputPathVec = std::move(pathVec);
            }

            void get_VertexFilePath(string & vertFile) {
                vector<string>& pathVec = *InputPathVec;
                pathVec.push_back("2box.V");
                vertFile = join(pathVec, PathDelimeter);
                pathVec.pop_back();
            }

            void get_EToVFilePath(string & eToVFile) {
                vector<string>& pathVec = *InputPathVec;
                pathVec.push_back("2box.E2V");
                eToVFile = join(pathVec, PathDelimeter);
                pathVec.pop_back();
            }

            void get_MshFile(string& mshFile) {
                vector<string>& pathVec = *InputPathVec;
                pathVec.push_back("coarse_box.msh");
                mshFile = join(pathVec, PathDelimeter);
                pathVec.pop_back();
            }

            It(Reads_Vertex_Files) {
                string vertexFile = "";
                get_VertexFilePath(vertexFile);
                cout << *ExePath << endl;

                cout << "MeshManager Reads Vertex File: " << vertexFile << endl;
                MeshManager & mgr = *meshManager;

                mgr.readVertices(vertexFile);

                Assert::That(mgr.get_NumVerts(), Equals(6));
                Assert::That(mgr.get_Dim(), Equals(2));

                const real_vector_type& verts = mgr.get_Vertices();

                Assert::That(verts(0), Equals(0.0));
                Assert::That(verts(1), Equals(0.0));

                Assert::That(verts(2), Equals(0.5));
                Assert::That(verts(3), Equals(0.0));

                Assert::That(verts(4), Equals(1.0));
                Assert::That(verts(5), Equals(0.0));

                Assert::That(verts(6), Equals(1.0));
                Assert::That(verts(7), Equals(1.0));

                Assert::That(verts(8), Equals(0.5));
                Assert::That(verts(9), Equals(1.0));

                Assert::That(verts(10), Equals(0.0));
                Assert::That(verts(11), Equals(1.0));
            }

            It(Reads_Element_Files) {
                string eToVFile = "";
                get_EToVFilePath(eToVFile);

                cout << "MeshManager Reads Elements Files: " << eToVFile << endl;
                MeshManager & mgr = *meshManager;

                mgr.readElements(eToVFile);

                Assert::That(mgr.get_NumElements(), Equals(2));
                Assert::That(mgr.get_NumFaces(), Equals(4));

                const index_vector_type& elements = mgr.get_Elements();

                Assert::That(elements(0), Equals(0));
                Assert::That(elements(1), Equals(1));
                Assert::That(elements(2), Equals(4));
                Assert::That(elements(3), Equals(5));
                Assert::That(elements(4), Equals(1));
                Assert::That(elements(5), Equals(2));
                Assert::That(elements(6), Equals(3));
                Assert::That(elements(7), Equals(4));
            }

            It(Can_Print_Vertices_And_DoesNotThrow) {
                string vertexFile = "";
                get_VertexFilePath(vertexFile);
                cout << "Can_Print_Vertices_And_DoesNotThrow" << endl;
                cout << "MeshManager Reads Vertex File: " << vertexFile << endl;

                MeshManager & mgr = *meshManager;   
                mgr.readVertices(vertexFile);
                cout << "Vertices:" << endl;
                mgr.printVertices();
            }

            It(Can_Print_Elements_And_DoesNotThrow) {
                string eToVFile = "";
                get_EToVFilePath(eToVFile);
                cout << "Can_Print_Elements_And_DoesNotThrow" << endl;
                cout << "MeshManager Reads Element File: " << eToVFile << endl;

                MeshManager & mgr = *meshManager;
                mgr.readElements(eToVFile);
                cout << endl << "Elements" << endl;
                mgr.printElements();
            }

            It(Can_Read_Gmsh_Mesh) {
                cout << "Can_Read_Gmsh_Mesh" << endl;
                string mshFile = "";
                get_MshFile(mshFile);
                MeshManager& mgr = *meshManager;
                mgr.readMesh(mshFile);

                const index_vector_type& elements = mgr.get_Elements();
                const index_vector_type& bcTable = mgr.get_BCType();
                const real_vector_type& verts = mgr.get_Vertices();
                const index_vector_type& EToE = mgr.get_EToE();
                const index_vector_type& EToF = mgr.get_EToF();

                index_type K = mgr.get_NumElements();
                index_type Nv = mgr.get_NumVerts();

                Assert::That(Nv, Equals(29));
                Assert::That(K, Equals(40));

                real_vector_type expectedVerts(Nv*3);
                index_vector_type expectedBcTable(K*3);
                index_vector_type expectedElements(K*3);
                index_vector_type expectedEToE(K*3);
                index_vector_type expectedEToF(K*3);

                expectedVerts = -1.,-1., 0., 1., -1., 0., 1., 1., 0., -1., 1., 0., -0.5, 1., 0., -2.7528e-12, 1., 0., 0.5, 1., 0., 1., 0.5, 0., 1., 2.7528e-12, 0., 1., -0.5, 0., 0.5, -1., 0., 2.7528e-12, -1., 0., -0.5, -1. , 0., -1., -0.5, 0, -1, -2.7528e-12, 0, -1, 0.5, 0, -5.41875e-18, -1.80625e-18, 0, 0.416667, 0.416667, 0, 0.416667, -0.416667, 0, -0.416667, -0.416667, 0, -0.416667, 0.416667, -0, 3.10552e-13, 0.559524, 0, 0.559524, -3.10552e-13, 0, -3.10552e-13, -0.559524, 0, -0.559524, 3.10437e-13, -0, -0.677083, 0.677083, -0, 0.677083, 0.677083, 0, 0.677083, -0.677083, 0, -0.677083, -0.677083, 0;
                expectedBcTable = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,3,0,0,3,0,0,3,0,0,3,0,0,3,0,0,3,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,3,0,0,3,0,0,3,0,0,0,0,3,0,0,3,0,0,3,0,0,3;
                expectedElements = 7,17,22,10,18,23,4,20,21,13,19,24,15,24,20,9,22,18,12,23,19,6,21,17,4,21,5,7,22,8,10,23,11,13,24,14,5,21,6,8,22,9,11,23,12,14,24,15,16,17,21,16,21,20,16,22,17,16,23,18,16,24,19,16,18,22,16,19,23,16,20,24,7,26,17,10,27,18,13,28,19,4,25,20,6,17,26,9,18,27,12,19,28,15,20,25,0,12,28,1,9,27,2,6,26,3,15,25,0,28,13,1,27,10,2,26,7,3,25,4;
                expectedEToE = 24,18,9,25,19,10,27,17,8,26,20,11,15,23,31,13,21,29,14,22,30,12,16,28,2,12,8,0,13,9,1,14,10,3,15,11,8,7,12,9,5,13,10,6,14,11,4,15,18,7,17,16,2,23,21,0,16,22,1,21,23,3,22,19,5,18,20,6,19,17,4,20,38,28,0,37,29,1,36,30,3,39,31,2,7,24,34,5,25,33,6,26,32,4,27,35,32,30,36,33,29,37,34,28,38,35,31,39,32,26,36,33,25,37,34,24,38,35,27,39;
                expectedEToF = 2,1,0,2,1,0,2,1,0,2,1,0,1,1,0,1,1,0,1,1,0,1,1,0,2,0,2,2,0,2,2,0,2,2,0,2,1,0,2,1,0,2,1,0,2,1,0,2,2,1,0,2,1,0,2,1,0,2,1,0,2,1,0,2,1,0,2,1,0,2,1,0,1,1,0,1,1,0,1,1,0,1,1,0,2,1,1,2,1,1,2,1,1,2,1,1,0,2,0,0,2,0,0,2,0,0,2,0,2,0,2,2,0,2,2,0,2,2,0,2;

                real_vector_type vertErr(Nv*3);
                real_vector_type bcErr(K*3);
                real_vector_type elemErr(K*3);
                real_vector_type EToEErr(K*3);
                real_vector_type EToFErr(K*3);

                vertErr = verts - expectedVerts;
                bcErr = bcTable - expectedBcTable;
                elemErr = elements - expectedElements;
                EToEErr = EToE - expectedEToE;
                EToFErr = EToF - expectedEToF;

                Assert::That(normInf(vertErr), IsLessThan(3.4e-7));
                Assert::That(normInf(elemErr), Equals(0.0));
                Assert::That(normInf(EToEErr), Equals(0.0));
                Assert::That(normInf(EToFErr), Equals(0.0));
                Assert::That(normInf(bcErr),   Equals(0.0));
            }

            It(Can_Partition_A_Mesh) {
                cout << "Can_Partition_A_Mesh" << endl;
                string eToVFile = ""; 
                string vertexFile = "" ;
                get_EToVFilePath(eToVFile);
                get_VertexFilePath(vertexFile);

                MeshManager & mgr = *meshManager;
                mgr.readVertices(vertexFile);
                mgr.readElements(eToVFile);

                cout << "K: " << mgr.get_NumElements() << endl;
                cout << "Nv: " << mgr.get_NumVerts() << endl;
                mgr.partitionMesh(2);

                const index_vector_type& epMap = mgr.get_ElementPartitionMap();
                Assert::That(epMap(0), Equals(1));
                Assert::That(epMap(1), Equals(0));

                const index_vector_type& vpMap = mgr.get_VertexPartitionMap();
                Assert::That(vpMap(0), IsGreaterThan(-1));
                Assert::That(vpMap(1), IsGreaterThan(-1));
                Assert::That(vpMap(2), IsGreaterThan(-1));
                Assert::That(vpMap(3), IsGreaterThan(-1));
                Assert::That(vpMap(4), IsGreaterThan(-1));
                Assert::That(vpMap(5), IsGreaterThan(-1));
            }
        };
    } // namespace MeshManagerTests
} // namespace blitzdg

