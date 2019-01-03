// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "MeshManager.hpp"
#include "CSVFileReader.hpp"
#include "Types.hpp"
#include "CSCMatrix.hpp"
#include <boost/algorithm/string.hpp>
#include <metis.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

using boost::algorithm::is_any_of;
using boost::algorithm::split;
using boost::algorithm::trim;
using boost::token_compress_on;
using std::cout;
using std::endl;
using std::getline;
using std::ifstream;
using std::istringstream;
using std::runtime_error;
using std::string;
using std::to_string;
using std::unique_ptr;
using std::vector;
using std::stoi;
using std::abs;

namespace blitzdg {
    namespace {
        template<typename T>
        void printArray(const vector_type<T>& arr, index_type numRows, index_type numCols) {
            for(index_type i = 0; i < numRows; ++i) {
                for (index_type j = 0; j < numCols; ++j) {
                    cout << arr(MeshManager::get_Index(i, j, numCols)) << " ";
                }
                cout << endl;
            }
        }
    } // anonymous namespace

    const real_type MeshManager::NodeTol = 1.e-10;

    MeshManager::MeshManager() 
        :  Dim{ 0 }, NumVerts{ 0 }, NumFaces{}, NumElements{ 0 },
        CsvDelimiters{ "\t " }, Vert{ nullptr }, EToV{ nullptr },
        ElementPartitionMap{ nullptr }, VertexPartitionMap{ nullptr }
    {}


    std::vector<index_type> MeshManager::parseElem(const std::vector<string>& input) {

        std::vector<index_type> elems;
        elems.reserve(input.size());

        for (const auto &s : input)
            elems.push_back(std::stoi(s));

        return elems;
    }

    void MeshManager::readMesh(const string& gmshInputFile) {
        // gmsh .msh files are csv's with special headers/separators that begin
        // with either "'s or $'s.

        CSVFileReader csvReader(gmshInputFile);

        string line;
        csvReader.readLine(line);

        // parse headers
        if (line != "$MeshFormat")
            throw runtime_error("Missing $MeshFormat header in .msh file!");

        csvReader.setNumCols(3);
        float vers;
        int fileType;
        int floatSize;
        csvReader.parseRowValues(vers, fileType, floatSize);

        if (vers < 2.0 || vers >= 3.0)
            throw runtime_error("Unsupported Gmsh version. Only 2.x is currently supported.");

        if(fileType != 0)
            throw runtime_error("Only ASCII-type Gmsh formats are supported.");

        if(floatSize != 8)
            throw runtime_error("Only 8-byte reals in Gmsh files are supported!");

        // This should skip the '$EndMeshFormat' line.
        csvReader.skipLines(1);

        csvReader.readLine(line);
        if (line != "$Nodes")
            throw runtime_error("Unexpected line marker in .msh file! Expected '$Nodes' but was:" + line + ".");
        
        csvReader.setNumCols(1);
        csvReader.parseRowValues(NumVerts);

        csvReader.setNumCols(4);

        // Allocate storage for Vertices.
        Vert = real_vec_smart_ptr(new real_vector_type(NumVerts*3));

        real_vector_type& Vref = *Vert;

        // Parse in the vertex coordinates. Gmsh format always has 
        // z-coordinate (though can be identically zero).
        for(index_type i=0; i < NumVerts; ++i) {
            index_type nodeNum;
            real_type vx, vy, vz;

            csvReader.parseRowValues(nodeNum, vx, vy, vz);
            Vref((nodeNum-1)*3)   = vx;
            Vref((nodeNum-1)*3+1) = vy;
            Vref((nodeNum-1)*3+2) = vz;
        }

        // This should skip the '$EndNodes' line.
        csvReader.skipLines(1);
        
        csvReader.readLine(line);
        if (line != "$Elements")
            throw runtime_error("Unexpected line marker in .msh file! Expected '$Elements' but was:" + line + ".");

        csvReader.setNumCols(1);
        index_type numElementRows = 0;
        csvReader.parseRowValues(numElementRows);

        // It is now assumed that Gmsh gives 'point' elements as Type 15,
        // 'Line' elements as Type 1, Triangles as Type 2, Quadrangles as Type 3
        // and so on as described at http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php.

        std::vector<string> elementInfo;

        std::vector<index_type> points;
        std::vector<std::vector<index_type>> lines;
        std::vector<std::vector<index_type>> tris;
        std::vector<std::vector<index_type>> quads;

        lines.reserve(numElementRows);
        tris.reserve(numElementRows);
        quads.reserve(numElementRows);

        for(index_type i =0; i < numElementRows; ++i) {
            csvReader.readLine(line);
            csvReader.tokenizeLine(line, elementInfo);

            index_type numCols = static_cast<index_type>(elementInfo.size());

            index_type elemType   = stoi(elementInfo[1]);
            index_type numTags    = stoi(elementInfo[2]);

            index_type numLocalVerts = numCols - numTags - 3;

            if (numLocalVerts == 1) {
                if (elemType != 15)
                    throw runtime_error("Incorrect Element Type for point element!");
                points.push_back(stoi(elementInfo[3]));
            }

            if (numLocalVerts == 2) {
                if (elemType !=1)
                    throw runtime_error("Incorrect Element Type for line element!");
                
                lines.push_back(parseElem(elementInfo));
            }

            if (numLocalVerts == 3) {
                if (elemType != 2)
                    throw runtime_error("Incorrect Element Type for triangle element!");

                tris.push_back(parseElem(elementInfo));
            }

            if (numLocalVerts == 4) {
                if (elemType != 3)
                    throw runtime_error("Incorrect Element Type for quadrangle element!");
                
                quads.push_back(parseElem(elementInfo));
            }
        }

        if (!quads.empty())
            throw runtime_error("Quadrangle elements currently not supported by blitzdg!");

        lines.shrink_to_fit();
        tris.shrink_to_fit();

        NumFaces = 3;

        // Allocate storage EToV and BC Table.
        // Note: we are doing this here as opposed to in the initializer list,
        // since prior to this point we did not know how many elements there are.
        index_type K = static_cast<index_type>(tris.size());
        EToV = index_vec_smart_ptr(new index_vector_type(K*3));
        BCType = index_vec_smart_ptr(new index_vector_type(K*3));
        EToE = index_vec_smart_ptr(new index_vector_type(K*3));
        EToF = index_vec_smart_ptr(new index_vector_type(K*3));

        index_vector_type& E2V = *EToV;

        NumElements = K;
        for (index_type k=0; k < K; ++k) {
            // Subtract one to go from 1-based (Gmsh) to 0-based (us).
            E2V(NumFaces*k)   = tris[k][5] - 1;
            E2V(NumFaces*k+1) = tris[k][6] - 1;
            E2V(NumFaces*k+2) = tris[k][7] - 1;
        }

        for (index_type k=0; k < K; ++k) {
            // Enforce counter-clockwise ordering of vertices in EToV table.
            real_type ax = Vref(E2V(NumFaces*k)*NumFaces),   ay = Vref(E2V(NumFaces*k)*NumFaces+1);
            real_type bx = Vref(E2V(NumFaces*k+1)*NumFaces), by = Vref(E2V(NumFaces*k+1)*NumFaces+1);
            real_type cx = Vref(E2V(NumFaces*k+2)*NumFaces), cy = Vref(E2V(NumFaces*k+2)*NumFaces+1);

            real_type det = (ax-cx)*(by-cy) - (bx-cx)*(ay-cy);

            if (det < 0) {
                using std::swap;
                swap(E2V(NumFaces*k+1), E2V(NumFaces*k+2));
            }
        }

        buildBCTable(lines);
        buildConnectivity();
    }

    void MeshManager::buildBCTable(index_type tagNumber) {
        index_vector_type& E2E = *EToE;
        index_vector_type& BCTable = *BCType;

        blitz::firstIndex ii;
        BCTable = 0*ii;
        for (index_type k=0; k < NumElements; ++k) {
            for (index_type f=0; f < 3; ++f) {
                index_type k2 = E2E(k*NumFaces + f);

                // self-referencing elements are boundary elements.
                if (k2 == k)
                    BCTable(3*k + f) = tagNumber;
            }
        }
    }

    void MeshManager::buildBCTable(std::vector<std::vector<index_type>>& edges) {
        index_vector_type& E2V = *EToV;
        real_vector_type& V = *Vert;
        index_vector_type& BCTable = *BCType;

        blitz::firstIndex ii;
        BCTable = 0*ii;
        for (index_type k=0; k < NumElements; ++k) {
            for (index_type f=0; f < 3; ++f) {

                index_type v1 = E2V(k*3 + f);
                index_type v2 = E2V(k*3 + ((f + 1) % 3));

                real_type v1x = V(v1*3);
                real_type v1y = V(v1*3 + 1);

                real_type v2x = V(v2*3);
                real_type v2y = V(v2*3 + 1);

                real_type midx = 0.5*(v1x + v2x);
                real_type midy = 0.5*(v1y + v2y);

                // check if midpoint of the triangle edge
                // lies on a boundary edge, if it does,
                // then the current face is on the boundary.
                for (index_type edge=0; edge < static_cast<index_type>(edges.size()); ++edge) {

                    // Look up node numbers from 'lines' table.
                    // (make node references 1-based).
                    index_type n1 = edges[edge][5] - 1;
                    index_type n2 = edges[edge][6] - 1;

                    index_type bcType = edges[edge][3];
                    // Assign default to something non-zero. 3 -> Wall, in the parlance of NUDG.
                    if (bcType == 0) {
                        bcType = 3;
                    }

                    real_type x1 = V(n1*3), y1 = V(n1*3+1);
                    real_type x2 = V(n2*3), y2 = V(n2*3+1);

                    // Eqn of line: y - y0 = (y2-y1)/(x2-x1)*(x-x0)
                    // => (midy-y2)*(x2-x1) - (y2-y1)(midx-x2) = 0
                    if (std::abs((y2-midy)*(x2-x1) - (y2-y1)*(x2-midx)) < NodeTol) {
                        BCTable(3*k + f) = bcType;
                        break;
                    }
                }
            }
        }
    }

    void MeshManager::buildConnectivity() {
        index_type totalFaces = NumFaces * NumElements;
        index_type numVertices = NumVerts;
        index_type localVertNum[3][2] = {{0, 1}, {1, 2}, {2, 0}};

        // Build global face-to-vertex array. Note: we actually
        // build its transpose for easy column access.
        CSCMat VToF(numVertices, totalFaces, 2*totalFaces);
        const index_vector_type& E2V = *EToV;
        index_type globalFaceNum = 0;
        index_type nnz=0;
        for (index_type k = 0; k < NumElements; ++k) {
            for (index_type f = 0; f < NumFaces; ++f) {
                VToF.colPtrs(globalFaceNum) = nnz; 
                index_type v1 = localVertNum[f][0];
                index_type v2 = localVertNum[f][1];

                index_type v1Global = E2V(k*NumFaces + v1);
                index_type v2Global = E2V(k*NumFaces + v2);

                // Entry for v1Global.
                VToF.rowInds(nnz) = v1Global;
                VToF.elems(nnz) = 1.0;
                ++nnz;

                // Entry for v2Global
                VToF.rowInds(nnz) = v2Global;
                VToF.elems(nnz) = 1.0;
                ++nnz;

                ++globalFaceNum;
            }
        }
        VToF.colPtrs(totalFaces) = nnz;
        CSCMat FToF = multiply(transpose(VToF), VToF);

        // Count the number of face-to-face connections.
        index_type connectionsCount = 0;
        for (index_type j = 0; j < totalFaces; ++j) {
            for (index_type k = FToF.colPtrs(j); k < FToF.colPtrs(j + 1); ++k) {
                index_type i = FToF.rowInds(k);
                if (i != j && std::abs(FToF.elems(k) - 2.0) < NodeTol) {
                    ++connectionsCount;
                }
            }
        }

        index_vector_type e1(connectionsCount);
        index_vector_type e2(connectionsCount);
        index_vector_type f1(connectionsCount);
        index_vector_type f2(connectionsCount);

        f1 = 0;
        f2 = 0;
        e1 = 0;
        e2 = 0;

        // store the faces in each connection.
        index_type c = 0;
        for (index_type j = 0; j < totalFaces; ++j) {
            for (index_type k = FToF.colPtrs(j); k < FToF.colPtrs(j + 1); ++k) {
                index_type i = FToF.rowInds(k);
                if (i != j && std::abs(FToF.elems(k) - 2.0) < NodeTol) {
                    f1(c) = i;
                    f2(c) = j;
                    ++c;
                }
            }
        }

        // Convert from global face number to element number with local face number.
        e1 = f1 / NumFaces;
        f1 = (f1 % NumFaces);
        e2 = f2 / NumFaces;
        f2 = (f2 % NumFaces);

        // Build connectivity matrices.
        index_vector_type& E2E = *EToE;
        index_vector_type& E2F = *EToF;

        // Populate arrays with self-connection defaults.
        for (index_type k = 0; k < NumElements; ++k) {
            for (index_type f = 0; f < NumFaces; ++f) {
                E2E(k*NumFaces + f) = k;
                E2F(k*NumFaces + f) = f;
            }
        }

        for (index_type i = 0; i < connectionsCount; ++i) {
            index_type ee1 = e1(i);
            index_type ee2 = e2(i);
            index_type ff1 = f1(i);
            index_type ff2 = f2(i);
            E2E(ee1*NumFaces + ff1) = ee2;
            E2F(ee1*NumFaces + ff1) = ff2;
        }
    }

    void MeshManager::partitionMesh(index_type numPartitions) {
        index_vector_type eptr(NumElements + 1);
        index_type objval = 0;
        index_type NE = NumElements;
        index_type NV = NumVerts;

        // set up mesh partitioning options
        idx_t metisOptions[METIS_NOPTIONS];
        METIS_SetDefaultOptions(metisOptions);

        metisOptions[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
        metisOptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // total communication volume minimization.
        metisOptions[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
        metisOptions[METIS_OPTION_IPTYPE] = METIS_IPTYPE_METISRB;
        metisOptions[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY;
        metisOptions[METIS_OPTION_NCUTS] = 1;
        metisOptions[METIS_OPTION_NITER] = 10; // default
        metisOptions[METIS_OPTION_SEED] = -1; // random number seed.
        metisOptions[METIS_OPTION_UFACTOR] = 30; // max load imbalance of "30" -> 1.03
        metisOptions[METIS_OPTION_NUMBERING] = 0; // 0-based numbering.
        metisOptions[METIS_OPTION_DBGLVL] = METIS_DBG_INFO; //debug level='info'. 0 for nothing.
        metisOptions[METIS_OPTION_MINCONN] = 1;
        metisOptions[METIS_OPTION_NOOUTPUT] = 0;
        metisOptions[METIS_OPTION_CONTIG] = 1;

        // output arrays
        ElementPartitionMap.reset(new index_vector_type(NumElements));
        VertexPartitionMap.reset(new index_vector_type(NumVerts));
        *ElementPartitionMap = 0;
        *VertexPartitionMap = 0;

        // Assume mesh with homogenous element type, then eptr 
        // dictates an equal stride of size NumFaces across EToV array.
        for (index_type i = 0; i <= NumElements; ++i) {
            eptr(i) = NumFaces * i;
        }

        cout << "About to call METIS_PartMeshNodal" << endl;
        index_type result =  METIS_PartMeshNodal(&NE, &NV, eptr.data(), EToV->data(),
            (idx_t*)NULL, (idx_t*)NULL, &numPartitions, (real_t*)NULL, 
            metisOptions, &objval, ElementPartitionMap->data(), 
            VertexPartitionMap->data());

        if (result == METIS_OK)
            cout << "METIS partitioning successful!" << endl;
        else if (result == METIS_ERROR_INPUT)
            cout << "METIS input error!" << endl;
        else if (result == METIS_ERROR_MEMORY)
            cout << "METIS could not allocate the required memory!" << endl;
        else
            cout << "Unknown METIS error: " << result << endl;

        cout << "total communication volume of partition: " << objval << endl;

        cout << "Element partitioning vector: " << endl;
        for (index_type i = 0; i < NumElements; ++i)
            cout << (*ElementPartitionMap)(i) << endl;

        cout << "Vertex partitioning vector: " << endl;
        for (index_type i = 0; i < NumVerts; ++i)
            cout << (*VertexPartitionMap)(i) << endl;
    }

    void MeshManager::readVertices(const string& vertFile) {
        Vert = csvread<real_type>(vertFile, NumVerts, Dim);
    }

    void MeshManager::readElements(const string& E2VFile) {
        EToV = csvread<index_type>(E2VFile, NumElements, NumFaces);

        BCType = index_vec_smart_ptr(new index_vector_type(NumElements*NumFaces));
        EToE = index_vec_smart_ptr(new index_vector_type(NumElements*NumFaces));
        EToF = index_vec_smart_ptr(new index_vector_type(NumElements*NumFaces));
        buildConnectivity();
        buildBCTable(3);
    }

    void MeshManager::printVertices() const {
        printArray(*Vert, NumVerts, Dim);
    }

    void MeshManager::printElements() const {
        printArray(*EToV, NumElements, NumFaces);
    }

    index_type MeshManager::get_Dim() const {
        return Dim;
    }

    const index_vector_type& MeshManager::get_BCType() const {
        return *BCType;
    }

    index_type MeshManager::get_NumVerts() const {
        return NumVerts;
    }

    index_type MeshManager::get_NumElements() const {
        return NumElements;
    }

    index_type MeshManager::get_NumFaces() const {
        return NumFaces;
    }

    const real_vector_type& MeshManager::get_Vertices() const {
        return *Vert;
    }

    const index_vector_type& MeshManager::get_Elements() const {
        return *EToV;
    }

    const index_vector_type& MeshManager::get_EToE() const {
        return *EToE;
    }

    const index_vector_type& MeshManager::get_EToF() const {
        return *EToF;
    }

    const index_vector_type& MeshManager::get_ElementPartitionMap() const {
        return *ElementPartitionMap;
    }

    const index_vector_type& MeshManager::get_VertexPartitionMap() const {
        return *VertexPartitionMap;
    }
} // namespace blitzdg
