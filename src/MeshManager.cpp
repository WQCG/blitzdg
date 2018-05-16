// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "MeshManager.hpp"
#include "Types.hpp"
#include <boost/algorithm/string.hpp>
#include <metis.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using boost::algorithm::is_any_of;
using boost::algorithm::split;
using boost::algorithm::trim;
using boost::token_compress_on;
using std::cout;
using std::endl;
using std::getline;
using std::ifstream;
using std::string;
using std::vector;

namespace blitzdg {
    /**
     * Constructor.
     */
    MeshManager::MeshManager() {
        Vert = nullptr;
        CsvDelimeters = " \t";
        NumVerts = 0;
        NumElements = 0;
        Dim = 0;

        ElementPartitionMap = nullptr;
        VertexPartitionMap = nullptr;
    }

    /**
     * Convert (row,col) index to integer index into a contiguous block of memory.
     */
    index_type MeshManager::get_Index(index_type row, index_type col, index_type numCols) {
        return col + row*numCols;
    }

    /**
     * Read a .msh file that was generated with Gmsh.
     */
    void MeshManager::readMesh(string gmshInputFile) {
        throw("Not implemented!");
    }

    template<typename T>
    void MeshManager::printArray(T * & arr, index_type numRows, index_type numCols) {
        for(index_type i=0; i < numRows; i++) {
            for (index_type j=0; j < numCols; j++) {
                cout << arr[get_Index(i, j, numCols)] << " ";
            }
            cout << endl;
        }
    }

        template<typename T>
        void  MeshManager::readCsvFile(string csvFile, string delimiters, T * & result, index_type * & dims) {
        ifstream fileStream(csvFile);

        string line;

        vector<string> splitVec;
        index_type numLines = 0;

        index_type numCols = -1;
        while(getline(fileStream, line)) {
            // Take first line as source of truth for number of columns.
            if (numLines == 0) {
                trim(line);
                split( splitVec, line, is_any_of(delimiters), token_compress_on );
                numCols = splitVec.size();
            }
            numLines++;
        }

        dims = new index_type[2];
        dims[0] = numLines;
        dims[1] = numCols;
        // roll-back stream.
        fileStream.clear();
        fileStream.seekg(0, std::ios::beg);

        result = new T[numLines*numCols];
        index_type count = 0;
        while(getline(fileStream, line)) {
            trim(line);
            vector<string> splitVec;
            split( splitVec, line, is_any_of(delimiters), token_compress_on ); 
            
            for(index_type i=0; i < numCols; i++) {
                result[count] = atof(splitVec[i].c_str());
                count++;
            }
        }
        fileStream.close();

    }

    /**
     * Partition a mesh into numPartitions partitions using METIS. Results can be obtained by calling
     * MeshManager.get_ElementPartitionMap() and MeshManager.get_VertexPartitionMap().
     */
    void MeshManager::partitionMesh(index_type numPartitions) {
        index_type * eind = EToV;
        index_type * eptr = new index_type[NumElements+1];
        index_type * objval = new index_type;
        index_type * numPartitionsPtr = new index_type;
        *numPartitionsPtr = numPartitions;

        // set up mesh partitioning options
        index_type * metisOptions = new index_type[METIS_NOPTIONS];
        metisOptions[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
        metisOptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // total communication volume minimization.
        metisOptions[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
        metisOptions[METIS_OPTION_IPTYPE] = METIS_IPTYPE_METISRB;
        metisOptions[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY;
        metisOptions[METIS_OPTION_NCUTS] = 1;
        metisOptions[METIS_OPTION_NITER] = 10; // default
        metisOptions[METIS_OPTION_SEED] = -1; // random number seed.
        metisOptions[METIS_OPTION_UFACTOR] = 1.030; // max load imbalance of 1.03
        metisOptions[METIS_OPTION_NUMBERING] = 0; // 0-based numbering.
        metisOptions[METIS_OPTION_DBGLVL] = METIS_DBG_INFO; //debug level='info'. 0 for nothing.
        metisOptions[METIS_OPTION_MINCONN] = 1;
        metisOptions[METIS_OPTION_NOOUTPUT] = 0;
        metisOptions[METIS_OPTION_CONTIG] = 1;

    // output arrays
        index_type * epart = new index_type[NumElements];
        index_type * npart = new index_type[NumVerts];

        for (index_type i=0; i < NumElements; i++)
            epart[i] = 0;

        for (index_type i=0; i < NumVerts; i++)
            npart[i] = 0;

        // Assume mesh with homogenous element type, then eptr 
        // dictates an equal stride of size ElementType across EToV array.
        for (index_type i=0; i <= NumElements; i++) {
            eptr[i] = ElementType*i;
            cout << eptr[i] << endl;
        }

        *objval = 0;

        index_type * NE = new index_type;
        index_type * NV = new index_type;

        *NE = NumElements;
        *NV = NumVerts;

        index_type * ncommon = new index_type;
        *ncommon = 1;

        cout << "About to call METIS_PartMeshNodal" << endl;
        index_type result =  METIS_PartMeshNodal( NE, NV, eptr, eind, (idx_t*)NULL, (idx_t*)NULL,
                        numPartitionsPtr, (real_t*)NULL, metisOptions, objval, epart, npart);

        if (result == METIS_OK)
            cout << "METIS partitioning successful!" << endl;
        else if (result == METIS_ERROR_INPUT)
            cout << "METIS input error!" << endl;
        else if (result == METIS_ERROR_MEMORY)
            cout << "METIS could not allocate the required memory!" << endl;
        else
            cout << "Unknown METIS error: " << result << endl;

        cout << "total communication volume of partition: " << *objval << endl;

        cout << "Element partitioning vector: " << endl;
        for (index_type i=0; i<NumElements; i++)
            cout << epart[i] << endl;

        cout << "Vertex partitioning vector: " << endl;
        for (index_type i=0; i<NumVerts; i++)
            cout << npart[i] << endl;

        ElementPartitionMap = epart;
        VertexPartitionMap = npart;

        delete NE;
        delete NV;
        delete ncommon;
        delete objval;
        delete numPartitionsPtr;
        delete[] eptr;
        delete[] metisOptions;
    }

    /**
     * Read a list of vertices from a file. x-, y-, (and z-) are coordinates delimited by spaces, e.g., 0.5 1.0.
     */
    void MeshManager::readVertices(string vertFile) {
        index_type * dims;
        readCsvFile<real_type>(vertFile, CsvDelimeters, Vert, dims);
        NumVerts = dims[0];
        Dim = dims[1];
    }

    /**
     * Read a list of elments from a file. Vertex numbers are written in a row and delimited by spaces, e.g., 1 2 3 4
     */
    void MeshManager::readElements(string E2VFile) {
        index_type * dims;
        readCsvFile<index_type>(E2VFile, CsvDelimeters, EToV, dims);
        NumElements = dims[0];
        ElementType = dims[1];
    }

    /**
     * Print the list of vertices to stdout.
     */
    void MeshManager::printVertices() {
        MeshManager::printArray<real_type>(Vert, NumVerts, Dim);
    }

    /**
     * Print the list of elements to stdout.
     */
    void MeshManager::printElements() {
        MeshManager::printArray<index_type>(EToV, NumElements, ElementType);
    }

    /**
     * Returns a reference to the list of vertices, a contiguous block of type real_type.
     */
    real_type* & MeshManager::get_Vertices() {
        return Vert;
    }

    /**
     * Returns the dimension of the vertex data. Usuallly will be 2 or 3.
     */
    index_type MeshManager::get_Dim() {
        return Dim;
    }

    /**
     * Returns the number of vertices.
     */
    index_type MeshManager::get_NumVerts() {
        return NumVerts;
    }

    /**
     * Returns the number of elements.
     */
    index_type MeshManager::get_NumElements() {
        return NumElements;
    }

    /**
     * Returns the type of element. 3 => triangles, 4 => quadrilaterals, etc.
     */
    index_type MeshManager::get_ElementType() {
        return ElementType;
    }

    /**
     * Returns a reference to the list of elements, a contiguous block of type index_type.
     */
    index_type * & MeshManager::get_Elements() {
        return EToV;
    }

    /**
     * Returns a reference to the element partition map, an array of type index_type.
     */
    index_type * & MeshManager::get_ElementPartitionMap() {
        return ElementPartitionMap;
    }

    /**
     * Returns a reference to the vertex partition map, an array of type index_type.
     */
    index_type * & MeshManager::get_VertexPartitionMap() {
        return VertexPartitionMap;
    }

    MeshManager::~MeshManager() {
        if (Vert != nullptr) delete[] Vert;
        if (EToV != nullptr) delete[] EToV;
        if (ElementPartitionMap != nullptr) delete[] ElementPartitionMap;
        if (VertexPartitionMap != nullptr) delete[] VertexPartitionMap;
    }
} // namespace blitzdg
