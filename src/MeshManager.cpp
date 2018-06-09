// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "MeshManager.hpp"
#include "CSVFileReader.hpp"
#include "Types.hpp"
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

    MeshManager::MeshManager() 
        :  Dim{ 0 }, NumVerts{ 0 }, ElementType{}, NumElements{ 0 },
        CsvDelimiters{ "\t " }, Vert{ nullptr }, EToV{ nullptr },
        ElementPartitionMap{ nullptr }, VertexPartitionMap{ nullptr }
    {}

    void MeshManager::readMesh(const string& gmshInputFile) {
        throw runtime_error("MeshManager::readMesh not implemented");
    }

    void MeshManager::partitionMesh(index_type numPartitions) {
        unique_ptr<index_type[]> eptr(new index_type[NumElements + 1]);
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
        // dictates an equal stride of size ElementType across EToV array.
        for (index_type i = 0; i <= NumElements; ++i) {
            eptr[i] = ElementType * i;
        }

        cout << "About to call METIS_PartMeshNodal" << endl;
        index_type result =  METIS_PartMeshNodal(&NE, &NV, eptr.get(), EToV->data(), 
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
        EToV = csvread<index_type>(E2VFile, NumElements, ElementType);
    }

    void MeshManager::printVertices() const {
        printArray(*Vert, NumVerts, Dim);
    }

    void MeshManager::printElements() const {
        printArray(*EToV, NumElements, ElementType);
    }

    index_type MeshManager::get_Dim() const {
        return Dim;
    }

    index_type MeshManager::get_NumVerts() const {
        return NumVerts;
    }

    index_type MeshManager::get_NumElements() const {
        return NumElements;
    }

    index_type MeshManager::get_ElementType() const {
        return ElementType;
    }

    const real_vector_type& MeshManager::get_Vertices() const {
        return *Vert;
    }

    const index_vector_type& MeshManager::get_Elements() const {
        return *EToV;
    }

    const index_vector_type& MeshManager::get_ElementPartitionMap() const {
        return *ElementPartitionMap;
    }

    const index_vector_type& MeshManager::get_VertexPartitionMap() const {
        return *VertexPartitionMap;
    }
} // namespace blitzdg
