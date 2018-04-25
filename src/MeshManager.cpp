// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#include <iostream>
#include <math.h>
#include <MeshManager.hpp>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <metis.h>

using namespace std;
using namespace boost;

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
int MeshManager::get_Index(int row, int col, int numCols) {
    return col + row*numCols;
}

/**
 * Read a .msh file that was generated with Gmsh.
 */
void MeshManager::readMesh(string gmshInputFile) {
    throw("Not implemented!");
}

template<typename T>
void MeshManager::printArray(T * & arr, int numRows, int numCols) {
    for(int i=0; i < numRows; i++) {
		for (int j=0; j < numCols; j++) {
			cout << arr[get_Index(i, j, numCols)] << " ";
		}
		cout << endl;
	}
}

    template<typename T>
    void  MeshManager::readCsvFile(string csvFile, string delimiters, T * & result, int * & dims) {
    ifstream fileStream(csvFile);

    string line("");

    vector<string> splitVec;
    int numLines = 0;

    int numCols = -1;
    while(getline(fileStream, line)) {
        // Take first line as source of truth for number of columns.
        if (numLines == 0) {
            trim(line);
            split( splitVec, line, is_any_of(delimiters), token_compress_on );
            numCols = splitVec.size();
        }
        numLines++;
    }

     dims = new int[2];
     dims[0] = numLines;
     dims[1] = numCols;
    // roll-back stream.
    fileStream.clear();
    fileStream.seekg(0, std::ios::beg);

    result = new T[numLines*numCols];
    int count = 0;
    while(getline(fileStream, line)) {
        trim(line);
        vector<string> splitVec;
        split( splitVec, line, is_any_of(delimiters), token_compress_on ); 
        
        for(int i=0; i < numCols; i++) {
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
void MeshManager::partitionMesh(int numPartitions) {
    int * metisOptions =  NULL;
    int * epart = NULL;
    int * npart = NULL;
    int * objval;

    int * numPartitionsPtr = new int;
    *numPartitionsPtr = 1;

/*
    Options ---------------------------------------------------------------------
 ptype=kway, objtype=cut, ctype=shem, rtype=greedy, iptype=metisrb
 dbglvl=0, ufactor=1.030, minconn=NO, contig=NO, nooutput=NO
 seed=-1, niter=10, ncuts=1
 gtype=dual, ncommon=1, niter=10, ncuts=1

Direct k-way Partitioning ---------------------------------------------------
 - Edgecut: 4. */

    // set up mesh partitioning options
    metisOptions = new int[METIS_NOPTIONS];
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
    metisOptions[METIS_OPTION_DBGLVL] = 0; //debug level='info'. 0 for nothing.
    metisOptions[METIS_OPTION_MINCONN] = 0;
    metisOptions[METIS_OPTION_NOOUTPUT] = 0;

    cout << "METIS_NOPTIONS: " << METIS_NOPTIONS << endl;
    cout << "METIS_OPTION_NUMBERING: " << METIS_OPTION_NUMBERING << endl;
    

    int * eind = new int [NumVerts*NumElements];
    int * eptr = new int[NumElements+1];

    eind[0] = 0;
    eind[1] = 1;
    eind[2] = 2;
    eind[3] = 3;
    eind[4] = 4;
    eind[5] = 5;

   // output arrays
    epart = new int[NumElements];
    npart = new int[6];

    epart[0] =0;
    epart[1] =0;

    npart[0] = 0;
    npart[1] = 0;
    npart[2] = 0;
    npart[3] = 0;
    npart[4] = 0;
    npart[5] = 0;

    objval = new int;
    // Assume mesh with homogenous element type, then eptr 
    // dictates an equal stride of size ElementType across EToV array.
    for (int i=0; i <= NumElements; i++) {
        eptr[i] = 3*i;
        cout << eptr[i] << endl;
    }

    cout << NumVerts << endl;
    for (int i=0; i < NumVerts; i++ ) {
        cout << eind[i] << endl;
    }
        
    *objval = 0;

    int * NE = new int;
    *NE = NumElements;
    int * NV = new int;
    *NV = 8; //NumVerts = 6

    int * ncommon = new int;
    *ncommon = 1;

    //cout << "About to call METIS_PartMeshNodal" << endl;
    //int result =  METIS_PartMeshNodal( NE, NV, eptr, eind, NULL, NULL,
    //                numPartitionsPtr, NULL, metisOptions, objval, epart, npart);

    cout << "About to call METIS_PartMeshDual" << endl;
    cout << metisOptions[METIS_OPTION_NUMBERING] << endl;
    int result = METIS_PartMeshDual(NE, NV, eptr, eind, NULL, NULL,
        ncommon, numPartitionsPtr, NULL, metisOptions, objval, epart, npart);

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
    for (int i=0; i<NumElements; i++)
        cout << epart[i] << endl;

    cout << "Vertex partitioning vector: " << endl;
    for (int i=0; i<NumVerts; i++)
        cout << npart[i] << endl;

    ElementPartitionMap = epart;
    VertexPartitionMap = npart;

}

/**
 * Read a list of vertices from a file. x-, y-, (and z-) are coordinates delimited by spaces, e.g., 0.5 1.0.
 */
void MeshManager::readVertices(string vertFile) {
    int * dims;
    readCsvFile<double>(vertFile, CsvDelimeters, Vert, dims);
    NumVerts = dims[0];
    Dim = dims[1];
}

/**
 * Read a list of elments from a file. Vertex numbers are written in a row and delimited by spaces, e.g., 1 2 3 4
 */
void MeshManager::readElements(string E2VFile) {
    int * dims;
    readCsvFile<int>(E2VFile, CsvDelimeters, EToV, dims);
    NumElements = dims[0];
    ElementType = dims[1];
}

/**
 * Print the list of vertices to stdout.
 */
void MeshManager::printVertices() {
    MeshManager::printArray<double>(Vert, NumVerts, Dim);
}

/**
 * Print the list of elements to stdout.
 */
void MeshManager::printElements() {
    MeshManager::printArray<int>(EToV, NumElements, ElementType);
}

/**
 * Returns a reference to the list of vertices, a contiguous block of type double.
 */
double * & MeshManager::get_Vertices() {
    return Vert;
}

/**
 * Returns the dimension of the vertex data. Usuallly will be 2 or 3.
 */
int MeshManager::get_Dim() {
    return Dim;
}

/**
 * Returns the number of vertices.
 */
int MeshManager::get_NumVerts() {
    return NumVerts;
}

/**
 * Returns the number of elements.
 */
int MeshManager::get_NumElements() {
    return NumElements;
}

/**
 * Returns the type of element. 3 => triangles, 4 => quadrilaterals, etc.
 */
int MeshManager::get_ElementType() {
    return ElementType;
}

/**
 * Returns a reference to the list of elements, a contiguous block of type int.
 */
int * & MeshManager::get_Elements() {
    return EToV;
}

/**
 * Returns a reference to the element partition map, an array of type int.
 */
int * & MeshManager::get_ElementPartitionMap() {
    return ElementPartitionMap;
}

/**
 * Returns a reference to the vertex partition map, an array of type int.
 */
int * & MeshManager::get_VertexPartitionMap() {
    return VertexPartitionMap;
}

MeshManager::~MeshManager() {
    if (Vert != nullptr) delete[] Vert;
    if (EToV != nullptr) delete[] EToV;
    if (ElementPartitionMap != nullptr) delete[] ElementPartitionMap;
    if (VertexPartitionMap != nullptr) delete[] VertexPartitionMap;
}
