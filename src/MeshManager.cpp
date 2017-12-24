#include <iostream>
#include <math.h>
#include <MeshManager.hpp>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <metis.h>

using namespace std;
using namespace boost;

MeshManager::MeshManager() {
    Vert = nullptr;
    CsvDelimeters = " \t";
    NumVerts = 0;
    NumElements = 0;
    Dim = 0;
}


int MeshManager::get_Index(int row, int col, int numCols) {
    return col + row*numCols;
}


void MeshManager::readMesh(string gmshInputFile) {
    throw("Not implemented!");
}

template<typename T>
vector<int> MeshManager::readCsvFile(string csvFile, string delimiters, T * & result) {

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

    vector<int> dims = {numLines,numCols}; 
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

    return dims;
}

void MeshManager::partitionMesh(int numPartitions) {
    int * metisOptions =  NULL;
    int objval =0 ;
    int * epart = NULL;
    int * npart = NULL;

    // set up mesh partitioning options
    metisOptions = new int[METIS_NOPTIONS];
    metisOptions[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // default
    metisOptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // total communication volume minimization.
    metisOptions[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
    metisOptions[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE;
    metisOptions[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY;
    metisOptions[METIS_OPTION_NCUTS] = 1;
    metisOptions[METIS_OPTION_NITER] = 10; // default
    metisOptions[METIS_OPTION_SEED] = 123; // random number seed.
    metisOptions[METIS_OPTION_UFACTOR] = 30; // max load imbalance of 1.03
    metisOptions[METIS_OPTION_NUMBERING] = 1; // 1-based numbering.
    metisOptions[METIS_OPTION_DBGLVL] = METIS_DBG_INFO; //debug level='info'. 0 for nothing.
    metisOptions[METIS_OPTION_CONTIG] = 1; // enforce a contiguous partition.

    int * eind = EToV;
    int * eptr = new int[NumElements+1];

    // output arrays
    epart = new int[NumElements];
    npart = new int[NumVerts];
    // Assume mesh with homogenous element type, then eptr 
    // dictates an equal stride of size ElementType across EToV array.
    for (int i=0; i <= NumElements; i++)
        eptr[i] = ElementType*i;

    int result =  METIS_PartMeshNodal( &NumElements, &NumVerts, eptr, eind, NULL, NULL,
                    &numPartitions, NULL, metisOptions, &objval, epart, npart);

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
    for (int i=0; i<NumElements; i++)
        cout << epart[i] << endl;

    cout << "Vertex partitioning vector: " << endl;
    for (int i=0; i<NumVerts; i++)
        cout << npart[i] << endl;

}

void MeshManager::readVertices(string vertFile) {
    cout << vertFile << endl;
    vector<int> dims = readCsvFile<double>(vertFile, CsvDelimeters, Vert);
    NumVerts = dims[0];
    Dim = dims[1];
}

void MeshManager::readElements(string E2VFile) {
    vector<int> dims = readCsvFile<int>(E2VFile, CsvDelimeters, EToV);
    NumElements = dims[0];
    ElementType = dims[1];
}

double * & MeshManager::get_Vertices() {
    return Vert;
}

int MeshManager::get_Dim() {
    return Dim;
}

int MeshManager::get_NumVerts() {
    return NumVerts;
}

int MeshManager::get_NumElements() {
    return NumElements;
}
int MeshManager::get_ElementType() {
    return ElementType;
}

int * & MeshManager::get_Elements() {
    return EToV;
}

MeshManager::~MeshManager() {
    if (Vert != nullptr) delete[] Vert;
}
