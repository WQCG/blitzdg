#include <iostream>
#include <math.h>
#include <MeshManager.hpp>
#include <gmsh/Gmsh.h>
#include <fstream>
#include <boost/algorithm/string.hpp>

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
