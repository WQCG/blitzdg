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
}


int MeshManager::GetIndex(int row, int col, int numCols) {
    return col + row*numCols;
}


void MeshManager::ReadMesh(string gmshInputFile) {
    
}

void MeshManager::ReadVertices(string vertFile) {
    string str1("hello world");
    to_upper(str1);
    
    ifstream vertFileStream(vertFile);

    string line("");

    vector<string> splitVec;
    int numLines=0;
    while(getline(vertFileStream, line)) {
        // Take first line as source of truth for number of columns.
        if (numLines == 0) {
            trim(line);
            split( splitVec, line, is_any_of(" \t"), token_compress_on );
            Dim = splitVec.size();
        }
        numLines++;
    }
    NumVerts = numLines;
    
    // roll-back stream.
    vertFileStream.clear();
    vertFileStream.seekg(0, std::ios::beg);

    // Allocate vertex storage.
    Vert = new double[numLines*Dim];

    int count = 0;
    while(getline(vertFileStream, line)) {
        trim(line);
        vector<string> splitVec;
        split( splitVec, line, is_any_of(" \t"), token_compress_on ); 
        
        for(int i=0; i < Dim; i++) {
            Vert[count] = atof(splitVec[i].c_str());
            count++;
        }
    }
    vertFileStream.close();
}

double * & MeshManager::GetVertices() {
    return Vert;
}

int MeshManager::GetDim() {
    return Dim;
}

int MeshManager::GetNumVerts() {
    return NumVerts;
}

MeshManager::~MeshManager() {
    if (Vert != nullptr) delete[] Vert;
}
