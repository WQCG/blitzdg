#include <iostream>
#include <math.h>
#include <MeshManager.hpp>
#include <gmsh/Gmsh.h>
#include <fstream>

using namespace std;

MeshManager::MeshManager() {

}


void MeshManager::ReadMesh(string gmshInputFile) {
    
}

void MeshManager::ReadMesh(string vertFile, string E2VFile) {
    string line;
    std::ifstream vertFileStream(vertFile);

    int idx=0;
    string token = "";
    while(getline(vertFileStream, line)) {
        //trim leading and trailing whitespace.
        line = line.erase(line.find_last_not_of(" \n\r\t")+1);
        line = line.substr(line.find_first_not_of(" \n\r\t"), line.length() );
        cout << line << endl;
//        while( (token =  line.substr(idx, line.find(' '))) != "" ) {
//            cout << token << endl;
//            idx = token.find(token);
//        }
    }
}



MeshManager::~MeshManager() {
}