// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include <blitz/array.h>
#include <boost/algorithm/string.hpp>
#include <vector>

using namespace std;
using namespace boost;

class MeshManager {
    double * Vert;
    int * EToV;
    int Dim;
    int NumVerts;
    int ElementType;
    int NumElements;
    string CsvDelimeters;

    int * ElementPartitionMap;
    int * VertexPartitionMap;

    template<typename T>
    vector<int> readCsvFile(string csvFile, string delimiters, T * & result);

    template<typename T>
    void printArray(T * & arr, int numRows, int numCols);

  public:
    MeshManager();

    int get_Index(int row, int col, int numCols);

    // Read gmsh .msh file.
    void readMesh(string gmshInputFile);
    
    void readVertices(string vertFile);

    void readElements(string E2VFile);

    // Split up with metis.
    void partitionMesh(int numPartitions);

    double * & get_Vertices();

    int get_Dim();

    int get_NumVerts();

    int * & get_Elements();

    int get_NumElements();
	  int get_ElementType();

    int * & get_ElementPartitionMap();
    int * & get_VertexPartitionMap();

    void printVertices();
    void printElements();
    
    ~MeshManager();
};

