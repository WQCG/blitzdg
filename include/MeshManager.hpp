// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include <blitz/array.h>
#include <string>

namespace blitzdg {
  class MeshManager {
      double * Vert;
      int * EToV;
      int Dim;
      int NumVerts;
      int ElementType;
      int NumElements;
      std::string CsvDelimeters;

      int * ElementPartitionMap;
      int * VertexPartitionMap;

      template<typename T>
      void  readCsvFile(std::string csvFile, std::string delimiters, T * & result, int * & dims);

      template<typename T>
      void printArray(T * & arr, int numRows, int numCols);

    public:
      MeshManager();

      int get_Index(int row, int col, int numCols);

      // Read gmsh .msh file.
      void readMesh(std::string gmshInputFile);
      
      void readVertices(std::string vertFile);

      void readElements(std::string E2VFile);

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
} // namespace blitzdg
