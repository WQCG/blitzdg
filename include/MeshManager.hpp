// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file MeshManager.hpp
 * @brief Defines the MeshManager class that reads mesh input files and partitions
 * meshes with metis. See http://glaros.dtc.umn.edu/gkhome/metis/metis/overview.
 */

#pragma once
#include "Types.hpp"
#include <string>

namespace blitzdg {
  class MeshManager {
      real_type* Vert;
      index_type* EToV;
      index_type Dim;
      index_type NumVerts;
      index_type ElementType;
      index_type NumElements;
      std::string CsvDelimeters;

      index_type* ElementPartitionMap;
      index_type* VertexPartitionMap;

      template<typename T>
      void  readCsvFile(std::string csvFile, std::string delimiters, T* & result, index_type* & dims);

      template<typename T>
      void printArray(T* & arr, index_type numRows, index_type numCols);

    public:
      MeshManager();

      index_type get_Index(index_type row, index_type col, index_type numCols);

      // Read gmsh .msh file.
      void readMesh(std::string gmshInputFile);
      
      void readVertices(std::string vertFile);

      void readElements(std::string E2VFile);

      // Split up with metis.
      void partitionMesh(index_type numPartitions);

      real_type* & get_Vertices();

      index_type get_Dim();

      index_type get_NumVerts();

      index_type* & get_Elements();

      index_type get_NumElements();
      index_type get_ElementType();

      index_type* & get_ElementPartitionMap();
      index_type* & get_VertexPartitionMap();

      void printVertices();
      void printElements();
      
      ~MeshManager();
  };
} // namespace blitzdg
