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
      std::string CsvDelimiters;
      index_type* ElementPartitionMap;
      index_type* VertexPartitionMap;

    public:
      MeshManager();

      static index_type get_Index(index_type row, index_type col, index_type numCols);

      // Read gmsh .msh file.
      void readMesh(const std::string& gmshInputFile);
      
      void readVertices(const std::string& vertFile);

      void readElements(const std::string& E2VFile);

      // Split up with metis.
      void partitionMesh(index_type numPartitions);

      const real_type* get_Vertices() const;

      index_type get_Dim() const;

      index_type get_NumVerts() const;

      const index_type* get_Elements() const;

      index_type get_NumElements() const;
      index_type get_ElementType() const;

      const index_type* get_ElementPartitionMap() const;
      const index_type* get_VertexPartitionMap() const;

      void printVertices() const;
      void printElements() const;
      
      ~MeshManager();
  };


} // namespace blitzdg
