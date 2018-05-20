// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file MeshManager.hpp
 * @brief Defines the MeshManager class that reads mesh input files and partitions
 * meshes with METIS. See http://glaros.dtc.umn.edu/gkhome/metis/metis/overview.
 */

#pragma once
#include "Types.hpp"
#include <string>

namespace blitzdg {
  /**
   * Reads mesh input files and partitions meshes with METIS.
   */
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
      /**
       * Default constructor.
       */
      MeshManager();

      /**
      * Converts a (row, col) index to an integer index into a contiguous block of memory.
      * @param[in] row Row index.
      * @param[in] col Column index.
      * @param[in] numCols The number of columns.
      * @return The linear index.
      */
      static index_type get_Index(index_type row, index_type col, index_type numCols) {
          return col + row*numCols;
      }

      /**
       * Reads a .msh file that was generated by Gmsh.
       * @param[in] gmshInputFile Full path to .msh file.
       */
      void readMesh(const std::string& gmshInputFile);
      
      /**
       * Reads a list of vertices from a file. 
       * 
       * The x-, y-, and optiontally z-coordinates are stored columnwise and
       * are delimited by spaces within a row. For example, the following
       * represents a valid file:
       * @verbatim
       * 0.5 1.0
       * 0.2 3.0
       * 0.1 1.0
       * @endverbatim
       * @note All rows in the file must contain the same number of coordinates, 
       * either two or three.
       * @param[in] vertFile Full path to the vertex file.
       */
      void readVertices(const std::string& vertFile);

      /**
       * Reads a list of elments from a file. 
       * 
       * An element is defined by its vertex indices. Each row of the file
       * defines an element by listing its vertex indices delimited by spaces.
       * For example, the following represents a valid file:
       * @verbatim
       * 1 2 3 4
       * 4 5 6 7
       * 8 7 9 5
       * @endverbatim
       * @note The vertex indices of an element should index the vertices read
       * from the vertex file.
       * @see readVertices()
       * @param[in] E2VFile Full path to the element index file.
       */
      void readElements(const std::string& E2VFile);

      /**
       * Partitions a mesh using METIS. Generates the vertex and
       * element partition maps.
       * 
       * Results can be obtained by calling MeshManager.get_ElementPartitionMap()
       * and MeshManager.get_VertexPartitionMap().
       * @param[in] numPartitions The number of partitions.
       */
      void partitionMesh(index_type numPartitions);

      /**
       * Returns a pointer to the array of vertex coordinates.
       * 
       * If vd is the integer denoting the vertex dimension and
       * nv is the number of vertices, then the jth coordinate
       * of the ith vertex may be accessed by
       * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       * const index_type* coords = mm.get_Vertices();
       * coords(MeshManager::get_Index(i, j, vd));
       * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       * for j = 0,...,vd-1 and i = 0,...,nv-1.
       */
      const real_type* get_Vertices() const;

      /**
       * Returns the dimension of a vertex, typically 2 or 3.
       */
      index_type get_Dim() const;

      /**
       * Returns the number of vertices.
       */
      index_type get_NumVerts() const;

      /**
       * Returns a pointer to the array of element indices.
       * 
       * If et is the integer denoting the element type and
       * ne is the number of elements, then the jth index
       * of the ith element may be accessed by
       * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       * const index_type* elems = mm.get_Elements();
       * elems(MeshManager::get_Index(i, j, et));
       * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       * for j = 0,...,et-1 and i = 0,...,ne-1.
       */
      const index_type* get_Elements() const;

      /**
       * Returns the number of elements.
       */
      index_type get_NumElements() const;
     
      /**
       * Returns an integer that defines the element type.
       * Returns 3 for triangles, 4 for a quadrilaterals, etc.
       */
      index_type get_ElementType() const;

      /**
       * Returns a pointer to the underlying array that stores
       * the element partition map.
       * 
       * Has length equal to the number of elements. It contains 
       * the processor rank of each element, which is used 
       * for distributed computing with MPI.
       */
      const index_type* get_ElementPartitionMap() const;
      
      /**
       * Returns a pointer to the underlying array that stores
       * the vertex partition map.
       * 
       * Has length equal to the number of vertices. It contains 
       * the processor rank of each vertex, which is used 
       * for distributed computing with MPI.
       */
      const index_type* get_VertexPartitionMap() const;

      /**
       * Write the vertices to the console.
       */
      void printVertices() const;

      /**
       * Write the elements to the console.
       */
      void printElements() const;
      
      /**
       * Destructor.
       */
      ~MeshManager();
  };
} // namespace blitzdg
