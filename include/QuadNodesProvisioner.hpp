// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file QuadNodesProvisioner.hpp
 * @brief Defines the QuadNodesProvisioner class that provides a construction
 * of 2D nodal points on the Quad.
 */

#pragma once
#include "Nodes1DProvisioner.hpp"
#include "DGContext2D.hpp"
#include "MeshManager.hpp"
#include "Types.hpp"
#include "LinAlgHelpers.hpp"
#include <memory>

namespace blitzdg {
  /**
   * Provides facilities for the construction of two-dimensional
   * nodes, operators, and geometric factors on Quads.
   * @note This class is moveable but not copyable.
   */ 
  class QuadNodesProvisioner {
      using index_mat_smart_ptr = std::unique_ptr<index_matrix_type>;
      using index_vec_smart_ptr = std::unique_ptr<index_vector_type>;

      index_type NumElements;
      index_type NOrder;
      index_type NumLocalPoints;
      index_type NumFacePoints;

      real_mat_smart_ptr xGrid;
      real_mat_smart_ptr yGrid;

      real_vec_smart_ptr rGrid;
      real_vec_smart_ptr sGrid;

      real_mat_smart_ptr V;
      real_mat_smart_ptr Dr;
      real_mat_smart_ptr Ds;
      real_mat_smart_ptr Drw;
      real_mat_smart_ptr Dsw;
      real_mat_smart_ptr Lift;
      real_mat_smart_ptr J;
      real_mat_smart_ptr rx;
      real_mat_smart_ptr sx;
      real_mat_smart_ptr ry;
      real_mat_smart_ptr sy;
      real_mat_smart_ptr nx;
      real_mat_smart_ptr ny;
      real_mat_smart_ptr Vinv;
      real_mat_smart_ptr Filter;

      index_mat_smart_ptr Fmask;
      real_mat_smart_ptr Fx;
      real_mat_smart_ptr Fy;

      real_mat_smart_ptr Fscale;

      index_vec_smart_ptr vmapM;
      index_vec_smart_ptr vmapP;
      index_vec_smart_ptr vmapB;
      index_vec_smart_ptr mapP;
      index_vec_smart_ptr mapB;
      std::unique_ptr<index_hashmap> BCmap;

      const MeshManager& Mesh2D;

      std::unique_ptr<Nodes1DProvisioner> Nodes1D;
	  JacobiBuilders Jacobi;
      VandermondeBuilders Vandermonde;
      DirectSolver LinSolver;
      DenseMatrixInverter Inverter;

  public:
      static const index_type NumFaces;
      static const real_type NodeTol;
    /**
     * Constructor.
     * @param[in] NOrder Order of the approximating polynomials.
     * @param[in] NumElements The number of elements.
     * @param[in] MeshManager The MeshManager storing the 2D mesh.
     * @note Assumes uniform elements.
     */
    QuadNodesProvisioner(index_type _NOrder, const MeshManager& _MeshManager);

    /**
     * Copy constructor (deleted).
     */
    QuadNodesProvisioner(const QuadNodesProvisioner&) = delete;

    /**
     * Copy assignment operator (deleted).
     */
    QuadNodesProvisioner& operator=(const QuadNodesProvisioner&) = delete;

    /**
     * Move constructor.
     */
    QuadNodesProvisioner(QuadNodesProvisioner&&) = default;

    /**
     * Move assignment operator.
     */
    QuadNodesProvisioner& operator=(QuadNodesProvisioner&&) = default;

    /**
     * Builds nodes and local operators.
     */
    void buildNodes();

    void buildMaps();

    void buildBCHash();

    void buildBCHash(const index_vector_type& bcType);

    void buildFilter(real_type Nc, index_type s);


   /**
    * Computes the 2D Vandermonde matrix which maps modal coefficients to nodal values.
    * @param[in] N Order of the approximating polynomials.
    * @param[in] r 1st coordinate for standard rectangle.
    * @param[in] s 2nd coordinate for standard rectangle.
    * @param[out] V Vandermonde matrix.
    */
    void computeVandermondeMatrix(index_type N, const real_vector_type & r, const real_vector_type & s, real_matrix_type & V) const;


    void computeGradVandermondeMatrix(index_type N,  const real_vector_type & r, const real_vector_type & s, real_matrix_type & V2Dr, real_matrix_type & V2Ds) const;

    void computeDifferentiationMatrices(const real_matrix_type & V2Dr, const real_matrix_type & V2Ds, const real_matrix_type & V, real_matrix_type & Dr, real_matrix_type & Ds, real_matrix_type& Drw, real_matrix_type& Dsw) const;

    void computeInterpMatrix(const real_vector_type& rout, const real_vector_type& sout, real_matrix_type& IM) const;


    /**
     * Build x and y coordinates and geometric factors.
     */
    void buildPhysicalGrid();

    /**
     * Builds the Surface-To-Volume Lifting Operator.
     */
    void buildLift();

    DGContext2D get_DGContext() const;

    const index_vector_type& get_mapB() const { return *mapB; };
    const index_vector_type& get_vmapB() const { return * vmapB; };

    const real_matrix_type& get_Filter() const { return *Filter; };

    index_type get_NumLocalPoints() const { return NumLocalPoints; };

    const index_vector_type& get_vmapM() const { return *vmapM; };
    const index_vector_type& get_vmapP() const { return *vmapP; };

    const real_matrix_type & get_xGrid() const { return *xGrid; };
    const real_matrix_type & get_yGrid() const { return *yGrid; };

    const real_matrix_type & get_Lift() const { return *Lift; };

    const real_matrix_type & get_Dr() const { return *Dr; };

    const real_matrix_type & get_Ds() const { return *Ds; };

    const real_matrix_type & get_Drw() const { return *Drw; };

    const real_matrix_type & get_Dsw() const { return *Dsw; };

    const real_matrix_type & get_rx() const { return *rx; };

    const real_matrix_type & get_ry() const { return *ry; };

    const real_matrix_type & get_sx() const { return *sx; };

    const real_matrix_type & get_sy() const { return *sy; };

    const index_hashmap& get_bcMap() const { return *BCmap; };

    const real_matrix_type& get_Fscale() const { return *Fscale; };

    const index_matrix_type& get_Fmask() const { return *Fmask; };

    const real_matrix_type& get_nx() const { return *nx; };

    const real_matrix_type& get_ny() const { return *ny; };

    int get_NumElements() const { return NumElements; };

  };
}
