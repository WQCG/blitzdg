// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file TriangleNodesProvisioner.hpp
 * @brief Defines the TriangleNodesProvisioner class that provides a construction
 * of 2D nodal points on the triangle.
 */

#pragma once
#include "NodesProvisioner2Dbase.hpp"
#include "Nodes1DProvisioner.hpp"
#include "DGContext2D.hpp"
#include "GaussFaceContext2D.hpp"
#include "TriangleCubatureRules.hpp"
#include "CubatureContext2D.hpp"
#include "CubatureContext2D.hpp"
#include "MeshManager.hpp"
#include "Types.hpp"
#include "LinAlgHelpers.hpp"
#include <memory>

namespace blitzdg {
  /**
   * Provides facilities for the construction of two-dimensional
   * nodes, operators, and geometric factors on triangles.
   * @note This class is moveable but not copyable.
   */ 
  class TriangleNodesProvisioner : public NodesProvisioner2DBase {
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
      real_mat_smart_ptr Fscale;
      real_mat_smart_ptr Fx;
      real_mat_smart_ptr Fy;

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
    TriangleNodesProvisioner(index_type _NOrder, const MeshManager& _MeshManager);

    /**
     * Copy constructor (deleted).
     */
    TriangleNodesProvisioner(const TriangleNodesProvisioner&) = delete;

    /**
     * Copy assignment operator (deleted).
     */
    TriangleNodesProvisioner& operator=(const TriangleNodesProvisioner&) = delete;

    /**
     * Move constructor.
     */
    TriangleNodesProvisioner(TriangleNodesProvisioner&&) = default;

    /**
     * Move assignment operator.
     */
    TriangleNodesProvisioner& operator=(TriangleNodesProvisioner&&) = default;


     /**
      * Compute (x,y) nodes in equilateral triangle for polynomial of order N.
      */
      void computeEquilateralNodes(real_vector_type & x, real_vector_type & y) const;

     /**
      * Computes the Jacobian (determinant), the geometric factor rx (dr/dx), 
      * and Fscale using nodes and differentiation matrix.
      */
      void computeJacobian();

     /**
      * Evaluates the 2D orthonormal polynomial on the simplex at (a,b) of order (i,j).
      */
      void evaluateSimplexPolynomial(const real_vector_type & a, const real_vector_type & b, index_type i, index_type j, real_vector_type & p) const;

     /**
      * Evaluates the gradient of the 2D orthonormal polynomial of index (id,jd) on the simplex at (a,b).
      */
      void evaluateGradSimplex(const real_vector_type & a, const real_vector_type & b, index_type id, index_type jd, real_vector_type & dpdr, real_vector_type & dpds) const;

     /**
      * Maps from (r,s) to (a,b) coordinates in the triangle.
      */     
      void rsToab(const real_vector_type & r, const real_vector_type & s, real_vector_type & a, real_vector_type & b) const;

     /**
      * Maps from (x,y) in the equilateral triangle to (r,s) coordinates in standard triangle
      */
      void xyTors(const real_vector_type & x, const real_vector_type & y, real_vector_type & r, real_vector_type & s) const;

     /**
      * Computes geometric warp function for interior nodes.
      */
      void computeWarpFactor(const real_vector_type & r, real_vector_type & warpFactor) const;

     /**
      * Computes the 2D Vandermonde matrix which maps modal coefficients to nodal values.
      * @param[in] N Order of the approximating polynomials.
      * @param[in] r 1st coordinate for standard triangle.
      * @param[in] s 2nd coordinate for standard triangle.
      * @param[out] V Vandermonde matrix.
      */
      void computeVandermondeMatrix(index_type N, const real_vector_type & r, const real_vector_type & s, real_matrix_type & V) const;

     /**
      * Computes the Gradient of the 2D Vandermonde matrix.
      * @param[in] N Order of the approximating polynomials;.
      * @param[in] r 1st coordinate for the standard triangle.
      * @param[in] s 2nd coordinate for the standard triangle.
      * @param[out] V2Dr r-component of the Gradient of the Vandermonde matrix.
      * @param[out] V2Ds s-component of the Gradient of the Vandermonde matrix.
      */
      void computeGradVandermondeMatrix(index_type N,  const real_vector_type & r, const real_vector_type & s, real_matrix_type & V2Dr, real_matrix_type & V2Ds) const;


     /**
      * Compute the two 2D differentiation matrices given the Gradient of the Vandermonde matrix and the Vandermonde matrix.
      * @param[in] V2Dr r-component of the Gradient of the Vandermonde matrix.
      * @param[in] V2Ds s-component of the Gradient of the Vandermonde matrix.
      * @param[in] V 2D Vandermonde matrix.
      * @param[out] Dr Differentiation matrix with respect to r coordinate.
      * @param[out] Ds Differentiation matrix with respect to s coordinate.
      */
      void computeDifferentiationMatrices(const real_matrix_type & V2Dr, const real_matrix_type & V2Ds, const real_matrix_type & V, real_matrix_type & Dr, real_matrix_type & Ds, real_matrix_type& Drw, real_matrix_type& Dsw) const;

     /**
      * Returns a reference to a matrix whose jth column contains the
      * x-component of the local grid points for the jth element.
      */
      const real_matrix_type & get_xGrid() const;

      /**
       * Returns a reference to a matrix whose jth column contains the
       * y-component of the local grid points for the jth element.
       */
      const real_matrix_type & get_yGrid() const;

      /**
       * Returns a reference to a vector that contains the 
       * r-component of the local grid points on the standard element.
       */
      const real_vector_type & get_rGrid() const;

      /**
       * Returns a reference to a vector that contains the 
       * s-component of the local grid points on the standard element.
       */
      const real_vector_type & get_sGrid() const;

      /**
       * Returns a reference to the differentiation matrix with respect
       * to the r component on the the standard element.
       */
      const real_matrix_type & get_Dr() const;

      /**
       * Returns a reference to the differentiation matrix with respect
       * to the s component on the the standard element.
       */
      const real_matrix_type & get_Ds() const;

      /**
       * Returns a reference to the weak differentiation matrix with respect
       * to the r component on the the standard element.
       */
      const real_matrix_type & get_Drw() const;

      /**
       * Returns a reference to the weak differentiation matrix with respect
       * to the s component on the the standard element.
       */
      const real_matrix_type & get_Dsw() const;

      /**
       * Returns a reference to the 2D generalized Vandermonde matrix
       * for 2D orthonormal basis functions on the triangle.
       */
      const real_matrix_type & get_V() const;

      /**
       * Returns a reference to the inverse of the 2D generalized Vandermonde matrix.
       */
      const real_matrix_type & get_Vinv() const;

      /**
       * Returns a reference to an exponential cutoff filter built with buildFilter().
       */
      const real_matrix_type & get_Filter() const;

      /**
       * Returns a reference to the Jacobian scaling matrix.
       */
      const real_matrix_type & get_J() const;

      /**
       * Returns a reference to the geometric scaling matrix rx.
       */
      const real_matrix_type & get_rx() const;

      /**
       * Returns a reference to the geometric scaling matrix ry.
       */
      const real_matrix_type & get_ry() const;

      /**
       * Returns a reference to the geometric scaling matrix sx.
       */
      const real_matrix_type & get_sx() const;

      /**
       * Returns a reference to the geometric scaling matrix sy.
       */
      const real_matrix_type & get_sy() const;

      /**
       * Returns a reference to a matrix whose jth column contains the
       * the x-coordinate of the unit normal vectors for each face point.
       */
      const real_matrix_type & get_nx() const;

      /**
       * Returns a reference to a matrix whose jth column contains the
       * the y-coordinate of the unit normal vectors for each face point.
       */
      const real_matrix_type & get_ny() const;

      /**
       * Returns a reference to the index-mask for the face nodes.
       */
      const index_vector_type & get_Fmask() const;

      /**
       * Returns a reference to the faces-only x-grid.
       */
      const real_matrix_type & get_Fx() const { return *Fx; }

      /**
       * Returns a reference to the faces-only y-grid.
       */
      const real_matrix_type & get_Fy() const { return *Fy; }

      /**
       * Returns a reference to the Face-scaling factor (inverse of Jacobian at face nodes).
       */
      const real_matrix_type & get_Fscale() const;

      /**
       * Returns a reference to the 2D lifting operator.
       */
      const real_matrix_type & get_Lift() const;

      /**
       * Returns a reference to the MeshManager.
       */
      const index_matrix_type & get_MeshManager() const;

      /**
       * Returns a refernce to the volume to surface map, 'minus' traces.
       */
      const index_vector_type & get_vmapM() const;

      /**
       * Returns a reference to the volume to surface map, 'plus' traces.
       */
      const index_vector_type & get_vmapP() const;

      /**
       * Returns a reference to the volume to surface map for boundary traces.
       */
      const index_vector_type & get_vmapB() const;

      /**
       * Returns a reference to the surface to surface map from 'minus' traces to boundary traces.
       */
      const index_vector_type & get_mapB() const;

      /**
       * Returns a reference to the hashmap mapping boundary types to list of corresponding boundary nodes.
       */
      const index_hashmap & get_bcMap() const;

      /**
       * Returns a new DG context, giving all fields necessary for implementing a solver.
       */
      DGContext2D get_DGContext() const;

      /**
       * Builds nodes and local operators.
       */
      void buildNodes();

      /**
       * Build x and y coordinates and geometric factors.
       */
      void buildPhysicalGrid();

      /**
       * Builds the Surface-To-Volume Lifting Operator.
       */
      void buildLift();

      /**
       * Build volume-to-surface maps for '-' and '+' traces, as well as for BC faces.
       */
      void buildMaps();

      /**
       *  Build hashmap mapping BCTypes to node indices that are on that boundary.
       */
      void buildBCHash();
      
      /**
       *  Build hashmap mapping BCTypes to node indices that are on that boundary using the
       *  input BCType table.
       */
      void buildBCHash(const index_vector_type& bcType);

      /**
       * Build exponential cut-off filter, where cut-off order is specified by Nc, and 
       * exponential power is determined by s.
       */
      void buildFilter(real_type Nc, index_type s);

      /**
       * Builds out the Gauss Face Nodes context.
       */
      GaussFaceContext2D buildGaussFaceNodes(index_type NGauss);

      // CubatureContext2D buildCubatureVolumeMesh(index_type NCubature);

      /**
       * Build an interpolation operator to a new set of nodes on the reference triangle.
       */
      void computeInterpMatrix(const real_vector_type& rout, const real_vector_type& sout, real_matrix_type& IM) const;

      /**
       * Splits high-order elements into smaller elements with linear basis functions.
       */
      void splitElements(const real_matrix_type& x, const real_matrix_type& y, const real_matrix_type& field, real_matrix_type& xnew, real_matrix_type& ynew, real_matrix_type& fieldnew) const;

      /**
       * Returns the number of nodes local to a 2D triangular element.
       */
      index_type get_NumLocalPoints() const;

      /**
       * Returns the number of nodes along a triangle face.
       */
      index_type get_NumFacePoints() const;

      /**
       * Returns the number of elements in the 2D grid.
       */
      index_type get_NumElements() const;

      /**
       * Returns order of the basis polynomials.
       */
      index_type get_NOrder() const { return NOrder; };
  };
} // namespace blitzdg

