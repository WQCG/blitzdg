// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Nodes1DProvisioner.hpp
 * @brief Class for computations relating to construction of
 * one-dimensional nodes, operators, and geometric factors.
 * Based in part on the implementations described in "Nodal Discontinuous Galerkin Methods"
 * by J.S. Hesthaven and T. Warburton, 2008.
 */

#pragma once
#include "DirectSolver.hpp"
#include "EigenSolver.hpp"
#include "JacobiBuilders.hpp"
#include "DenseMatrixInverter.hpp"
#include "Types.hpp"

namespace blitzdg {
  /**
   * Provides facilities for the construction of one-dimensional
   * nodes, operators, and geometric factors.
   * 
   * Facilities provided include:
   * <ul>
   * <li> Construction of one-dimensional nodes over the computational domain based
   * on Gauss-Jacobi-Lobotto quadrature rules. </li>
   * 
   * <li> </li>
   * 
   * <li> </li>
   * 
   * <li> </li>
   * </ul>
   */ 
  class Nodes1DProvisioner {
      real_type Min_x;
      real_type Max_x;
      index_type NumElements;
      index_type NOrder;
      index_type NumLocalPoints;

      index_type mapI;
      index_type mapO;
      index_type vmapI;
      index_type vmapO;
      
      real_matrix_type* xGrid;
      real_vector_type* rGrid;

      real_matrix_type* V;
      real_matrix_type* Dr;
      real_matrix_type* Lift;
      real_matrix_type* J;
      real_matrix_type* rx;
      real_matrix_type* nx;

      real_matrix_type Vinv;

      index_vector_type* Fmask;
      real_matrix_type* Fx;

      real_matrix_type* Fscale;

      index_matrix_type* EToV;
      index_matrix_type* EToE;
      index_matrix_type* EToF;

      index_vector_type* vmapM;
      index_vector_type* vmapP;

      EigenSolver EigSolver;
      DirectSolver LinSolver;
      JacobiBuilders Jacobi;
      DenseMatrixInverter Inverter;

    public:
      /**
       * The number of nodes per face of a 1D element.
       * 
       * In one dimension a face is just a single point, so this
       * variable is always set to 1.
       */
      static const index_type NumFacePoints;
      /**
       * The number of faces for a 1D element.
       * 
       * In one dimension an element is just a line segment
       * with two endpoints, so this variable is always set to 2.
       */
      static const index_type NumFaces;
      /**
       * Tolerance that is used to determine when two nodes are
       * sufficiently close to be considered equal.
       * 
       * Two nodes with coordinates \f$x\f$ and \f$y\f$ are considered
       * equal if \f$|x - y| < \tau\f$, where \f$\tau\f$ is the node
       * tolerance, in this case always set to 1e-5.
       */
      static const real_type NodeTol; 

      /**
       * Constructor.
       * @param[in] NOrder Order of the Gauss-Lobotto quadrature rule.
       * @param[in] NumElements The number of elements.
       * @param[in] xmin Left endpoint of the computational domain.
       * @param[in] xmax Right endpoint of the computational domain.
       * @note Assumes uniform elements.
       */
      Nodes1DProvisioner(index_type NOrder, index_type NumElements, real_type xmin, real_type xmax);

      /**
       * Builds the nodes and all geometric factors for each element.
       */
      void buildNodes();

      /**
       * Builds the global connectivity matrices (EToE, EToF) for 1D grid
       * based on the EToV (Element-to-Vertex) matrix.
       */
      void buildConnectivityMatrices();

      /**
       *  Builds the face mask that when applied to volume nodes gives the surface nodes.
       */
      void buildFaceMask();

      /**
       * Builds the differentiation matrix Dr on the standard element.
       */
      void buildDr();

      /**
       * Computes the Vandermonde matrix which maps modal coefficients to nodal values.
       */
      void buildVandermondeMatrix();

      /**
       * Builds the lifting operator.
       */
      void buildLift();

      /**
       * Builds the volume to surface maps.
       */
      void buildMaps();

      /**
       * Builds the unit normal vectors to the faces of the
       * 1D elements.
       */
      void buildNormals();

      /**
       * Builds the elementwise derivative of the Vandermonde matrix.
       * @param[out] DVr Elementwise derivative of Vandermonde matrix.
       */
      void computeGradVandermonde(real_matrix_type& DVr) const;

     /**
      * Computes the Jacobian (determinant), the geometric factor rx (dr/dx), 
      * and Fscale using nodes and differentiation matrix.
      */
      void computeJacobian();
      
      /**
       * Returns a reference to a matrix whose jth column contains the
       * local grid points for the jth element.
       */
      const real_matrix_type & get_xGrid() const;

      /**
       * Returns a reference to a vector that contains the 
       * Gauss-Jacobi-Lobotto points on the standard element,
       * i.e., on the interval [-1,1].
       */
      const real_vector_type & get_rGrid() const;

      /**
       * Returns a reference to the differentiation matrix on the
       * the standard element, i.e., the interval [-1,1].
       */
      const real_matrix_type & get_Dr() const;

      /**
       * Returns a reference to the generalized Vandermonde matrix.
       */
      const real_matrix_type & get_V() const;

      /**
       * Returns a reference to the inverse of the generalized Vandermonde matrix.
       */
      const real_matrix_type & get_Vinv() const;

      /**
       * Returns a reference to the Jacobian scaling matrix.
       */
      const real_matrix_type & get_J() const;

      /**
       * Returns a reference to the geometric scaling matrix.
       */
      const real_matrix_type & get_rx() const;

      /**
       * Returns a reference to a matrix whose jth column contains the
       * the x-coordinate of the unit normal vectors for each face point.
       */
      const real_matrix_type & get_nx() const;

      /**
       * Returns a reference to the index-mask for the face nodes.
       */
      const index_vector_type & get_Fmask() const;

      /**
       * Returns a reference to the faces only x-grid.
       */
      const real_matrix_type & get_Fx() const;

      /**
       * Returns a reference to the Face-scaling factor (inverse of Jacobian at face nodes).
       */
      const real_matrix_type & get_Fscale() const;

      /**
       * Returns a reference to the 1D lifting operator.
       */
      const real_matrix_type & get_Lift() const;

      /**
       * Returns a reference to the Element-To-Vertex connectivity table.
       */
      const index_matrix_type & get_EToV() const;
    
      /**
       * Returns a reference to the Element-To-Element connectivity table.
       */
      const index_matrix_type & get_EToE() const;

      /**
       * Returns a reference to the Element-to-Face connectivity table.
       */
      const index_matrix_type & get_EToF() const;

      /**
       * Returns a refernce to the volume to surface map, 'minus' traces.
       */
      const index_vector_type & get_vmapM() const;

      /**
       * Returns a reference to the volume to surface map, 'plus' traces.
       */
      const index_vector_type & get_vmapP() const;

      /**
       * Returns the surface index of the inflow boundary.
       */
      index_type get_mapI() const;

      /**
       * Returns the surface index of the outflow boundary.
       */
      index_type get_mapO() const;

      /**
       * Returns the volume index of the inflow boundary.
       */
      index_type get_vmapI() const;

      /**
       * Returns the volume index of the outflow boundary.
       */
      index_type get_vmapO() const;
      
      /**
       * Returns the number of nodes local to a 1D element.
       */
      index_type get_NumLocalPoints() const;

      /**
       * Returns the number of elements in the 1D grid.
       */
      index_type get_NumElements() const;
      
      /**
       * Destructor.
       */
      ~Nodes1DProvisioner();
  };
} // namespace blitzdg
