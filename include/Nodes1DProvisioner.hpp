// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Nodes1DProvisioner.hpp
 * @brief Defines a class for computations related to the construction of
 * one-dimensional nodes, operators, and geometric factors.
 * Based in part on the implementations described in "Nodal Discontinuous Galerkin Methods"
 * by J.S. Hesthaven and T. Warburton, 2008.
 */
#pragma once
#include "DirectSolver.hpp"
#include "JacobiBuilders.hpp"
#include "VandermondeBuilders.hpp"
#include "Types.hpp"
#include <memory>
#include <boost/python/numpy.hpp>

namespace blitzdg {
  /**
   * Provides facilities for the construction of one-dimensional
   * nodes, operators, and geometric factors.
   * @note This class is moveable but not copyable.
   */ 
  class Nodes1DProvisioner {
      using real_mat_smart_ptr = std::unique_ptr<real_matrix_type>;
      using index_mat_smart_ptr = std::unique_ptr<index_matrix_type>;
      using real_vec_smart_ptr = std::unique_ptr<real_vector_type>;
      using index_vec_smart_ptr = std::unique_ptr<index_vector_type>;

      real_type Min_x;
      real_type Max_x;
      index_type NumElements;
      index_type NOrder;
      index_type NumLocalPoints;

      index_type mapI;
      index_type mapO;
      index_type vmapI;
      index_type vmapO;
      
      real_mat_smart_ptr xGrid;
      real_vec_smart_ptr rGrid;

      real_mat_smart_ptr V;
      real_mat_smart_ptr Dr;
      real_mat_smart_ptr Lift;
      real_mat_smart_ptr J;
      real_mat_smart_ptr rx;
      real_mat_smart_ptr nx;
      real_mat_smart_ptr Vinv;

      index_vec_smart_ptr Fmask;
      real_mat_smart_ptr Fx;

      real_mat_smart_ptr Fscale;

      index_mat_smart_ptr EToV;
      index_mat_smart_ptr EToE;
      index_mat_smart_ptr EToF;

      index_vec_smart_ptr vmapM;
      index_vec_smart_ptr vmapP;

      DirectSolver LinSolver;
      JacobiBuilders Jacobi;
      VandermondeBuilders Vandermonde;

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
       * Copy constructor (deleted).
       */
      Nodes1DProvisioner(const Nodes1DProvisioner&) = delete;

      /**
       * Copy assignment operator (deleted).
       */
      Nodes1DProvisioner& operator=(const Nodes1DProvisioner&) = delete;

      /**
       * Move constructor.
       */
      Nodes1DProvisioner(Nodes1DProvisioner&&) = default;

      /**
       * Move assignment operator.
       */
      Nodes1DProvisioner& operator=(Nodes1DProvisioner&&) = default;
      
      /**
       * Builds the nodes and all geometric factors for each element.
       */
      void buildNodes();

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
       * Python bindings.
       */
      boost::python::numpy::ndarray get_xGrid_numpy() const;
      boost::python::numpy::ndarray get_Dr_numpy() const;
      boost::python::numpy::ndarray get_Fscale_numpy() const;
      boost::python::numpy::ndarray get_Lift_numpy() const;
      boost::python::numpy::ndarray get_rx_numpy() const;
      boost::python::numpy::ndarray get_vmapM_numpy() const;
      boost::python::numpy::ndarray get_vmapP_numpy() const;
      boost::python::numpy::ndarray get_nx_numpy() const;

  };
} // namespace blitzdg
