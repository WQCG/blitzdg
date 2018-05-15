// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Nodes1DProvisioner.hpp
 * @brief Class for computations relating to construction of
 * one-dimensional nodes, operators, and geometric factors.
 * based in part on the implementations described in "Nodal Discontinuous Galerkin Methods"
 * by J.S. Hesthaven and T. Warburton, 2008.
 */

#pragma once
#include "SparseMatrixConverter.hpp"
#include "EigenSolver.hpp"
#include "DirectSolver.hpp"
#include "Types.hpp"

namespace blitzdg {
  class Nodes1DProvisioner {
      real_type Min_x;
      real_type Max_x;
      index_type NumElements;
      index_type NOrder;
      index_type NumLocalPoints;
      
      matrix_type* xGrid;
      vector_type* rGrid;

      matrix_type* V;
      matrix_type* Dr;
      matrix_type* Lift;
      matrix_type* J;
      matrix_type* rx;

      index_vector_type* Fmask;
      matrix_type* Fx;

      index_matrix_type* EToV;
      index_matrix_type* EToE;
      index_matrix_type* EToF;

      index_vector_type * vmapM;
      index_vector_type * vmapP;

      SparseMatrixConverter * MatrixConverter;
      EigenSolver * EigSolver;
      DirectSolver * LinSolver;

    public:
      static const index_type NumFacePoints;
      static const index_type NumFaces;
      static const real_type NodeTol;

      Nodes1DProvisioner(index_type NOrder, index_type NumElements, real_type xmin, real_type xmax, SparseMatrixConverter & converter, EigenSolver & eigenSolver, DirectSolver & directSolver);

      void buildNodes();
      void buildConnectivityMatrices();
      void buildFaceMask();
      void buildDr();
      void buildVandermondeMatrix();
      void buildLift();
      void buildMaps();
      void computeGradVandermonde(matrix_type& DVr);
      void computeJacobian();
      
      matrix_type& get_xGrid();
      vector_type& get_rGrid();
      matrix_type& get_Dr();
      matrix_type& get_V();
      matrix_type& get_J();
      matrix_type& get_rx();

      index_vector_type& get_Fmask();
      matrix_type& get_Fx();

      index_matrix_type& get_EToV();
      matrix_type& get_Lift();
    
      index_matrix_type& get_EToE();
      index_matrix_type& get_EToF();

      const index_vector_type & get_vmapM();
      const index_vector_type & get_vmapP();

      index_type get_NumLocalPoints();

      // these can be moved to a helper (polynomials) class or made private within this class.
      void computeJacobiPolynomial(vector_type const & x, const real_type alpha, const real_type beta, const index_type N, vector_type & p);
      void computeJacobiQuadWeights(real_type alpha, real_type beta, index_type N, vector_type & x, vector_type & w);
      void computeGaussLobottoPoints(real_type alpha, real_type beta, index_type N, vector_type & x);
      void computeGradJacobi(vector_type const & x, const real_type alpha, const real_type beta, const index_type N, vector_type & dp);

      ~Nodes1DProvisioner();
  };
} // namespace blitzdg