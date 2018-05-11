// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#pragma once
#include <blitz/array.h>
#include <SparseMatrixConverter.hpp>
#include <EigenSolver.hpp>
#include <DirectSolver.hpp>
#include <Types.hpp>

using namespace std;
using namespace blitz;

class Nodes1DProvisioner {
    double Min_x;
    double Max_x;
    int NumElements;
    int NOrder;
    int NumLocalPoints;

    blitzdg::index_type mapI;
    blitzdg::index_type mapO;
    blitzdg::index_type vmapI;
    blitzdg::index_type vmapO;
    
    Array<double, 2> * xGrid;
    Array<double, 1> * rGrid;

    Array<double, 2> * V;
    Array<double, 2> * Dr;
    Array<double, 2> * Lift;
    Array<double, 2> * J;
    Array<double, 2> * rx;

    Array<int, 1> * Fmask;
    Array<double, 2> * Fx;

    Array<int, 2> * EToV;
    Array<int, 2> * EToE;
    Array<int, 2> * EToF;

    blitzdg::index_vector_type * vmapM;
    blitzdg::index_vector_type * vmapP;

    SparseMatrixConverter * MatrixConverter;
    EigenSolver * EigSolver;
    DirectSolver * LinSolver;

  public:
    static const int NumFacePoints;
    static const int NumFaces;
    static const double NodeTol;

    Nodes1DProvisioner(int NOrder, int NumElements, double xmin, double xmax, SparseMatrixConverter & converter, EigenSolver & eigenSolver, DirectSolver & directSolver);

    void buildNodes();
    void buildConnectivityMatrices();
    void buildFaceMask();
    void buildDr();
    void buildVandermondeMatrix();
    void buildLift();
    void buildMaps();
    void computeGradVandermonde(Array<double,2> & DVr);
    void computeJacobian();
    
    Array<double, 2> & get_xGrid();
    Array<double, 1> & get_rGrid();
    Array<double, 2> & get_Dr();
    Array<double, 2> & get_V();
    Array<double, 2> & get_J();
    Array<double, 2> & get_rx();

    Array<int, 1> & get_Fmask();
    Array<double, 2> & get_Fx();

    Array<int, 2> & get_EToV();
    Array<double, 2> & get_Lift();
  
    Array<int, 2> & get_EToE();
    Array<int, 2> & get_EToF();

    const blitzdg::index_vector_type & get_vmapM();
    const blitzdg::index_vector_type & get_vmapP();

    const blitzdg::index_type get_mapI();
    const blitzdg::index_type get_mapO();
    const blitzdg::index_type get_vmapI();
    const blitzdg::index_type get_vmapO();

    int get_NumLocalPoints();

    // these can be moved to a helper (polynomials) class or made private within this class.
    void computeJacobiPolynomial(Array<double,1> const & x, const double alpha, const double beta, const int N, Array<double,1> & p);
    void computeJacobiQuadWeights(double alpha, double beta, int N, Array<double,1> & x, Array<double,1> & w);
    void computeGaussLobottoPoints(double alpha, double beta, int N, Array<double,1> & x);
    void computeGradJacobi(Array<double,1> const & x, const double alpha, const double beta, const int N, Array<double,1> & dp);

    ~Nodes1DProvisioner();
};