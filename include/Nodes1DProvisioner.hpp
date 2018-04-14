#pragma once
#include <blitz/array.h>
#include <SparseMatrixConverter.hpp>
#include <EigenSolver.hpp>
#include <DirectSolver.hpp>

using namespace std;
using namespace blitz;

class Nodes1DProvisioner {
    
    double Min_x;
    double Max_x;
    int NumElements;
    int NOrder;
    
    Array<double, 2> * xGrid;
    Array<double, 1> * rGrid;

    Array<double, 2> * V;
    Array<double, 2> * Dr;

    SparseMatrixConverter * MatrixConverter;
    EigenSolver * EigSolver;
    DirectSolver * LinSolver;

  public:
    Nodes1DProvisioner(int NOrder, int NumElements, double xmin, double xmax, SparseMatrixConverter & converter, EigenSolver & eigenSolver, DirectSolver & directSolver);

    void buildNodes();
    void buildDr();
    void buildVandermondeMatrix();
    void computeGradVandermonde(Array<double,2> & DVr);
    void computeJacobian(Array<double,2> & J, Array<double,2> & rx);

    Array<double, 2> & get_xGrid();
    Array<double, 1> & get_rGrid();
    Array<double, 2> & get_Dr();
    Array<double, 2> & get_V();

    // these can be moved to a helper (polynomials) class or made private within this class.
    void computeJacobiPolynomial(Array<double,1> const & x, const double alpha, const double beta, const int N, Array<double,1> & p);
    void computeJacobiQuadWeights(double alpha, double beta, int N, Array<double,1> & x, Array<double,1> & w);
    void computeGaussLobottoPoints(double alpha, double beta, int N, Array<double,1> & x);
    void computeGradJacobi(Array<double,1> const & x, const double alpha, const double beta, const int N, Array<double,1> & dp);

    ~Nodes1DProvisioner();
};