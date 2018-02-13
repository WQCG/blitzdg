#pragma once
#include <blitz/array.h>
#include <SparseMatrixConverter.hpp>
#include <EigenSolver.hpp>

using namespace std;
using namespace blitz;

class Nodes1DProvisioner {
    
    double Min_x;
    double Max_x;
    int NumElements;
    int NOrder;
    
    Array<double, 1> * xGrid;
    Array<double, 1> * rGrid;

    Array<double, 2> Dr;

    SparseMatrixConverter * MatrixConverter;
    EigenSolver * EigSolver;

  public:
    Nodes1DProvisioner(int NOrder, int NumElements, double xmin, double xmax, SparseMatrixConverter & converter, EigenSolver & eigenSolver);

    void buildNodes();

    void buildDr();

    Array<double, 1> & get_xGrid();
    Array<double, 1> & get_rGrid();
    Array<double, 2> & get_Dr();

    void computeJacobiPolynomial(Array<double,1> const & x, const double alpha, const double beta, const int N,  Array<double,1> & p);
    void computeJacobiQuadWeights(double alpha, double beta, int N, Array<double,1> & x, Array<double,1> & w);
    void computeGaussLobottoPoints(double alpha, double beta, int N, Array<double,1> & x);

    ~Nodes1DProvisioner();
};