#include <iostream>
#include <math.h>
#include <Nodes1DProvisioner.hpp>

using namespace std;

/**
 * Constructor. Takes order of polynomials, number of elements, and dimensions of the domain.
 * Assumes equally-spaced elements.
 */
Nodes1DProvisioner::Nodes1DProvisioner(int _NOrder, int _NumElements, double _xmin, double _xmax, SparseMatrixConverter & converter, EigenSolver & eigenSolver) {
    NOrder = _NOrder;
    NumElements = _NumElements;
    Min_x = _xmin;
    Max_x = _xmax;
    MatrixConverter = &converter;
    EigSolver = &eigenSolver;

    rGrid = new Array<double,1>(NOrder+1);
}

/**
 * Build nodes and geometric factors for all elements.
 */
void Nodes1DProvisioner::buildNodes() {
    const double alpha = 0.0;
    const double beta = 0.0;

    computeGaussLobottoPoints(alpha, beta, NOrder, *rGrid);
}

/**
 * Compute Vandermonde matrix which maps modal coefficients to nodal values.
 */
void Nodes1DProvisioner::buildVandermondeMatrix() {
    V = new Array<double, 2>(NOrder+1, NOrder+1);

    Array<double, 2> & Vref = *V;

    Array<double, 1> p(NOrder+1);

    for (int j=1; j <= NOrder+1; j++) {
        computeJacobiPolynomial(get_rGrid(), 0.0, 0.0, j-1, p);
        Vref(Range::all(), j) = p;
    }
}

/**
 * Build differentiation matrix Dr on the standard element.
 */
void Nodes1DProvisioner::buildDr() {
    throw("Not implemented.");
} 

/**
 * Get reference to physical x-grid.
 */
Array<double, 1> & Nodes1DProvisioner::get_xGrid() {
    return *xGrid;
}

/**
 * Get reference to r-grid on the standard element.
 */
Array<double, 1> & Nodes1DProvisioner::get_rGrid() {
    return *rGrid;
}


/**
 * Get reference to differentiation matrix Dr on the standard element.
 */
Array<double, 2> & Nodes1DProvisioner::get_Dr() {
    throw("Not implemented.");
}

/**
 * Get reference to generalized Vandermonde matrix V.
 */
Array<double, 2> & Nodes1DProvisioner::get_V() {
    return *V;
}

/**
 * Destructor
 */
Nodes1DProvisioner::~Nodes1DProvisioner() {
    if (rGrid == nullptr) { delete[] rGrid; }
}

/**  Compute the Nth Jacobi polynomial of type (alpha,beta) > -1 ( != -0.5)
  *   and weights, w, associated with the Jacobi polynomial at the points x.
  */
void Nodes1DProvisioner::computeJacobiPolynomial(Array<double,1> const & x, const double alpha, const double beta, const int N, Array<double,1> & p) {
    Range all = Range::all();
    int Np = (x.length())(0);

    Array<double,2> pStorage(N+1, Np);

    double gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);

    p = 1/sqrt(gamma0);
    pStorage(0, all) = p;

    if (N==0) 
        return;

    double gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
    p = ((alpha+beta+2)*x/2 + (alpha-beta)/2)/sqrt(gamma1);

    pStorage(1, all) = p;

    if (N==1) 
        return;

    double aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

    // Forward recurrence using the symmetry of the recurrence.
    for(int i=1; i <= N-1; i++) {
        double h1 = 2*i+alpha+beta;
        double anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3));
        double bnew = - (alpha*alpha-beta*beta)/h1/(h1+2);
        pStorage(i+1,all) = 1/anew*( -aold*pStorage(i-1,all) + (x-bnew)*pStorage(i,all));
        aold = anew;
    }
    p = pStorage(N, all);
}

/**  Compute the Nth order Gauss quadrature points, x,
  *   and weights, w, associated with the Jacobi polynomial, of type (alpha,beta) > -1 ( != -0.5).
  */
void Nodes1DProvisioner::computeJacobiQuadWeights(double alpha, double beta, int N, Array<double,1> & x, Array<double,1> & w) {

    if ( N == 0) {
        x(0) = -(alpha-beta)/(alpha+beta+2);
        w(0) = 2.0;
        return;
    }

    const double eps = numeric_limits<double>::epsilon();

    firstIndex ii;
    secondIndex jj;

    // Form symmetric matrix.
    Array<double, 2> J(N+1,N+1);
    for (int i=0; i < N+1; i++) {
        double h1 = 2.*i+alpha+beta;
        J(i,i)   = -0.5*(alpha*alpha-beta*beta)/(h1+2.)/h1;
        if (i < N) {
            J(i,i+1) = 2./(h1+2)*sqrt((i+1)*((i+1)+alpha+beta)*((i+1)+alpha)*((i+1)+beta)/(h1+1)/(h1+3));
        }
        if ( i > 0) {
            J(i, i-1) = 0.0;
        }
    }

    if ((alpha + beta) < 10*eps ) {
        J(0,0) = 0.0;
    }
    J = J(ii,jj) + J(jj,ii);

    //SparseMatrixConverter & matConverter = *MatrixConverter;
    EigenSolver & eigSolver = *EigSolver;
    
    Array<double, 2> eigenvectors(N+1, N+1);

    eigSolver.solve(J, x, eigenvectors);

    // The eigenvalues give the x points.
    
    // The weights are given by:
    Array<double, 1> v1(N+1);
    v1 = eigenvectors( 0, Range::all() ); 

    double gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);

    w = (v1*v1)*gamma0;
}

/**  Compute the Nth order Gauss Lobatto quadrature points, x,
  *  associated with the Jacobi polynomial, of type (alpha,beta) > -1 ( != -0.5).
  */
void Nodes1DProvisioner::computeGaussLobottoPoints(double alpha, double beta, int N, Array<double,1> & x) {
    if (N==1) {
        x(0) = -1.0;
        x(1) = 1.0;
        return;
    }

    x(0) = -1.0;
    x(N) = 1.0;

    Array<double, 1> xJG(N-1);
    Array<double, 1> w(N-1);

    computeJacobiQuadWeights(alpha+1., beta+1., N-2, xJG, w);
    
    for(int i=1; i < N; i++)
        x(i) = xJG(i-1);
}                