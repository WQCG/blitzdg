// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include <iostream>
#include <math.h>
#include <Nodes1DProvisioner.hpp>
#include <SparseTriplet.hpp>
#include <Types.hpp>

using namespace std;
using namespace blitz;
using namespace blitzdg;

const int Nodes1DProvisioner::NumFacePoints = 1;
const int Nodes1DProvisioner::NumFaces = 2;
const double Nodes1DProvisioner::NodeTol = 1.e-5;

/**
 * Constructor. Takes order of polynomials, number of elements, and dimensions of the domain.
 * Assumes equally-spaced elements.
 */
Nodes1DProvisioner::Nodes1DProvisioner(int _NOrder, int _NumElements, double _xmin, double _xmax, SparseMatrixConverter & converter, EigenSolver & eigenSolver, DirectSolver & directSolver) {
    NOrder = _NOrder;
    NumElements = _NumElements;
    Min_x = _xmin;
    Max_x = _xmax;
    MatrixConverter = &converter;
    EigSolver = &eigenSolver;
    LinSolver = &directSolver;

    // This is true in 1D only.
    NumLocalPoints = NOrder + 1;
    mapI = 0;
    mapO = NumFacePoints*NumFaces*NumElements - 1;
    vmapI = 0;
    vmapO = NumLocalPoints*NumElements - 1;

    rGrid = new Array<double, 1>(NumLocalPoints);
    xGrid = new Array<double, 2>(NumLocalPoints, NumElements);
    J =  new Array<double, 2>(NumLocalPoints, NumElements);
    rx = new Array<double, 2>(NumLocalPoints, NumElements);

    Lift = new Array<double, 2>(NumLocalPoints, NumFacePoints*NumFaces);
    EToV = new Array<int, 2>(NumElements, NumFaces);
    EToE = new Array<int, 2>(NumElements, NumFaces);
    EToF = new Array<int, 2>(NumElements, NumFaces);
    Fmask = new Array<int, 1> (NumFacePoints*NumFaces);
    Fx = new Array<double, 2>(NumFacePoints*NumFaces, NumElements);
    Fscale = new matrix_type(NumFacePoints*NumFaces, NumElements);
    nx = new matrix_type(NumFacePoints*NumFaces, NumElements);

    vmapM = new index_vector_type(NumFacePoints*NumFaces*NumElements);
    vmapP = new index_vector_type(NumFacePoints*NumFaces*NumElements);
}

/**
 * Build nodes and geometric factors for all elements.
 */
void Nodes1DProvisioner::buildNodes() {
    const double alpha = 0.0;
    const double beta = 0.0;

    Array<double,1> & r = *rGrid;

    computeGaussLobottoPoints(alpha, beta, NOrder, r);
    
    buildVandermondeMatrix();
    buildDr();
    buildLift();

    double L = Max_x - Min_x;
    double width = L / NumElements;

    Array<double, 2> & x = *xGrid;
    for (int k=0; k < NumElements; k++) {
        x(Range::all(), k) = Min_x + width*(k + 0.5*(r+1.));
    }

    Array<int, 2> & E2V = *EToV;

    // Create Element-to-Vertex connectivity table.
    for (int k=0; k < NumElements; k++) {
        E2V(k, 0) = k;
        E2V(k, 1) = k+1;
    }

    buildConnectivityMatrices();
    buildFaceMask();
	buildMaps();
    buildNormals();
}

/**
 * Build unit normals at element faces. Trivial in 1D.
 */
void Nodes1DProvisioner::buildNormals() {
    matrix_type & nxref = *nx;

    real_type mult = -1;

    for (index_type k=0; k < NumElements; k++) {
        for (index_type f=0; f < NumFaces*NumFacePoints; f++) {
            nxref(f,k) = mult;
            mult *= -1;
        }
    }
}

/**
 * Build volume to surface maps.
 */
void Nodes1DProvisioner::buildMaps() {
    firstIndex ii;
    secondIndex jj;

    index_matrix_type nodeIds(NumLocalPoints, NumElements);

    // Set up reference to the objects we need to interact with.
    SparseMatrixConverter & matConverter = *MatrixConverter;
    
    index_vector_type & Fmsk = *Fmask;
    index_vector_type & vmM = *vmapM;
    index_vector_type & vmP = *vmapP;

    index_matrix_type & E2E = *EToE;
    index_matrix_type & E2F = *EToF;
    matrix_type & xmat = *xGrid;

    matrix_type xmatTrans(NumElements, NumLocalPoints);
    xmatTrans = xmat(jj,ii);

    double * x = new double[NumElements*NumLocalPoints];
    matConverter.fullToPodArray(xmatTrans, x);

    // Assemble global volume node numbering.
    nodeIds = ii + NumLocalPoints*jj;

    vmM = 0*ii;
    vmP = 0*ii;

    index_type count=0;
    for (index_type k = 0; k < NumElements; k++) {
        for( index_type f = 0; f < NumFaces; f++ ) {
            vmM(count) = nodeIds(Fmsk(f), k);
            count++;
        }
    }

    count = 0;
    for (index_type k1=0; k1 < NumElements; k1++) {
        for (index_type f1=0; f1 < NumFaces; f1++) {
            index_type k2 = E2E(k1, f1);
            index_type f2 = E2F(k1, f1);

            index_type vidM = vmM(k1*NumFaces + f1);
            index_type vidP = vmM(k2*NumFaces + f2);

            real_type dx = x[vidM] - x[vidP];
            real_type dist = sqrt(dx * dx);

            if ( dist < NodeTol) {
                vmP(count) = vidP;
            }
            count++;
        }
    }

    delete[] x;
}

/**
 *  Build Fmask. Mask that when applied to volume nodes gives the surface nodes.
 */
void Nodes1DProvisioner::buildFaceMask() {
    Array<double, 2> & x = *xGrid;
    Array<double, 2> & Fxref = *Fx;
    Array<int, 1> & Fmaskref = *Fmask;

    Fmaskref = 0, (NumLocalPoints - 1);

    for (int k = 0;  k < NumElements; k++) {
        for (int f = 0; f < NumFacePoints*NumFaces; f++) {
            Fxref(f, k) = x(Fmaskref(f), k);
        }
    }
}


/**
 * Build global connectivity matrices (EToE, EToF) for 1D grid
 * based using EToV (Element-to-Vertex) matrix.
 */
void Nodes1DProvisioner::buildConnectivityMatrices() {

    firstIndex ii;
    secondIndex jj;
    thirdIndex kk;

    int totalFaces = NumFaces*NumElements;
    int numVertices = NumElements + 1;

    int localVertNum[2];
    localVertNum[0] = 0; localVertNum[1] = 1;

    // Build global face-to-vertex array. (should be sparse matrix in 2D/3D).
    Array<double, 2> FToV(totalFaces, numVertices);
    FToV = 0*jj;

    Array<int, 2> & E2V = *EToV;

    int globalFaceNum = 0;
    for (int k=0; k < NumElements; k++) {
        for (int f=0; f < NumFaces; f++) {
            int v = localVertNum[f];
            int vGlobal = E2V(k,v);
            FToV(globalFaceNum, vGlobal) = 1;
            globalFaceNum++;
        }
    }

    Array<double, 2> FToF(totalFaces, totalFaces);
    Array<double, 2> I(totalFaces, totalFaces);

    FToF = 0*jj;
    I = 0*jj;

    for (int f=0; f < totalFaces; f++)
        I(f,f) = 1;

    // Global Face-to-Face connectivity matrix.
    FToF = sum(FToV(ii,kk)*FToV(jj,kk), kk) - I;

    Array<int,1> f1(totalFaces - 2); // '- 2' => for physical boundaries.
    Array<int,1> f2(totalFaces - 2);

    f1 = 0*ii;
    f2 = 0*ii;

    int connectionsCount = 0;
    for (int i=0; i < totalFaces; i++) {
        for (int j=0; j < totalFaces; j++) {
            if (FToF(i,j) == 1) {
                f1(connectionsCount) = i;
                f2(connectionsCount) = j;
                connectionsCount++;
            }
        }
    }

    Array<int, 1> e1(totalFaces - 2);
    Array<int, 1> e2(totalFaces - 2);

    // Convert face global number to local element and face numbers.
    e1 = floor(f1 / NumFaces);
    f1 = (f1 % NumFaces);
    e2 = floor(f2 / NumFaces);
    f2 = (f2 % NumFaces);

    // Build connectivity matrices.
    Array<int, 2> & E2E = *EToE;
    Array<int, 2> & E2F = *EToF;
    for (int k = 0; k < NumElements; k++) {
        for (int f = 0; f < NumFaces; f++) {
            E2E(k, f) = k;
            E2F(k, f) = f;
        }
    }

    for (int i=0; i < totalFaces - 2; i++) {
        int ee1 = e1(i);
        int ee2 = e2(i);
        int ff1 = f1(i);
        int ff2 = f2(i);
        E2E(ee1, ff1) = ee2;
        E2F(ee1, ff1) = ff2;
    }
}

void Nodes1DProvisioner::buildLift() {
    int Np = NumLocalPoints;
    firstIndex ii;
    secondIndex jj;
    thirdIndex kk;

    Array<double, 2> E(Np, NumFaces*NumFacePoints);
    E = 0*jj;
    E(0, 0)  = 1.;
    E(Np-1, 1) = 1.;

    Array<double, 2> & Vref = *V;
    Array<double, 2> & L = *Lift;

    Array<double, 2> Vtrans(Np, Np);
    Vtrans = Vref(jj,ii);
    Array<double, 2> temp(Np, NumFaces*NumFacePoints);
    temp = sum(Vtrans(ii,kk)*E(kk,jj), kk);
    L = sum(Vref(ii,kk)*temp(kk,jj), kk);
}

/**
 * Compute Jacobian (determinant) J and geometric factor rx (dr/dx), and Fscale using nodes and differentiation matrix.
 */
void Nodes1DProvisioner::computeJacobian() {
    firstIndex ii;
    secondIndex jj;
    thirdIndex kk;

    Array<double,2> & x = get_xGrid();
    Array<double,2> & Dr = get_Dr();
    Array<double,2> & Jref = get_J();
    Array<double,2> & rxref = get_rx();

    matrix_type & Fscaleref = *Fscale;

    // is Fmask not intiialized at this point?
    index_vector_type & Fmsk = get_Fmask();

    Jref = sum(Dr(ii,kk)*x(kk,jj), kk);
    rxref = 1/Jref;

    for(index_type f=0; f < NumFaces; f++) {
        Fscaleref(f, Range::all()) = 1/Jref(Fmsk(f), Range::all());
    }
}

/**
 * Compute Vandermonde matrix which maps modal coefficients to nodal values.
 */
void Nodes1DProvisioner::buildVandermondeMatrix() {
    V = new Array<double, 2>(NOrder+1, NOrder+1);

    Array<double, 2> & Vref = *V;

    for (int j=1; j <= NOrder+1; j++) {
        Array<double, 1> p(NOrder+1);
        computeJacobiPolynomial(*rGrid, 0.0, 0.0, j-1, p);
        Vref(Range::all(), j-1) = p;
    }
}

/**
 * Build differentiation matrix Dr on the standard element.
 */
void Nodes1DProvisioner::buildDr() {
    firstIndex ii;
    secondIndex jj;

    Dr = new Array<double, 2>(NOrder+1, NOrder+1);

    Array<double, 2> & Vref = *V;
    Array<double, 2> & Drref = *Dr;

    Array<double, 2> DVr(NOrder+1, NOrder+1);
    DVr = 0.*jj;

    computeGradVandermonde(DVr);

    // Dr = DVr / V;

    Array<double, 2> Vtrans(NOrder+1, NOrder+1);
    Array<double, 2> DVrtrans(NOrder+1, NOrder+1);
    Array<double, 2> Drtrans(NOrder+1, NOrder+1);

    Vtrans = Vref(jj, ii);
    DVrtrans = DVr(jj, ii);

    DirectSolver & linSolver = *LinSolver;
    linSolver.solve(Vtrans, DVrtrans,  Drtrans);

    Drref = Drtrans(jj, ii);
} 

/**
 * Get number of elements.
 */
int Nodes1DProvisioner::get_NumElements() {
    return NumElements;
}

/**
 * Get reference to 1D Lifting Operator.
 */
Array<double, 2> & Nodes1DProvisioner::get_Lift() {
    return *Lift;
}

/**
 * Get reference to physical x-grid.
 */
Array<double, 2> & Nodes1DProvisioner::get_xGrid() {
    return *xGrid;
}

/**
 * Get reference to r-grid on the standard element.
 */
Array<double, 1> & Nodes1DProvisioner::get_rGrid() {
    return *rGrid;
}

/**
 * Get reference to Element-to-Vertex connectivity table.
 */
Array<int, 2> & Nodes1DProvisioner::get_EToV() {
    return *EToV;
}

/**
 * Get reference to Element-to-Element connectivity table.
 */
Array<int, 2> & Nodes1DProvisioner::get_EToE() {
    return *EToE;
}

/**
 * Get reference to Element-to-Face connectivity table.
 */
Array<int, 2> & Nodes1DProvisioner::get_EToF() {
    return *EToF;
}

/**
 * Get reference to differentiation matrix Dr on the standard element.
 */
Array<double, 2> & Nodes1DProvisioner::get_Dr() {
    return *Dr;
}

/**
 * Get reference to generalized Vandermonde matrix V.
 */
Array<double, 2> & Nodes1DProvisioner::get_V() {
    return *V;
}

/**
 * Get reference to Jacobian scaling array J.
 */
Array<double, 2> & Nodes1DProvisioner::get_J() {
    return *J;
}

/**
 * Get reference to geometric scaling array rx.
 */
Array<double, 2> & Nodes1DProvisioner::get_rx() {
    return *rx;
}

/**
 * Get reference to normals array nx.
 */
const matrix_type & Nodes1DProvisioner::get_nx() {
    return *nx;
} 

int Nodes1DProvisioner::get_NumLocalPoints() {
    return NumLocalPoints;
}

/**
 * Get the faces-only x-grid.
 */
Array<double, 2> & Nodes1DProvisioner::get_Fx() {
    return *Fx;
}

/**
 * Get the Face-scaling factor (Inverse of Jacobian at Face nodes).
 */

const matrix_type & Nodes1DProvisioner::get_Fscale() {
    return *Fscale;
}

/**
 * Get the index-mask for the face nodes.
 */
Array<int, 1> & Nodes1DProvisioner::get_Fmask() {
    return *Fmask;
}

/**
 * Get the volume to surface map, 'minus' traces.
 */
const index_vector_type & Nodes1DProvisioner::get_vmapM() {
    return *vmapM;
}

/**
 * Get the volume to surface map, 'plus' traces.
 */
const index_vector_type & Nodes1DProvisioner::get_vmapP() {
    return *vmapP;
}

/**
 * Get the surface index of the inflow boundary.
 */
const index_type Nodes1DProvisioner::get_mapI() {
    return mapI;
}

const index_type Nodes1DProvisioner::get_mapO() {
    return mapO;
}

const index_type Nodes1DProvisioner::get_vmapI() {
    return vmapI;
}

const index_type Nodes1DProvisioner::get_vmapO() {
    return vmapO;
}

/**
 * Destructor.
 */
Nodes1DProvisioner::~Nodes1DProvisioner() {
    delete rGrid;
    delete xGrid;
    delete J;
    delete rx;
    delete Lift;
    delete EToV;
    delete EToE;
    delete EToF;
    delete Fmask;
    delete Fx;
    delete Fscale;
    delete nx;
    delete vmapM;
    delete vmapP;
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
    J = 0.;
    for (int i=0; i < N+1; i++) {
        double h1 = 2.*i+alpha+beta;
        J(i,i)   = -0.5*(alpha*alpha-beta*beta)/(h1+2.)/h1;
        if (i < N) {
            J(i,i+1) = 2./(h1+2)*sqrt((i+1)*((i+1)+alpha+beta)*((i+1)+alpha)*((i+1)+beta)/(h1+1)/(h1+3));
        }
    }

    if ((alpha + beta) < 10*eps ) {
        J(0,0) = 0.0;
    }
    J = J(ii,jj) + J(jj,ii);

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

void Nodes1DProvisioner::computeGradJacobi(Array<double,1> const & x, const double alpha, const double beta, const int N, Array<double,1> & dp) {
    if (N == 0) {
        dp = 0.0;
        return;
    }

    Array<double, 1> p(x.length());

    computeJacobiPolynomial(x, alpha+1, beta+1, N-1, p);
    dp = sqrt(N*(N+alpha+beta+1))*p;
}

void Nodes1DProvisioner::computeGradVandermonde(Array<double,2> & DVr) {

    firstIndex ii;
    for (int i=0; i<=NOrder; i++) {
        Array<double, 1> dp(NOrder+1);
        dp = 0.*ii;
        computeGradJacobi(*rGrid, 0.0, 0.0, i, dp);
        DVr(Range::all(), i) = dp;
    }
}