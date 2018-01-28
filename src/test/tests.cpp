#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <MeshManager.hpp>
#include <Nodes1DProvisioner.hpp>
#include <LUSolver.hpp>
#include <EigenSolver.hpp>


using namespace igloo;
using namespace blitz;
using namespace std;

const int N=5;
const double eps=10*numeric_limits<double>::epsilon();

Array<double,2> A(N,N), B(N,N), C(N,N), D(N,N), Adiag(N,N), Asymmetric(N,N);
Array<double,1> b(N), soln(N), d(N), e(N), x(N);

firstIndex ii;
secondIndex jj;

LUSolver * luSolver = nullptr;
MeshManager * meshManager=nullptr;
SparseMatrixConverter * matrixConverter = nullptr;
Nodes1DProvisioner * nodes1DProvisioner = nullptr;
EigenSolver * eigenSolver = nullptr;


Describe(Simple_blitz_array_operations)
{
  void SetUp() {

    A = 2,3,0,0,0,
		    3,0,4,0,6,
		    0,-1,-3,2,0,
		    0,0,1,0,0,
		    0,4,2,0,1;

	  B = jj;

    b =  8,
        45,
        -3,
         3,
        19;

    x = 1,
        2,
        3,
        4,
        5;

    d = 1,2,3,4;
    e = 2,3,4,5;
  }

  It(Properly_Multiplies_Pointwise)
  {
    Array<double, 1> result(N);
    result = d(ii)*e(ii);

    Assert::That(result(0), Equals(2.));
    Assert::That(result(1), Equals(6.));
    Assert::That(result(2), Equals(12.));
    Assert::That(result(3), Equals(20.));
  }

  It(Properly_Does_Dot_Product)
  {
    double result;
    result = sum(d * e);

    Assert::That(result, Equals(40.));
  }
  It(Properly_Does_Matrix_Vector_Product)
  {
    Array <double, 1> result(N);
    result = sum(A(ii,jj)*x(jj), jj);

    Assert::That(result(0), Equals(b(0)));
    Assert::That(result(1), Equals(b(1)));
    Assert::That(result(2), Equals(b(2)));
    Assert::That(result(3), Equals(b(3)));
    Assert::That(result(4), Equals(b(4)));
  }
};

Describe(LUSolver_Object)
{
  void SetUp() {

    A = 2,3,0,0,0,
		    3,0,4,0,6,
		    0,-1,-3,2,0,
		    0,0,1,0,0,
		    0,4,2,0,1;

    b =  8,
        45,
        -3,
         3,
        19;

    x = 1,
        2,
        3,
        4,
        5;

    matrixConverter = new SparseMatrixConverter();
    luSolver = new LUSolver(&A, *matrixConverter);

  }
  It(Solves_Ax_equals_b) 
  {
    LUSolver & solver = *luSolver;
    //LUSolver luSolver(&A);
    Array <double, 1> soln(N);

	  // Compute LU factors.
    cout << endl;
	  solver.factorize();
    solver.solve(b, soln);
    
    Assert::That(abs(soln(0)-x(0)), IsLessThan(eps));
    Assert::That(abs(soln(1)-x(1)), IsLessThan(eps));
    Assert::That(abs(soln(2)-x(2)), IsLessThan(eps));
    Assert::That(abs(soln(3)-x(3)), IsLessThan(eps));
    Assert::That(abs(soln(4)-x(4)), IsLessThan(eps));
  }
};

Describe(MeshManager_Object) {
  void SetUp() {
    meshManager = new MeshManager();
  }

  It(Reads_Vertex_Files) {
    MeshManager & mgr = *meshManager;

    mgr.readVertices("input/2box.V");

    Assert::That(mgr.get_NumVerts(), Equals(6));
    Assert::That(mgr.get_Dim(), Equals(2));

    double * verts = mgr.get_Vertices();

    Assert::That(verts[0], Equals(0.0));
    Assert::That(verts[1], Equals(0.0));

    Assert::That(verts[2], Equals(0.5));
    Assert::That(verts[3], Equals(0.0));

    Assert::That(verts[4], Equals(1.0));
    Assert::That(verts[5], Equals(0.0));

    Assert::That(verts[6], Equals(1.0));
    Assert::That(verts[7], Equals(1.0));

    Assert::That(verts[8], Equals(0.5));
    Assert::That(verts[9], Equals(1.0));

    Assert::That(verts[10], Equals(0.0));
    Assert::That(verts[11], Equals(1.0));
  }

  It(Reads_Element_Files) {
    MeshManager & mgr = *meshManager;

    mgr.readElements("input/2box.E2V");

    Assert::That(mgr.get_NumElements(), Equals(2));
    Assert::That(mgr.get_ElementType(), Equals(4));

    int * elements = mgr.get_Elements();

    Assert::That(elements[0], Equals(1));
    Assert::That(elements[1], Equals(2));
    Assert::That(elements[2], Equals(5));
    Assert::That(elements[3], Equals(6));

    Assert::That(elements[4], Equals(2));
    Assert::That(elements[5], Equals(3));
    Assert::That(elements[6], Equals(4));
    Assert::That(elements[7], Equals(5));
  }

  It(Can_Print_Vertices_And_DoesNotThrow) {
    MeshManager & mgr = *meshManager;
    mgr.readVertices("input/2box.V");
    cout << endl << "Vertices:" << endl;
    mgr.printVertices();
  }

  It(Can_Print_Elements_And_DoesNotThrow) {
    MeshManager & mgr = *meshManager;
    mgr.readElements("input/2box.E2V");
    cout << endl << "Elements" << endl;
    mgr.printElements();
  }

  It(Can_Partition_A_Mesh) {
    MeshManager & mgr = *meshManager;
    mgr.readVertices("input/2box.V");
    mgr.readElements("input/2box.E2V");

    mgr.partitionMesh(2);

    int * & epMap = mgr.get_ElementPartitionMap();
    Assert::That(epMap[0], Equals(1));
    Assert::That(epMap[1], Equals(2));

    int * & vpMap = mgr.get_VertexPartitionMap();
    Assert::That(vpMap[0], Equals(1));
    Assert::That(vpMap[1], Equals(1));
    Assert::That(vpMap[2], Equals(2));
    Assert::That(vpMap[3], Equals(2));
    Assert::That(vpMap[4], Equals(1));
    Assert::That(vpMap[5], Equals(2));
  }
};

Describe(EigenSolver_Object) {
  void SetUp() {

    Adiag = 1,0,0,0,0,
        0,2,0,0,0,
        0,0,3,0,0,
        0,0,0,4,0,
        0,0,0,0,5;

    Asymmetric = 0,0.57735,0,0,0,
                 0.57735,-0,0.516398,0,0,
                 0,0.516398,-0,0.507093,0,
                 0,0,0.507093,-0,0.503953,
                 0,0,0,0.503953,-0;

    matrixConverter = new SparseMatrixConverter();
  }

  It(Should_Solve_Trivial_Symmetric_Eigenproblem) {
    eigenSolver = new EigenSolver(*matrixConverter);
    EigenSolver & solver = *eigenSolver;

    Array<double, 1> eigenvalues(5);
    eigenvalues = 0,0,0,0,0;
    Array<double, 2> eigenvectors(5,5);
    eigenvectors = 0,0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0;

    solver.solve(Adiag, eigenvalues, eigenvectors);
    
    Array<double, 2> expectedEvecs(5,5);
    
    expectedEvecs = 1,0,0,0,0,
                    0,1,0,0,0,
                    0,0,1,0,0,
                    0,0,0,1,0,
                    0,0,0,0,1;

    Assert::That(eigenvalues(0), Equals(1.));
    Assert::That(eigenvalues(1), Equals(2.));
    Assert::That(eigenvalues(2), Equals(3.));
    Assert::That(eigenvalues(3), Equals(4.));
    Assert::That(eigenvalues(4), Equals(5.));

    Array <double, 2> res(5,5);
    res = eigenvectors - expectedEvecs;
    Assert::That(sum(res(ii)*res(ii)), IsLessThan(eps));
  }

  It(Should_Solve_NonTrivial_Symmetric_Eigenproblem) {
    eigenSolver = new EigenSolver(*matrixConverter);
    EigenSolver & solver = *eigenSolver;

    Array<double, 1> eigenvalues(5);
    eigenvalues = 0,0,0,0,0;
    Array<double, 2> eigenvectors(5,5);
    eigenvectors = 0,0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0,
                   0,0,0,0,0;

    cout << "Solving eigenproblem!";
    solver.solve(Asymmetric, eigenvalues, eigenvectors);
    cout << "eigenvectors: " << eigenvectors << endl;

    Array<double, 2> expectedEvecs(5,5);

    expectedEvecs = 0.344185,-0.540215,0.563165,-0.456254,0.253736,
                   -0.489198,0.456254,0.0711849,-0.540215,0.505587,
                    0.533334,2.06417e-16,-0.596285,-2.67144e-16,0.6,
                   -0.489198,-0.456254,0.0711849,0.540215,0.505587,
                    0.344185,0.540215,0.563165,0.456254,0.253736;

    const float epsf = 5.e-7;

    Assert::That(eigenvalues(0) - -0.90618,    IsLessThan(epsf));
    Assert::That(eigenvalues(1) - -0.538469,   IsLessThan(epsf));
    Assert::That(eigenvalues(2) - 9.62592e-17, IsLessThan(epsf));
    Assert::That(eigenvalues(3) - 0.538469,    IsLessThan(epsf));
    Assert::That(eigenvalues(4) - 0.90618,     IsLessThan(epsf));

    Array<double, 2> res(5,5);
    res = eigenvectors - expectedEvecs;
    Assert::That(sum(res(ii)*res(ii)), IsLessThan(epsf));
  }
};

Describe(Nodes1DProvisioner_Object) {
  void SetUp() {
    const int NOrder = 3;
    const int NumElements = 5;
    const double xmin = -1.0;
    const double xmax = 1.0;

    matrixConverter = new SparseMatrixConverter();
    eigenSolver = new EigenSolver(*matrixConverter);

    nodes1DProvisioner = new Nodes1DProvisioner(NOrder, NumElements, xmin, xmax, *matrixConverter, *eigenSolver);
  }

  It(Should_Generate_0th_Order_Legendre_Polynomial) {
    Array<double, 1> x(3);
    x = -1.,0.,1.;
    Array<double, 1>  p(3);

    Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

    nodes1D.computeJacobiPolynomial(x, 0.0, 0.0, 0, p);

    Assert::That(p(0), Equals(1./sqrt(2.)));
    Assert::That(p(1), Equals(1./sqrt(2.)));
    Assert::That(p(2), Equals(1./sqrt(2.)));
  }

  It(Should_Generate_1st_Order_Legendre_Polynomial) {
    Array<double, 1> x(3);
    x = -1.,0.,1.;
    Array<double, 1>  p(3);

    Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

    nodes1D.computeJacobiPolynomial(x, 0.0, 0.0, 1, p);

    Assert::That(p(0), Equals(-sqrt(3./2.)));
    Assert::That(p(1), Equals(0.0));
    Assert::That(p(2), Equals(sqrt(3./2.)));
  }


  It(Should_Generate_2nd_Order_Legendre_Polynomial) {
    Array<double, 1> x(3);
    x = -1.,0.,1.;
    Array<double, 1>  p(3);

    Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

    nodes1D.computeJacobiPolynomial(x, 0.0, 0.0, 2, p);

    Assert::That(abs(p(0)-sqrt(5./2.)), IsLessThan(eps));
    Assert::That(abs(p(1)-(-sqrt(5./8.))), IsLessThan(eps));
    Assert::That(abs(p(2)-sqrt(5./2.)), IsLessThan(eps));
  }
};

int main(const int argc, const char *argv[])
{
  return TestRunner::RunAllTests(argc, argv);
  // exit code returns number of failed tests.
}
