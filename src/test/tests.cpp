#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <LUSolver.hpp>
#include <MeshManager.hpp>

using namespace igloo;
using namespace blitz;
using namespace std;

const int N=5;
const double eps=10*numeric_limits<double>::epsilon();

Array<double,2> A(N,N), B(N,N), C(N,N), D(N,N);
Array<double,1> b(N), soln(N), d(N), e(N), x(N);

firstIndex ii;
secondIndex jj;

LUSolver * luSolver = nullptr;
MeshManager * meshManager=nullptr;
SparseMatrixConverter *matrixConverter = nullptr;


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

int main(const int argc, const char *argv[])
{
  return TestRunner::RunAllTests(argc, argv);
  // exit code returns number of failed tests.
}
