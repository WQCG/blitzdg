#include <iostream>
#include <blitz/array.h>
#include <MeshManager.hpp>

using namespace std;
using namespace blitz;

int main(int argc, char **argv) {
	MeshManager mgr;

	mgr.readVertices("input/2box.V");

	int dim = mgr.get_Dim();
	int numVerts = mgr.get_NumVerts();

	cout << "dim: " << dim << endl;
	cout << "numVerts: " << numVerts << endl;

	mgr.printVertices();

	mgr.readElements("input/2box.E2V");
	int numElements = mgr.get_NumElements();
	int elementType = mgr.get_ElementType();

	cout << "numElements: " << numElements << endl;
	cout << "elementType: " << elementType << endl;

	mgr.printElements();

	mgr.partitionMesh(2);
    return 0;
}