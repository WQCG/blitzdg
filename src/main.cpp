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


	double * & verts = mgr.get_Vertices();

	// can move this into a printverts function or something.
	for(int i=0; i < numVerts; i++) {
		for (int j=0; j < dim ; j++) {
			cout << verts[mgr.get_Index(i, j, dim)] << " ";
		}
		cout << endl;
	}

	mgr.readElements("input/2box.E2V");
	int numElements = mgr.get_NumElements();
	int elementType = mgr.get_ElementType();

	cout << "numElements: " << numElements << endl;
	cout << "elementType: " << elementType << endl;

	int * & elements = mgr.get_Elements();

	// can move this into a printElements function or something.
	for(int i=0; i < numElements; i++) {
		for (int j=0; j < elementType; j++) {
			cout << elements[mgr.get_Index(i, j, elementType)] << " ";
		}
		cout << endl;
	}

	mgr.partitionMesh(2);
    return 0;
}