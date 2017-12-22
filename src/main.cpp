#include <iostream>
#include <blitz/array.h>
#include <MeshManager.hpp>

using namespace std;
using namespace blitz;

int main(int argc, char **argv) {
	MeshManager mgr;

	//mgr.ReadMesh("foo");
	mgr.ReadVertices("/home/dsteinmo/blitzdg/input/2box.V");

	int dim = mgr.GetDim();
	int numVerts = mgr.GetNumVerts();

	cout << "dim: " << dim << endl;
	cout << "numVerts: " << numVerts << endl;

	double * & verts = mgr.GetVertices();

	// can move this into a printverts function or something.
	for(int i=0; i < numVerts; i++) {
		for (int j=0; j < dim ; j++) {
			cout << verts[mgr.GetIndex(i, j, dim)] << " ";
		}
		cout << endl;
	}


    return 0;
}