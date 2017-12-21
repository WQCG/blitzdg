#include <iostream>
#include <blitz/array.h>
#include <MeshManager.hpp>

using namespace std;
using namespace blitz;

int main(int argc, char **argv) {
	MeshManager mgr;

	mgr.ReadMesh("foo");
	mgr.ReadMesh("/home/dsteinmo/blitzdg/input/2box.V", "foo");


    return 0;
}