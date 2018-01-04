/** \mainpage blitzdg documentation
  *
  * \section intro_sec Introduction
  *
  * This is the auto-generated API docs page for the blitzdg project.
  *
  * Click on 'Classes' to start familarizing yourself with the API
  */

#include <iostream>
#include <blitz/array.h>
#include <MeshManager.hpp>
#include <Nodes1DProvisioner.hpp>
#include <SparseMatrixConverter.hpp>

using namespace std;
using namespace blitz;

int main(int argc, char **argv) {
	int N = 4;
	int K = 10;
	double xmin =-1.0;
	double xmax = 1.0;

  SparseMatrixConverter matrixConverter;
	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax, matrixConverter);
	
  nodes1DProvisioner.buildNodes();
  
  return 0;
}