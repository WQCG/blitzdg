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
#include <EigenSolver.hpp>
#include <DirectSolver.hpp>

using namespace std;
using namespace blitz;

int main(int argc, char **argv) {
	int N = 4;
	int K = 10;
	double xmin =-1.0;
	double xmax = 1.0;

  SparseMatrixConverter matrixConverter;
  EigenSolver eigenSolver(matrixConverter);
  DirectSolver directSolver(matrixConverter);

	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax, matrixConverter, eigenSolver, directSolver);
	
  nodes1DProvisioner.buildNodes();

  Array<double,1> rGrid = nodes1DProvisioner.get_rGrid();
  cout << rGrid << endl;

  nodes1DProvisioner.buildVandermondeMatrix();
  
  cout << nodes1DProvisioner.get_V() << endl;

  Array<double,1> dp(N+1);
  nodes1DProvisioner.computeGradJacobi(rGrid, 0.0, 0.0, N, dp);
  cout << dp << endl;

  Array<double, 2> DVr(N+1, N+1);

  nodes1DProvisioner.computeGradVandermonde(DVr);

  cout << DVr << endl;

  nodes1DProvisioner.buildDr();

  return 0;
}