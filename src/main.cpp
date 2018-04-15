// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

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

  cout << nodes1DProvisioner.get_V() << endl;

  Array<double, 2> DVr(N+1, N+1);

  nodes1DProvisioner.computeGradVandermonde(DVr);

  cout << DVr << endl;
  cout << nodes1DProvisioner.get_Dr() << endl;

  cout << nodes1DProvisioner.get_xGrid() << endl;

  Array<double,2> J(N+1, K);
  Array<double,2> rx(N+1, K);
  nodes1DProvisioner.computeJacobian(J, rx);
  cout << J << endl;
  cout << rx << endl;

  return 0;
}