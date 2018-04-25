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

  int Np = nodes1DProvisioner.get_NumLocalPoints();

  Array<double,1> rGrid = nodes1DProvisioner.get_rGrid();
  cout << rGrid << endl;

  cout << nodes1DProvisioner.get_V() << endl;

  Array<double, 2> DVr(Np, Np);

  nodes1DProvisioner.computeGradVandermonde(DVr);

  cout << DVr << endl;
  cout << nodes1DProvisioner.get_Dr() << endl;

  Array<double,2> x = nodes1DProvisioner.get_xGrid();

  cout << x << endl;

  Array<double,2> J(Np, K);
  Array<double,2> rx(Np, K);

  nodes1DProvisioner.computeJacobian(J, rx);

  cout << "Can_Partition_A_Mesh" << endl;
  MeshManager * manager = new MeshManager();
  MeshManager & mgr = *manager;
  mgr.readVertices("input/2box.V");
  mgr.readElements("input/2box.E2V");

  mgr.partitionMesh(2);

  //int * & epMap = mgr.get_ElementPartitionMap();
  //int * & vpMap = mgr.get_VertexPartitionMap();

  delete manager;
  return 0;
}