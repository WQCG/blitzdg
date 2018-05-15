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

#include "MeshManager.hpp"
#include "Nodes1DProvisioner.hpp"
#include "SparseMatrixConverter.hpp"
#include "EigenSolver.hpp"
#include "DirectSolver.hpp"
#include "Types.hpp"
#include <iostream>

using std::cout;
using std::endl;

int main(int argc, char **argv) {
  using namespace blitzdg;
	index_type N = 4;
	index_type K = 10;
	real_type xmin =-1.0;
	real_type xmax = 1.0;

  SparseMatrixConverter matrixConverter;
  EigenSolver eigenSolver(matrixConverter);
  DirectSolver directSolver(matrixConverter);

	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax, matrixConverter, eigenSolver, directSolver);
	
  nodes1DProvisioner.buildNodes();
  nodes1DProvisioner.computeJacobian();

  index_type Np = nodes1DProvisioner.get_NumLocalPoints();

  vector_type rGrid = nodes1DProvisioner.get_rGrid();
  cout << rGrid << endl;

  cout << nodes1DProvisioner.get_V() << endl;

  matrix_type DVr(Np, Np);

  nodes1DProvisioner.computeGradVandermonde(DVr);

  cout << DVr << endl;
  cout << nodes1DProvisioner.get_Dr() << endl;

  matrix_type x = nodes1DProvisioner.get_xGrid();

  cout << x << endl;

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