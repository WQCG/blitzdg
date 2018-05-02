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
	
  // Physical parameters
  double xmin =-1.0;
	double xmax = 1.0;
  double c = 0.1;

  const double finalTime = 10.0;
  double t = 0.0;

  // Numerical parameters:
  int N = 4;
	int K = 20;
  double CFL = 0.5;

  //Build dependencies.
  SparseMatrixConverter matrixConverter;
  EigenSolver eigenSolver(matrixConverter);
  DirectSolver directSolver(matrixConverter);
	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax, matrixConverter, eigenSolver, directSolver);
	
   
  nodes1DProvisioner.buildNodes();
  int Np = nodes1DProvisioner.get_NumLocalPoints();
  
  Array<double,2> J(Np, K);
  Array<double,2> rx(Np, K);
  nodes1DProvisioner.computeJacobian(J, rx);


  Array<double,2> & x = nodes1DProvisioner.get_xGrid();
  Array<double,2> & Dr = nodes1DProvisioner.get_Dr();
  Array<double,2> & V = nodes1DProvisioner.get_V();


  cout << Dr << endl;

  cout << V << endl;

  double min_dx = x(1,0) - x(0,0);

  double dt = CFL*min_dx/c;

  double u = 1;
  while (t < finalTime) {

    u = u - dt*c*u;

    cout << u << endl;
    t += dt;
  }
  
  return 0;
}