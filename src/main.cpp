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


void computeRHS(const Array<double,2> & u, const double c, Nodes1DProvisioner & nodes1D, Array<double,2> & RHS) {
  //something like: RHS = -c*rx*(Dr*u) + c*LIFT*Fscale*numFlux
  Array<double,2> & Dr = nodes1D.get_Dr();
}

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

  // Blitz indices
  firstIndex ii;
  secondIndex jj;

  // Build dependencies.
  SparseMatrixConverter matrixConverter;
  EigenSolver eigenSolver(matrixConverter);
  DirectSolver directSolver(matrixConverter);
	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax, matrixConverter, eigenSolver, directSolver);
	
  // Pre-processing. Build grid, and get data we need for initialization.
  nodes1DProvisioner.buildNodes();
  nodes1DProvisioner.computeJacobian();

  int Np = nodes1DProvisioner.get_NumLocalPoints();

  Array<double,2> & x = nodes1DProvisioner.get_xGrid();

  double min_dx = x(1,0) - x(0,0);

  double dt = CFL*min_dx/c;

  Array<double, 2> u(Np, K);
  Array<double, 2> RHS(Np, K);

  u = exp(-5*(x(ii)*x(ii)));
  while (t < finalTime) {

    computeRHS(u, c, nodes1DProvisioner, RHS);
    // Forward Euler time-step for now, to be replaced.
    u = u + dt*RHS;

    cout << u << endl;
    t += dt;
  }
  
  return 0;
}