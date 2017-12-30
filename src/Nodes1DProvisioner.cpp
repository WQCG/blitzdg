#include <iostream>
#include <math.h>
#include <Nodes1DProvisioner.hpp>

using namespace std;

/**
 * Constructor. Takes order of polynomials, number of elements, and dimensions of the domain.
 * Assumes equally-spaced elements.
 */
Nodes1DProvisioner::Nodes1DProvisioner(int NOrder, int NumElements, double xmin, double xmax) {
    throw("Not implemented.");
}

/**
 * Build nodes and geometric factors for all elements.
 */
void Nodes1DProvisioner::buildNodes() {
    throw("Not implemented.");
}

/**
 * Build differentiation matrix Dr on the standard element.
 */
void Nodes1DProvisioner::buildDr() {
    throw("Not implemented.");
} 

/**
 * Get reference to physical x-grid.
 */
Array<double, 2> & Nodes1DProvisioner::get_xGrid() {
    throw("Not implemented.");
}

/**
 * Get reference to r-grid on the standard element.
 */
Array<double, 1> & Nodes1DProvisioner::get_rGrid() {
    throw("Not implemented.");
}


/**
 * Get reference to differentiation matrix Dr on the standard element.
 */
Array<double, 2> & Nodes1DProvisioner::get_Dr() {
    throw("Not implemented.");
}

/**
 * Destructor
 */
Nodes1DProvisioner::~Nodes1DProvisioner() {
    throw("Not implemented.");
}