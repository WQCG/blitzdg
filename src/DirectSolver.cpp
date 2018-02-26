#include <DirectSolver.hpp>
#include <blitz/array.h>

using namespace blitz;

/**
 * Constructor. Takes a reference to a SparseMatrixConverter.
 */
DirectSolver::DirectSolver(SparseMatrixConverter const & _matrixConverter) {
    MatrixConverter = _matrixConverter;
}

extern "C" {
    void dsgesv_( int* n, int* nrhs, double* a, int* lda,
                int* ipiv, double* b, int* ldb, double* x, int* ldx, 
                double* work, float* swork, int* iter, int* info );
}

/**
 * Solve AX=B using LAPACK. Here, B and X are allowed to have multiple columns.
 */
void DirectSolver::solve(const Array<double,2> & A, const Array<double, 2> & B, Array<double, 2> & X) {

    int sz = A.rows();
    int Nrhs = B.cols();

    int dim = sz*Nrhs;

    int lda = sz;
    int ldb = Nrhs; 
    int ldx = ldb;

    int ipiv[sz];

    double work[dim];
    float swork[dim];

    int info;
    int iter;

    double * Apod = new double[sz*lda];
    double * Bpod = new double[sz*ldb];
    double * Xpod = new double[sz*ldx];

    MatrixConverter.fullToPodArray(A, Apod);
    MatrixConverter.fullToPodArray(B, Bpod);

    MatrixConverter.podArrayToFull(Xpod, X);
}

DirectSolver::~DirectSolver() {

}