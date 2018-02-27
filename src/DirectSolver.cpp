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
 * Solve A*X=B using LAPACK. Here, B and X are allowed to have multiple columns.
 */
void DirectSolver::solve(const Array<double,2> & A, const Array<double, 2> & B, Array<double, 2> & X) {

    firstIndex ii;
    secondIndex jj;

    int sz = A.rows();
    int Nrhs = B.cols();

    int dim = sz*Nrhs;

    int lda = sz;
    int ldb = sz; 
    int ldx = sz;

    int ipiv[sz];

    double work[sz*Nrhs];
    float swork[sz*(sz+Nrhs)];

    int info;
    int iter;

    double Apod[sz*lda];
    double Bpod[dim];
    double Xpod[dim];

    Array<double, 2> Atrans(sz, sz);
    Array<double, 2> Btrans(Nrhs, sz);
    Array<double, 2> Xtrans(Nrhs, sz);

    Atrans = A(jj,ii);
    Btrans = B(jj,ii);

    MatrixConverter.fullToPodArray(Atrans, Apod);
    MatrixConverter.fullToPodArray(Btrans, Bpod);

    dsgesv_(&sz, &Nrhs, Apod, &lda,
             ipiv, Bpod, &ldb, Xpod, &ldx, 
             work, swork, &iter, &info);

    MatrixConverter.podArrayToFull(Xpod, Xtrans);

    X = Xtrans(jj,ii);
}

DirectSolver::~DirectSolver() {

}