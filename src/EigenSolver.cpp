// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include <EigenSolver.hpp>
#include <blitz/array.h>

using namespace blitz;

/**
 * Constructor. Takes a reference to a SparseMatrixConverter.
 */
EigenSolver::EigenSolver(SparseMatrixConverter const & _matrixConverter) {
    MatrixConverter = _matrixConverter;
}

extern "C" {
    void dsyevd_( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* iwork, int* liwork, int* info );
}

/**
 * Solve Ax=λx using LAPACK. Eigenvalues are stored in reference 'eigenvalues' and eigenvectors are stored column-wise
 * in reference 'eigenvectors.'
 */
void EigenSolver::solve(const Array<double,2> & A, Array<double,1> & eigenvalues, Array<double, 2> & eigenvectors) {

    int sz = A.rows();
    int lda = sz;
    int iwkopt;

    double ww[sz];
    double wkopt;
    int lwork = -1;
    int liwork = -1;
    int info;

    char JOBZ = 'V';
    char UPLO[] = "UP";

    double * Apod = new double[sz*lda];

    MatrixConverter.fullToPodArray(A, Apod);

    /* Determining optimal workspace parameters */
    dsyevd_( &JOBZ, UPLO, &sz, Apod, &lda, ww, &wkopt, &lwork, &iwkopt, &liwork, &info );

    lwork = (int)wkopt;
    double * work = new double[lwork];
    liwork = iwkopt;
    int * iwork = new int[liwork];

    /* Solve eigenproblem */
    dsyevd_( &JOBZ, UPLO, &sz, Apod, &lda, ww, work, &lwork, iwork, &liwork, &info );

    MatrixConverter.podArrayToFull(Apod, eigenvectors);

    for (int i=0; i < sz; i++)
        eigenvalues(i) = ww[i];

    firstIndex ii;
    secondIndex jj;
    eigenvectors = eigenvectors(jj,ii);
}

EigenSolver::~EigenSolver() {

}