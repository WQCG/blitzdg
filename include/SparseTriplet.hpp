#pragma once
struct SparseTriplet {
    int nz;
    int *row;
    int *col;
    double *val;
};