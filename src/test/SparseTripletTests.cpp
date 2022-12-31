// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "SparseTriplet.hpp"
#include "LinAlgHelpers.hpp"
#include <blitz/array.h>
#include <igloo/igloo_alt.h>
#include <iostream>
#include <sstream>
#include <string>

using blitz::firstIndex;
using blitz::secondIndex;
using std::cout;
using std::endl;
using std::ostringstream;
using std::string;

namespace blitzdg {
    namespace SparseTripletTests {
        using namespace igloo;

        real_matrix_type tripToFull(const SparseTriplet& trip) {
            real_matrix_type ret(trip.rows(), trip.cols());
            ret = real_type(0);
            for (index_type k = 0; k < trip.nnz(); ++k) {
                ret(trip.row(k), trip.col(k)) = trip.elem(k);
            }
            return ret;
        }

        Describe(SparseTriplet_Object) {
            It(Should_Build_From_A_Dense_Matrix) {
                real_matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                SparseTriplet trip(full);
                real_matrix_type ret = tripToFull(trip);
                ret -= full;
                real_type diff = normMax(ret);
                cout << "SparseTriplet: build from full matrix" << endl;
                Assert::That(trip.rows(), Equals(5));
                Assert::That(trip.cols(), Equals(5));
                Assert::That(trip.nnz(), Equals(11));
                Assert::That(trip.nzmax(), Equals(11));
                Assert::That(trip.empty(), Equals(false));
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Build_By_Inserting_Elements) {
                real_matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                SparseTriplet trip(5, 5, 4);
                cout << "SparseTriplet: build by insertion" << endl;
                trip.insert(0, 0, 1);
                trip.insert(0, 3, 3);
                trip.insert(0, 4, 2);
                trip.insert(1, 2, 1);
                trip.insert(1, 3, 1); // should reallocate
                Assert::That(trip.nzmax(), Equals(6));
                trip.insert(2, 1, 5);
                trip.insert(3, 1, 2); // should reallocate
                Assert::That(trip.nzmax(), Equals(9));
                trip.insert(3, 4, 4);
                trip.insert(4, 0 ,1);
                trip.insert(4, 1, 1); // should reallocate
                Assert::That(trip.nzmax(), Equals(13));
                trip.insert(4, 2, 1);
                Assert::That(trip.nnz(), Equals(11));
                real_matrix_type ret = tripToFull(trip);
                ret -= full;
                real_type diff = normMax(ret);
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Build_From_Array_Input) {
                real_matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                index_vector_type r(11), c(11);
                real_vector_type e(11);
                r = 0, 0, 0, 1, 1, 2, 3, 3, 4, 4, 4;
                c = 0, 3, 4, 2, 3, 1, 1, 4, 0, 1, 2;
                e = 1, 3, 2, 1, 1, 5, 2, 4, 1, 1, 1;
                SparseTriplet trip(11, r.data(), c.data(), e.data());
                real_matrix_type ret = tripToFull(trip);
                ret -= full;
                real_type diff = normMax(ret);
                cout << "SparseTriplet: build from arrays" << endl;
                Assert::That(trip.rows(), Equals(5));
                Assert::That(trip.cols(), Equals(5));
                Assert::That(trip.nnz(), Equals(11));
                Assert::That(trip.nzmax(), Equals(11));
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Expand_Its_Capacity) {
                SparseTriplet trip(5, 5, 5);
                cout << "SparseTriplet: expand its capacity" << endl;
                trip.expand(4);
                Assert::That(trip.nzmax(), Equals(5));
                trip.expand(10);
                Assert::That(trip.nzmax(), Equals(10));
            }

            It(Should_Clear_All_Its_Elements) {
                real_matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                SparseTriplet trip(full);
                trip.clear();
                real_matrix_type ret = tripToFull(trip);
                real_type diff = normMax(ret);
                cout << "SparseTriplet: clear all elements" << endl;
                Assert::That(trip.nnz(), Equals(0));
                Assert::That(trip.nzmax(), Equals(11));
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Swap_Two_Matrices) {
                real_matrix_type full(5,5), dens(3,3);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                dens = 1, 2, 3,
                       4, 5, 6,
                       7, 8, 9;
                SparseTriplet t1(full), t2(dens);
                using std::swap;
                swap(t1, t2);
                real_matrix_type ret1(3,3), ret2(5,5);
                ret1 = tripToFull(t1);
                ret2 = tripToFull(t2);
                ret1 -= dens;
                ret2 -= full;
                cout << "SparseTriplet: swap" << endl;
                Assert::That(normMax(ret1), Equals(0.0));
                Assert::That(normMax(ret2), Equals(0.0));
            }

            It(Should_Print_To_An_Output_Stream) {
                real_matrix_type full(2,2);
                full = 1,2,
                       3,4;
                string s = "rows = 2, cols = 2, nnz = 4\n\n0 0 1\n0 1 2\n1 0 3\n1 1 4\n";
                SparseTriplet t(full);
                ostringstream oss;
                oss << t;
                Assert::That((oss && oss.str() == s), Equals(true));
            }
        };
    } // namespace SparseTripletTests
} // namespace blitzdg