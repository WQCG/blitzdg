// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "SparseTriplet.hpp"
#include "LinAlgHelpers.hpp"
#include <blitz/array.h>
#include <igloo/igloo_alt.h>
#include <iostream>

using blitz::firstIndex;
using blitz::secondIndex;
using std::cout;
using std::endl;

namespace blitzdg {
    namespace SparseTripletTests {
        using namespace igloo;

        matrix_type tripToFull(const SparseTriplet& trip) {
            matrix_type ret(trip.rows(), trip.cols());
            ret = real_type(0);
            for (index_type k = 0; k < trip.nnz(); ++k) {
                ret(trip.row(k), trip.col(k)) = trip.elem(k);
            }
            return ret;
        }

        Describe(SparseTriplet_Object) {
            It(Should_Build_From_A_Dense_Matrix) {
                matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                SparseTriplet trip(full);
                matrix_type ret = tripToFull(trip);
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
                matrix_type full(5,5);
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
                matrix_type ret = tripToFull(trip);
                ret -= full;
                real_type diff = normMax(ret);
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Build_From_Array_Input) {
                matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                index_vector_type r(11), c(11);
                vector_type e(11);
                r = 0, 0, 0, 1, 1, 2, 3, 3, 4, 4, 4;
                c = 0, 3, 4, 2, 3, 1, 1, 4, 0, 1, 2;
                e = 1, 3, 2, 1, 1, 5, 2, 4, 1, 1, 1;
                SparseTriplet trip(11, r.data(), c.data(), e.data());
                matrix_type ret = tripToFull(trip);
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
                matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                SparseTriplet trip(full);
                trip.clear();
                matrix_type ret = tripToFull(trip);
                real_type diff = normMax(ret);
                cout << "SparseTriplet: clear all elements" << endl;
                Assert::That(trip.nnz(), Equals(0));
                Assert::That(trip.nzmax(), Equals(11));
                Assert::That(diff, Equals(0.0));
            }
        };
    } // namespace SparseTripletTests
} // namespace blitzdg