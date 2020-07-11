// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "Traits.hpp"
#include <igloo/igloo_alt.h>
#include <iostream>
#include <iterator>
#include <forward_list>
#include <list>
#include <vector>

using std::cout;
using std::endl;
using std::forward_list;
using std::istream_iterator;
using std::list;
using std::ostream_iterator;
using std::vector;

namespace blitzdg {
    namespace TraitsTests {
        using namespace igloo;

        Describe(Traits_Object) {
            It(Identifies_Integral_Types) {
                cout << "Traits: identify integral types" << endl;
                Assert::That(isIntegral<char>(), Equals(true));
                Assert::That(isIntegral<short>(), Equals(true));
                Assert::That(isIntegral<int>(), Equals(true));
                Assert::That(isIntegral<long>(), Equals(true));
                Assert::That(isIntegral<long long>(), Equals(true));
                Assert::That(isIntegral<unsigned char>(), Equals(true));
                Assert::That(isIntegral<signed char>(), Equals(true));
                Assert::That(isIntegral<unsigned short>(), Equals(true));
                Assert::That(isIntegral<unsigned int>(), Equals(true));
                Assert::That(isIntegral<unsigned long>(), Equals(true));
                Assert::That(isIntegral<unsigned long long>(), Equals(true));
                Assert::That(isIntegral<float>(), Equals(false));
                Assert::That(isIntegral<double>(), Equals(false));
            }

            It(Identifies_Real_Types) {
                cout << "Traits: identify real types" << endl;
                Assert::That(isFloatingPoint<float>(), Equals(true));
                Assert::That(isFloatingPoint<double>(), Equals(true));
            }

            It(Identifies_Numeric_Types) {
                cout << "Traits: identify numeric types" << endl;
                Assert::That(isArithmetic<char>(), Equals(true));
                Assert::That(isArithmetic<short>(), Equals(true));
                Assert::That(isArithmetic<int>(), Equals(true));
                Assert::That(isArithmetic<long>(), Equals(true));
                Assert::That(isArithmetic<long long>(), Equals(true));
                Assert::That(isArithmetic<unsigned char>(), Equals(true));
                Assert::That(isArithmetic<signed char>(), Equals(true));
                Assert::That(isArithmetic<unsigned short>(), Equals(true));
                Assert::That(isArithmetic<unsigned int>(), Equals(true));
                Assert::That(isArithmetic<unsigned long>(), Equals(true));
                Assert::That(isArithmetic<unsigned long long>(), Equals(true));
                Assert::That(isArithmetic<float>(), Equals(true));
                Assert::That(isArithmetic<double>(), Equals(true));
            }

            It(Identifies_Iterator_Value_Types_Same) {
                cout << "Traits: iterator value types are the same" << endl;
                using dItr = vector<double>::iterator;
                using const_dItr = vector<double>::const_iterator;
                using iItr = vector<int>::iterator;
                Assert::That(iteratorValueTypesSame<dItr, dItr>(), Equals(true));
                Assert::That(iteratorValueTypesSame<const_dItr, dItr>(), Equals(true));
                Assert::That(iteratorValueTypesSame<iItr, dItr>(), Equals(false));
            }

            It(Identifies_Iterator_Value_Type_Same_As_Type) {
                cout << "Traits: iterator value type is same as type" << endl;
                using dItr = vector<double>::iterator;
                using const_dItr = vector<double>::const_iterator;
                using iItr = vector<int>::iterator;
                Assert::That(iteratorValueTypeSameAs<dItr, double>(), Equals(true));
                Assert::That(iteratorValueTypeSameAs<const_dItr, double>(), Equals(true));
                Assert::That(iteratorValueTypeSameAs<iItr, double>(), Equals(false));
            }

            It(Identifies_Iterator_Categories) {
                using raItr = vector<int>::iterator;
                using bdItr = list<int>::iterator;
                using fwItr = forward_list<int>::iterator;
                using inItr = istream_iterator<int>::iterator;
                using otItr = ostream_iterator<int>::iterator;
                cout << "Traits: is a random-access iterator" << endl;
                Assert::That(randomAccessIterator<raItr>(), Equals(true));
                Assert::That(randomAccessIterator<bdItr>(), Equals(false));
                Assert::That(randomAccessIterator<fwItr>(), Equals(false));
                Assert::That(randomAccessIterator<inItr>(), Equals(false));
                Assert::That(randomAccessIterator<otItr>(), Equals(false));
                cout << "Traits: is a bidirectional iterator" << endl;
                Assert::That(bidirectionalIterator<bdItr>(), Equals(true));
                Assert::That(bidirectionalIterator<raItr>(), Equals(true));
                Assert::That(bidirectionalIterator<fwItr>(), Equals(false));
                Assert::That(bidirectionalIterator<inItr>(), Equals(false));
                Assert::That(bidirectionalIterator<otItr>(), Equals(false));
                cout << "Traits: is a forward iterator" << endl;
                Assert::That(forwardIterator<fwItr>(), Equals(true));
                Assert::That(forwardIterator<raItr>(), Equals(true));
                Assert::That(forwardIterator<bdItr>(), Equals(true));
                Assert::That(forwardIterator<inItr>(), Equals(false));
                Assert::That(forwardIterator<otItr>(), Equals(false));
                cout << "Traits: is a input iterator" << endl;
                Assert::That(inputIterator<inItr>(), Equals(true));
                Assert::That(inputIterator<raItr>(), Equals(true));
                Assert::That(inputIterator<bdItr>(), Equals(true));
                Assert::That(inputIterator<fwItr>(), Equals(true));
                Assert::That(inputIterator<otItr>(), Equals(false));
                cout << "Traits: is a output iterator" << endl;
                Assert::That(outputIterator<otItr>(), Equals(true));
                Assert::That(outputIterator<raItr>(), Equals(true));
                Assert::That(outputIterator<bdItr>(), Equals(true));
                Assert::That(outputIterator<fwItr>(), Equals(true));
                Assert::That(outputIterator<inItr>(), Equals(false));
            }
        };
    } // namespace TraitsTests
} // namespace blitzdg