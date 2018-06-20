// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "LinAlgHelpers.hpp"
#include "TriangleNodesProvisioner.hpp"
#include "MeshManager.hpp"
#include "Types.hpp"
#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <iostream>
#include <limits>

using blitz::firstIndex;
using blitz::secondIndex;
using std::cout;
using std::endl;
using std::numeric_limits;

namespace blitzdg {
    namespace TriangleNodesProvisionerTests {
        using namespace igloo;

        firstIndex ii;
        secondIndex jj;

		MeshManager * meshManager = nullptr;
        Nodes1DProvisioner * nodes1DProvisioner = nullptr;
		TriangleNodesProvisioner * triangleNodesProvisioner = nullptr;

        Describe(Nodes1DProvisioner_Object) {
			const real_type eps = 5.*numeric_limits<double>::epsilon();
			const float epsf = 5.8e-5;
			const index_type NOrder = 3;
			const index_type NumElements = 5;
			const int NumFaces = 2;

            void SetUp() {
				meshManager = new MeshManager();
				triangleNodesProvisioner = new TriangleNodesProvisioner(NOrder, NumElements, meshManager);
			}

			void TearDown() {
				delete triangleNodesProvisioner;
				delete meshManager;
			}

			It(Should_Evaluate_Orthonormal_Simplex2D_Polynomial) {
                cout << "Should_Evaluate_Orthonormal_Simplex2D_Polynomial" << endl;
                TriangleNodesProvisioner & triangleNodes = *triangleNodesProvisioner;

				real_vector_type a(3);
				real_vector_type b(3);

				a = .1,.2,.3;
				b = .2,.3,.4;

				real_vector_type p(3);
				p = .0,.0,.0;

				triangleNodes.evaluateSimplexPolynomial(a, b, 1, 2, p);

				cout << p << endl;

				Assert::That(abs(p(0) - 0.133252242007405), IsLessThan(eps));
				Assert::That(abs(p(1) - 0.355359724434270), IsLessThan(eps));
				Assert::That(abs(p(2) - 0.637112282097905), IsLessThan(eps));
            }
		};
   } // namespace Nodes1DProvisionerTests
} // namespace blitzdg
