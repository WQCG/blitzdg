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
#include <cmath>

using blitz::firstIndex;
using blitz::secondIndex;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::abs;

namespace blitzdg {
    namespace TriangleNodesProvisionerTests {
        using namespace igloo;

        firstIndex ii;
        secondIndex jj;

		MeshManager * meshManager = nullptr;
        Nodes1DProvisioner * nodes1DProvisioner = nullptr;
		TriangleNodesProvisioner * triangleNodesProvisioner = nullptr;

        Describe(Nodes1DProvisioner_Object) {
			const real_type eps = 20.*numeric_limits<double>::epsilon();
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

				Assert::That(abs(p(0) - 0.133252242007405), IsLessThan(eps));
				Assert::That(abs(p(1) - 0.355359724434270), IsLessThan(eps));
				Assert::That(abs(p(2) - 0.637112282097905), IsLessThan(eps));
            }

			It(Should_Map_rs_Coords_To_ab) {
                cout << "Should_Map_rs_Coords_To_ab" << endl;
                TriangleNodesProvisioner & triangleNodes = *triangleNodesProvisioner;

				real_vector_type r(3);
				real_vector_type s(3);

				r = -.1,.1,.2;
				s = .2,.3,.5;

				real_vector_type a(3);
				real_vector_type b(3);


				triangleNodes.rsToab(r, s, a, b);

				cout << "a: " << a << endl;
				cout << "b: " << b << endl;

				Assert::That(abs(b(0) - 0.2), IsLessThan(eps));
				Assert::That(abs(b(1) - 0.3), IsLessThan(eps));
				Assert::That(abs(b(2) - 0.5), IsLessThan(eps));

				Assert::That(abs(a(0) - 1.25), IsLessThan(eps));
				Assert::That(abs(a(1) - 2.14285714285714), IsLessThan(eps));
				Assert::That(abs(a(2) - 3.8), IsLessThan(eps));
            }

			It(Should_Compute_Warp_Factor) {
				cout << "Should_Compute_warpFactor" << endl;
                TriangleNodesProvisioner & triangleNodes = *triangleNodesProvisioner;

				real_vector_type r(3);
				real_vector_type warp(3);
				r = -.1,.1,.2;
				warp = 0.0,0.0,0.0;

				triangleNodes.computeWarpFactor(r, warp);

				Assert::That(abs(warp(0) - -0.0384345884812357), IsLessThan(eps));
				Assert::That(abs(warp(1) -  0.0384345884812359), IsLessThan(eps));
				Assert::That(abs(warp(2) -  0.0768691769624717), IsLessThan(eps));

				cout << warp(0) << endl;
				cout << warp(1) << endl;
				cout << warp(2) << endl;
			}
		};
   } // namespace Nodes1DProvisionerTests
} // namespace blitzdg
