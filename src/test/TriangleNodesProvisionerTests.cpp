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
			const real_type eps = 5.e1*numeric_limits<double>::epsilon();
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

			It(Should_Map_xy_Coords_To_rs) {
                cout << "Should_Map_xy_Coords_To_rs" << endl;
                TriangleNodesProvisioner & triangleNodes = *triangleNodesProvisioner;

				real_vector_type x(3), y(3), r(3), s(3);

				x = 0.5,0.6,0.7;
				y = 0.2,0.3,0.4;

				triangleNodes.xyTors(x, y, r, s);

				Assert::That(abs(r(0) -  0.0511966128287416), IsLessThan(eps));
				Assert::That(abs(r(1) -  0.0934615859097789), IsLessThan(eps));
				Assert::That(abs(r(2) -  0.1357265589908162), IsLessThan(eps));

				Assert::That(abs(s(0) - -0.1023932256574831), IsLessThan(eps));
				Assert::That(abs(s(1) -  0.0130768281804420), IsLessThan(eps));
				Assert::That(abs(s(2) -  0.1285468820183672), IsLessThan(eps));
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

			It(Should_Compute_Vandermode_Matrix) {
				cout << "Should_Compute_Vandermonde_Matrix" << endl;
                TriangleNodesProvisioner & triangleNodes = *triangleNodesProvisioner;

				const int Np = (NOrder+1)*(NOrder+2)/2;

				real_vector_type r(Np);
				real_vector_type s(Np);

				r = .1,.2,.3,.4,.5,.6,.7,.8,.9,1.;
				s = .1,.2,.3,.4,.5,.6,.7,.8,.9,1.;

				real_matrix_type V(Np,Np), V_expected(Np,Np);
				V = 0*jj;

				triangleNodes.computeVandermondeMatrix(NOrder, r, s, V);

				V_expected =    0.70711,0.65000,-0.45928,-0.76279,1.12583,2.41300,1.19811,1.45831,4.79915,1.83014,
								0.70711,0.80000,-0.24495,-0.90510,1.38564,3.39411,2.66504,2.40998,8.90497,4.07092,
								0.70711,0.95000,0.03062,-0.92012,1.64545,4.53432,4.82274,3.53966,14.50972,7.36686,
								0.70711,1.10000,0.36742,-0.77075,1.90526,5.83363,7.78693,4.84734,21.82920,11.89473,
								0.70711,1.25000,0.76547,-0.41984,2.16506,7.29204,11.67335,6.33304,31.07926,17.83134,
								0.70711,1.40000,1.22474,0.16971,2.42487,8.90955,16.59774,7.99675,42.47571,25.35347,
								0.70711,1.55000,1.74526,1.03503,2.68468,10.68615,22.67585,9.83847,56.23439,34.63793,
								0.70711,1.70000,2.32702,2.21324,2.94449,12.62186,30.02340,11.85819,72.57111,45.86149,
								0.70711,1.85000,2.97001,3.74148,3.20429,14.71666,38.75613,14.05593,91.70170,59.20097,
								0.70711,2.00000,3.67423,5.65685,-0.00000,-0.00000,-0.00000,0.00000,0.00000,-0.00000;

				real_matrix_type res(Np,Np);
				res = V-V_expected;
				Assert::That(normFro(res), IsLessThan(epsf));
			}

			It(Should_Compute_Equilateral_Nodes) {
				cout << "Should_Compute_Equilateral_Nodes" << endl;

				TriangleNodesProvisioner & triangleNodes = *triangleNodesProvisioner;

				const int Np = (NOrder+1)*(NOrder+2)/2;

				real_vector_type x(Np), y(Np), expectedx(Np), expectedy(Np);
				triangleNodes.computeEquilateralNodes(x, y);

				expectedx = -1.000000000000000,
							-0.447213595499958,
							 0.447213595499958,
							 1.000000000000000,
							-0.723606797749979,
							-0.000000000000000,
							 0.723606797749979,
							-0.276393202250021,
							 0.276393202250021,
							 0.000000000000000;

				expectedy = -5.77350269189626e-01,
							-5.77350269189626e-01,
							-5.77350269189626e-01,
							-5.77350269189626e-01,
							-9.86232000259289e-02,
							-9.67499238300623e-17,
							-9.86232000259288e-02,
							 6.75973469215554e-01,
							 6.75973469215555e-01,
							 1.15470053837925e+00;

				real_vector_type resx(Np), resy(Np);

				resx = x - expectedx;
				resy = y - expectedy;

				Assert::That(normInf(resx), IsLessThan(eps));
				Assert::That(normInf(resy), IsLessThan(eps));
			}

			It(Should_Compute_Gradient_Of_Simplex_Polynomials) {
				cout << "Should_Compute_Gradient_Of_Simplex_Polynomials" << endl;

				TriangleNodesProvisioner & triangleNodes = *triangleNodesProvisioner;

				const int Np = (NOrder+1)*(NOrder+2)/2;

				real_vector_type a(Np), b(Np), dpdr(Np), dpds(Np), x(Np), y(Np), r(Np), s(Np);

				triangleNodes.computeEquilateralNodes(x, y);

				triangleNodes.xyTors(x, y, r, s);

				triangleNodes.rsToab(r, s, a, b);

				triangleNodes.evaluateGradSimplex(a, b, 1, 2, dpdr, dpds);

				Assert::That(dpdr(0) -  2.44948974278318, IsLessThan(eps));
				Assert::That(dpdr(1) -  2.44948974278318, IsLessThan(eps));
				Assert::That(dpdr(2) -  2.44948974278318, IsLessThan(eps));
				Assert::That(dpdr(3) -  2.44948974278318, IsLessThan(eps));
				Assert::That(dpdr(4) - -1.74516635192836, IsLessThan(eps));
				Assert::That(dpdr(5) - -1.63299316185545, IsLessThan(eps));
				Assert::That(dpdr(6) - -1.74516635192836, IsLessThan(eps));
				Assert::That(dpdr(7) -  8.11383968316462, IsLessThan(eps));
				Assert::That(dpdr(8) -  8.11383968316462, IsLessThan(eps));
				Assert::That(dpdr(9) -  24.49489742783178, IsLessThan(eps));

				Assert::That(dpds(0) -  15.921683328090657, IsLessThan(eps));
				Assert::That(dpds(1) -  7.797415561453585, IsLessThan(eps));
				Assert::That(dpds(2) - -5.347925818670407, IsLessThan(eps));
				Assert::That(dpds(3) - -13.472193585307480, IsLessThan(eps));
				Assert::That(dpds(4) - -0.525635522272996, IsLessThan(eps));
				Assert::That(dpds(5) - -0.816496580927726, IsLessThan(eps));
				Assert::That(dpds(6) - -1.219530829655363, IsLessThan(eps));
				Assert::That(dpds(7) - -2.168803194788500, IsLessThan(eps));
				Assert::That(dpds(8) -  10.282642877953123, IsLessThan(eps));
				Assert::That(dpds(9) -  12.247448713915887, IsLessThan(eps));
			}
		};
   } // namespace Nodes1DProvisionerTests
} // namespace blitzdg
