// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "LinAlgHelpers.hpp"
#include "Nodes1DProvisioner.hpp"
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
    namespace Nodes1DProvisionerTests {
        using namespace igloo;

        firstIndex ii;
        secondIndex jj;

        Nodes1DProvisioner * nodes1DProvisioner = nullptr;

        Describe(Nodes1DProvisioner_Object) {
			const real_type eps = 2.*numeric_limits<double>::epsilon();
			const float epsf = 5.8e-5;
			const index_type NOrder = 3;
			const index_type NumElements = 5;
			const int NumFaces = 2;

            void SetUp() {
                const real_type xmin = -1.0;
                const real_type xmax = 1.0;

                nodes1DProvisioner = new Nodes1DProvisioner(NOrder, NumElements, xmin, xmax);
            }

			void TearDown() {
				delete nodes1DProvisioner;
			}

			It(Should_Build_3rd_Order_Vandermonde_Matrix) {
                cout << "Should_Build_3rd_Order_Vandermonde_Matrix" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();
                const real_matrix_type & V = nodes1D.get_V();

                real_matrix_type expectedV(NOrder+1,NOrder+1);
                expectedV = 0.70711,-1.22474, 1.58114,-1.87083,
                            0.70711,-0.54772,-0.31623, 0.83666,
                            0.70711, 0.54772,-0.31623,-0.83666,
                            0.70711, 1.22474, 1.58114, 1.87083;

                real_matrix_type res(NOrder+1,NOrder+1);
                res  = V - expectedV;
                Assert::That(normFro(res), IsLessThan(epsf));
            }

            It(Should_Build_3rd_Order_Differentiation_Matrix) {
                cout << "Should_Build_3rd_Order_Differentiation_Matrix" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
                nodes1D.buildNodes();

                const real_matrix_type & Dr = nodes1D.get_Dr();

                real_matrix_type expectedDr(NOrder+1,NOrder+1);
                expectedDr = -3.0000e+00, 4.0451e+00,-1.5451e+00, 5.0000e-01,
                            -8.0902e-01,-4.0540e-16, 1.1180e+00,-3.0902e-01,
                            3.0902e-01,-1.1180e+00, 6.2804e-16, 8.0902e-01,
                            -5.0000e-01, 1.5451e+00,-4.0451e+00, 3.0000e+00; 

                real_matrix_type res(NOrder+1,NOrder+1);
                res = Dr - expectedDr;
                Assert::That(normFro(res), IsLessThan(epsf));
            }

            It(Should_Build_A_1D_X_Grid) {
                cout << "Should_Build_A_1D_X_Grid" << endl;

                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();

                const real_matrix_type & x = nodes1D.get_xGrid();

                real_matrix_type expectedx(NOrder+1, NumElements);
                expectedx = -1.000000,-0.600000,-0.200000,0.200000,0.600000,
                            -0.889443,-0.489443,-0.089443,0.310557,0.710557,
                            -0.710557,-0.310557, 0.089443,0.489443,0.889443,
                            -0.600000,-0.200000, 0.200000,0.600000,1.000000;

                real_matrix_type res(NOrder+1, NumElements);
                res = x - expectedx;
                Assert::That(normFro(res), IsLessThan(epsf));
            }

            It(Should_Build_Element_To_Vertex_Connectivity) {
                cout << "Should_Build_Element_To_Vertex_Connectivity" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();

                const index_matrix_type & EToV = nodes1D.get_EToV();

                Assert::That(EToV(0,0), Equals(0)); Assert::That(EToV(0,1), Equals(1));
                Assert::That(EToV(1,0), Equals(1)); Assert::That(EToV(1,1), Equals(2));
                Assert::That(EToV(2,0), Equals(2)); Assert::That(EToV(2,1), Equals(3));
                Assert::That(EToV(3,0), Equals(3)); Assert::That(EToV(3,1), Equals(4));
                Assert::That(EToV(4,0), Equals(4)); Assert::That(EToV(4,1), Equals(5));
            }

            It(Should_Compute_Jacobian) {
                cout << "Should_Compute_Jacobian" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();

                cout << "Nodes built.\n";

                nodes1D.computeJacobian();

                cout << "Jacobian computed.\n";


                const real_matrix_type & J = nodes1D.get_J();
                const real_matrix_type & rx = nodes1D.get_rx();
                const real_matrix_type & Fscale = nodes1D.get_Fscale();

                real_matrix_type expectedJ(NOrder+1, NumElements), expectedrx(NOrder+1, NumElements), expectedFscale(2, NumElements);
                expectedJ = 0.20000,0.20000,0.20000,0.20000,0.20000,
                            0.20000,0.20000,0.20000,0.20000,0.20000,
                            0.20000,0.20000,0.20000,0.20000,0.20000,
                            0.20000,0.20000,0.20000,0.20000,0.20000;

                expectedrx =5,5,5,5,5,
                            5,5,5,5,5,
                            5,5,5,5,5,
                            5,5,5,5,5;
                
                expectedFscale = 5,5,5,5,5,
                                 5,5,5,5,5;

                real_matrix_type resJ(NOrder+1, NumElements), resrx(NOrder+1, NumElements), resFscale(2, NumElements);
                resJ = J - expectedJ;
                resrx = rx - expectedrx;
                resFscale = Fscale - expectedFscale;

                Assert::That(normFro(resJ), IsLessThan(epsf));
                Assert::That(normFro(resrx), IsLessThan(epsf));
                Assert::That(normFro(resFscale), IsLessThan(epsf));
            }

            It(Should_Build_1D_Lift_Operator) {
                cout << "Should_Build_1D_Lift_Operator" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();
                real_matrix_type Lift = nodes1D.get_Lift();

                real_matrix_type expectedLift(NOrder+1,NumFaces);
                expectedLift =  8.00000,-2.00000,
                               -0.89443, 0.89443,
                                0.89443,-0.89443,
                               -2.00000, 8.00000;


                real_matrix_type resLift(NOrder+1,2);

                resLift = Lift - expectedLift;
                Assert::That(normFro(resLift), IsLessThan(epsf));
            }

            It(Should_Build_1D_Connectivity_Matrices) {
                cout << "Should_Build_1D_Connectivity_Matrices" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();
                const index_matrix_type & EToE = nodes1D.get_EToE();
                const index_matrix_type & EToF = nodes1D.get_EToF();

                index_matrix_type expectedEToE(NumElements,NumFaces);
                expectedEToE = 0,1,
                            0,2,
                            1,3,
                            2,4,
                            3,4;

                index_matrix_type expectedEToF(NumElements,NumFaces);
                expectedEToF = 0,0,
                            1,0,
                            1,0,
                            1,0,
                            1,1;

                index_matrix_type resEToE(NumElements,NumFaces), resEToF(NumElements,NumFaces);

                resEToE = EToE - expectedEToE;
                resEToF = EToF - expectedEToF;

                Assert::That(normFro(resEToE), Equals(0));
                Assert::That(normFro(resEToF), Equals(0));
            }

            It(Should_Build_Face_Mask) {
                cout << "Should_Build_Face_Mask" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
                nodes1D.buildNodes();

                const index_vector_type & Fmask = nodes1D.get_Fmask();
                const real_matrix_type & Fx = nodes1D.get_Fx();

                Assert::That(Fmask(0), Equals(0));
                Assert::That(Fmask(1), Equals(3));

                real_matrix_type expectedFx(NumFaces, NumElements);
                expectedFx = -1,-0.6,-0.2,0.2,0.6,
                            -0.6,-0.2,0.2,0.6,1;

				cout << "Fx: " << Fx << endl;
				cout << "expectedFx: " << expectedFx << endl;

                real_matrix_type resFx(NumFaces, NumElements);
                resFx = Fx - expectedFx;

                Assert::That(normFro(resFx), IsLessThan(eps));
            }

            It(Should_Build_Volume_Maps) {
                cout << "Should_Build_Volume_Maps" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
                nodes1D.buildNodes();

                index_vector_type vmapM = nodes1D.get_vmapM();
                index_vector_type vmapP = nodes1D.get_vmapP();

                index_vector_type expectedVmapM(NumFaces*NumElements);
                index_vector_type expectedVmapP(NumFaces*NumElements);

                expectedVmapM = 0,3,4,7,8,11,12,15,16,19;
                expectedVmapP = 0,4,3,8,7,12,11,16,15,19;

                index_vector_type resVmapM(NumFaces*NumElements);
                index_vector_type resVmapP(NumFaces*NumElements);

                resVmapM = vmapM - expectedVmapM;
                resVmapP = vmapP - expectedVmapP;

                Assert::That(norm2(resVmapM), Equals(0));
                Assert::That(norm2(resVmapP), Equals(0));
            }

            It(Should_Build_Normals) {
                cout << "Should_Build_Normals" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
                nodes1D.buildNodes();

                const real_matrix_type & nx = nodes1D.get_nx();

                real_matrix_type expectednx(NumFaces,NumElements);
                expectednx = -1,-1,-1,-1,-1,
                              1, 1, 1, 1, 1;

                real_matrix_type resnx(NumFaces,NumElements);

                resnx = nx - expectednx;

                Assert::That(normFro(resnx), Equals(0));
            }
		};
   } // namespace Nodes1DProvisionerTests
} // namespace blitzdg
