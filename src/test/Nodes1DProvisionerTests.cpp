// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

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
        const float epsf = 5.7e-6;

        firstIndex ii;
        secondIndex jj;

        Nodes1DProvisioner * nodes1DProvisioner = nullptr;

        Describe(Nodes1DProvisioner_Object) {
            void SetUp() {
                const index_type NOrder = 3;
                const index_type NumElements = 5;
                const real_type xmin = -1.0;
                const real_type xmax = 1.0;

                nodes1DProvisioner = new Nodes1DProvisioner(NOrder, NumElements, xmin, xmax);
            }

            It(Should_Build_3rd_Order_Vandermonde_Matrix) {
                cout << "Should_Build_3rd_Order_Vandermonde_Matrix" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();
                nodes1D.buildVandermondeMatrix();
                const real_matrix_type & V = nodes1D.get_V();

                real_matrix_type expectedV(4,4);
                expectedV = 0.70711,-1.22474, 1.58114,-1.87083,
                            0.70711,-0.54772,-0.31623, 0.83666,
                            0.70711, 0.54772,-0.31623,-0.83666,
                            0.70711, 1.22474, 1.58114, 1.87083;

                real_matrix_type res(4,4);
                res  = V - expectedV;
                Assert::That(sqrt(sum(res(ii)*res(ii))), IsLessThan(epsf));
            }

            It(Should_Build_4th_Order_GradVandermonde_Matrix) {
                cout << "Should_Build_4th_Order_GradVandermonde_Matrix" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();

                real_matrix_type DVr(5,5);
                nodes1D.computeGradVandermonde(DVr);

                real_matrix_type expectedDVr(5,5);
                expectedDVr = 0.00000,1.22474,-4.74342,11.22497,-21.21320,
                            0.00000,1.22474,-3.10530, 3.20713, -0.00000,
                            0.00000,1.22474,-0.00000,-2.80624,  0.00000,
                            0.00000,1.22474, 3.10530, 3.20713,  0.00000,
                            0.00000,1.22474 ,4.74342,11.22497, 21.21320;

                real_matrix_type res(4,4);
                res = DVr - expectedDVr;
                Assert::That(sqrt(sum(res(ii)*res(ii))), IsLessThan(epsf));
            }

            It(Should_Build_3rd_Order_Differentiation_Matrix) {
                cout << "Should_Build_3rd_Order_Differentiation_Matrix" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
                nodes1D.buildNodes();

                const real_matrix_type & Dr = nodes1D.get_Dr();

                real_matrix_type expectedDr(4,4);
                expectedDr = -3.0000e+00, 4.0451e+00,-1.5451e+00, 5.0000e-01,
                            -8.0902e-01,-4.0540e-16, 1.1180e+00,-3.0902e-01,
                            3.0902e-01,-1.1180e+00, 6.2804e-16, 8.0902e-01,
                            -5.0000e-01, 1.5451e+00,-4.0451e+00, 3.0000e+00; 


                real_matrix_type res(4,4);
                res = Dr - expectedDr;
                Assert::That(sqrt(sum(res(ii)*res(ii))), IsLessThan(epsf));
            }

            It(Should_Build_A_1D_X_Grid) {
                cout << "Should_Build_A_1D_X_Grid" << endl;

                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();

                const real_matrix_type & x = nodes1D.get_xGrid();

                real_matrix_type expectedx(4,5);
                expectedx = -1.000000,-0.600000,-0.200000,0.200000,0.600000,
                            -0.889443,-0.489443,-0.089443,0.310557,0.710557,
                            -0.710557,-0.310557, 0.089443,0.489443,0.889443,
                            -0.600000,-0.200000, 0.200000,0.600000,1.000000;

                real_matrix_type res(4,5);
                res = x - expectedx;
                Assert::That(sqrt(sum(res(ii)*res(ii))), IsLessThan(epsf));
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

                nodes1D.buildDr();
                nodes1D.computeJacobian();

                const real_matrix_type & J = nodes1D.get_J();
                const real_matrix_type & rx = nodes1D.get_rx();
                const real_matrix_type & Fscale = nodes1D.get_Fscale();

                real_matrix_type expectedJ(4,5), expectedrx(4,5), expectedFscale(2,5);
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

                real_matrix_type resJ(4,5), resrx(4,5), resFscale(2,5);
                resJ = J - expectedJ;
                resrx = rx - expectedrx;
                resFscale = Fscale - expectedFscale;

                Assert::That(sqrt(sum(resJ(ii)*resJ(ii))), IsLessThan(epsf));
                Assert::That(sqrt(sum(resrx(ii)*resrx(ii))), IsLessThan(epsf));
                Assert::That(sqrt(sum(resFscale(ii)*resFscale(ii))), IsLessThan(epsf));
            }

            It(Should_Build_1D_Lift_Operator) {
                cout << "Should_Build_1D_Lift_Operator" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();
                real_matrix_type Lift = nodes1D.get_Lift();

                real_matrix_type expectedLift(4,2);
                expectedLift =  8.00000,-2.00000,
                               -0.89443, 0.89443,
                                0.89443,-0.89443,
                               -2.00000, 8.00000;


                real_matrix_type resLift(4,2);

                resLift = Lift - expectedLift;
                Assert::That(sqrt(sum(resLift*resLift)), IsLessThan(epsf));
            }

            It(Should_Build_1D_Connectivity_Matrices) {
                cout << "Should_Build_1D_Connectivity_Matrices" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

                nodes1D.buildNodes();
                const index_matrix_type & EToE = nodes1D.get_EToE();
                const index_matrix_type & EToF = nodes1D.get_EToF();

                index_matrix_type expectedEToE(5,2);
                expectedEToE = 0,1,
                            0,2,
                            1,3,
                            2,4,
                            3,4;

                index_matrix_type expectedEToF(5,2);
                expectedEToF = 0,0,
                            1,0,
                            1,0,
                            1,0,
                            1,1;

                index_matrix_type resEToE(5,2), resEToF(5,2);

                resEToE = EToE - expectedEToE;
                resEToF = EToF - expectedEToF;

                Assert::That(sqrt(sum(resEToE*resEToE)), Equals(0));
                Assert::That(sqrt(sum(resEToF*resEToF)), Equals(0));
            }

            It(Should_Build_Face_Mask) {
                cout << "Should_Build_Face_Mask" << endl;
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
                nodes1D.buildNodes();

                index_vector_type Fmask = nodes1D.get_Fmask();
                real_matrix_type Fx = nodes1D.get_Fx();

                Assert::That(Fmask(0), Equals(0));
                Assert::That(Fmask(1), Equals(3));

                real_matrix_type expectedFx(2, 5);
                expectedFx = -1,-0.6,-0.2,0.2,0.6,
                            -0.6,-0.2,0.2,0.6,1;

                real_matrix_type resFx;
                resFx = Fx - expectedFx;

                Assert::That(sqrt(sum(resFx*resFx)), Equals(0));
            }

            It(Should_Build_Volume_Maps) {
                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
                nodes1D.buildNodes();

                index_vector_type vmapM = nodes1D.get_vmapM();
                index_vector_type vmapP = nodes1D.get_vmapP();

                index_vector_type expectedVmapM(10);
                index_vector_type expectedVmapP(10);

                expectedVmapM = 0,3,4,7,8,11,12,15,16,19;
                expectedVmapP = 0,4,3,8,7,12,11,16,15,19;

                index_vector_type resVmapM(10);
                index_vector_type resVmapP(10);

                resVmapM = vmapM - expectedVmapM;
                resVmapP = vmapP - expectedVmapP;

                Assert::That(sqrt(sum(resVmapM*resVmapM)), Equals(0));
                Assert::That(sqrt(sum(resVmapP*resVmapP)), Equals(0));
            }

            It(Should_Build_Normals) {

                Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
                nodes1D.buildNodes();

                const real_matrix_type & nx = nodes1D.get_nx();

                real_matrix_type expectednx(2,5);
                expectednx = -1,-1,-1,-1,-1,
                              1, 1, 1, 1, 1;

                real_matrix_type resnx(2,5);

                resnx = nx - expectednx;

                Assert::That(sqrt(sum(resnx*resnx)), Equals(0));
            }
        };
   } // namespace Nodes1DProvisionerTests
} // namespace blitzdg
