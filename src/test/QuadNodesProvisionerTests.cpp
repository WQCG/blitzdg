// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "LinAlgHelpers.hpp"
#include "QuadNodesProvisioner.hpp"
#include "MeshManager.hpp"
#include "Types.hpp"
#include "PathResolver.hpp"
#include <whereami.h>
#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <iostream>
#include <limits>
#include <cmath>
#include <memory>
#include <vector>

using boost::algorithm::find_all;
using boost::algorithm::join;
using boost::algorithm::replace_last;
using boost::algorithm::trim_right;
using boost::iterator_range;
using blitz::firstIndex;
using blitz::secondIndex;
using std::string;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::abs;
using std::unique_ptr;
using std::shared_ptr;
using std::vector;

namespace blitzdg {
    namespace QuadNodesProvisionerTests {
        using namespace igloo;

        firstIndex ii;
        secondIndex jj;

		unique_ptr<QuadNodesProvisioner> quadNodesProvisioner = nullptr;
		shared_ptr<MeshManager> meshManager =  nullptr;

        Describe(QuadNodesProvisioner_Object) {
			const real_type eps = 50*numeric_limits<double>::epsilon();
			const float epsf = 5.8e-5;
			const index_type NOrder = 3;
			const index_type NumFaces = 2;
            string PathDelimeter = "/";

			using find_vector_type = vector<iterator_range<string::iterator>>;

            void SetUp() {
				PathResolver resolver;
				meshManager = shared_ptr<MeshManager>(new MeshManager());

				string root = resolver.get_RootPath();
				string inputPath = resolver.joinPaths(root, "input");
				string meshPath = resolver.joinPaths(inputPath, "coarse_box_quads.msh");

				meshManager->readMesh(meshPath);
				quadNodesProvisioner = unique_ptr<QuadNodesProvisioner>(new QuadNodesProvisioner(NOrder, *meshManager)); 
			}

			It(Should_Build_2d_rs_grid) {
				QuadNodesProvisioner & quadNodes = *quadNodesProvisioner;
				quadNodes.buildNodes();

				const DGContext2D& ctx = quadNodes.get_DGContext();
				const real_vector_type& r = ctx.r();
				const real_vector_type& s = ctx.s();

				index_type numLocalPoints = (NOrder+1)*(NOrder+1);
				real_vector_type rExpected(numLocalPoints), sExpected(numLocalPoints);

				rExpected = -1.00e+00,-1.00e+00,-1.00e+00,-1.00e+00,-4.47e-01,-4.47e-01,-4.47e-01,-4.47e-01,4.47e-01,4.47e-01,4.47e-01,4.47e-01,1.00e+00,1.00e+00,1.00e+00,1.00e+00;
				sExpected = -1.00e+00,-4.47e-01,4.47e-01,1.00e+00,-1.00e+00,-4.47e-01,4.47e-01,1.00e+00,-1.00e+00,-4.47e-01,4.47e-01,1.00e+00,-1.00e+00,-4.47e-01,4.47e-01,1.00e+00;

				rExpected -= r;
				sExpected -= s;

				Assert::That(normInf(rExpected), IsLessThan(1.e-3));
				Assert::That(normInf(sExpected), IsLessThan(1.e-3));
			}

			It(Should_Compute_GradVandermonde_Properly) {
                cout << "Should_Compute_GradVandermonde_Properly" << endl;
                QuadNodesProvisioner & quadNodes = *quadNodesProvisioner;

				quadNodes.buildNodes();

				const DGContext2D& ctx = quadNodes.get_DGContext();
				const real_vector_type& r = ctx.r();
				const real_vector_type& s = ctx.s();

				index_type numLocalPoints = (NOrder+1)*(NOrder+1);

				real_matrix_type V2Dr(numLocalPoints, numLocalPoints),
					V2Ds(numLocalPoints, numLocalPoints);
				V2Dr = 0.0; V2Ds = 0.0;

				quadNodes.computeGradVandermondeMatrix(NOrder, r, s, V2Dr, V2Ds);

				real_matrix_type V2DrExpected(numLocalPoints, numLocalPoints),
					V2DsExpected(numLocalPoints, numLocalPoints);

				V2DrExpected = 0.00000,0.86603,-3.35410,7.93725,-0.00000,-1.50000,5.80948,-13.74773,0.00000,1.93649,-7.50000,17.74824,-0.00000,-2.29129,8.87412,-21.00000
					,0.00000,0.86603,-3.35410,7.93725,-0.00000,-0.67082,2.59808,-6.14817,-0.00000,-0.38730,1.50000,-3.54965,0.00000,1.02470,-3.96863,9.39149
					,0.00000,0.86603,-3.35410,7.93725,0.00000,0.67082,-2.59808,6.14817,-0.00000,-0.38730,1.50000,-3.54965,-0.00000,-1.02470,3.96863,-9.39149
					,0.00000,0.86603,-3.35410,7.93725,0.00000,1.50000,-5.80948,13.74773,0.00000,1.93649,-7.50000,17.74824,0.00000,2.29129,-8.87412,21.00000
					,0.00000,0.86603,-1.50000,0.00000,-0.00000,-1.50000,2.59808,-0.00000,0.00000,1.93649,-3.35410,0.00000,-0.00000,-2.29129,3.96863,-0.00000
					,0.00000,0.86603,-1.50000,0.00000,-0.00000,-0.67082,1.16190,-0.00000,-0.00000,-0.38730,0.67082,-0.00000,0.00000,1.02470,-1.77482,0.00000
					,0.00000,0.86603,-1.50000,0.00000,0.00000,0.67082,-1.16190,0.00000,-0.00000,-0.38730,0.67082,-0.00000,-0.00000,-1.02470,1.77482,-0.00000
					,0.00000,0.86603,-1.50000,0.00000,0.00000,1.50000,-2.59808,0.00000,0.00000,1.93649,-3.35410,0.00000,0.00000,2.29129,-3.96863,0.00000
					,0.00000,0.86603,1.50000,0.00000,-0.00000,-1.50000,-2.59808,-0.00000,0.00000,1.93649,3.35410,0.00000,-0.00000,-2.29129,-3.96863,-0.00000
					,0.00000,0.86603,1.50000,0.00000,-0.00000,-0.67082,-1.16190,-0.00000,-0.00000,-0.38730,-0.67082,-0.00000,0.00000,1.02470,1.77482,0.00000
					,0.00000,0.86603,1.50000,0.00000,0.00000,0.67082,1.16190,0.00000,-0.00000,-0.38730,-0.67082,-0.00000,-0.00000,-1.02470,-1.77482,-0.00000
					,0.00000,0.86603,1.50000,0.00000,0.00000,1.50000,2.59808,0.00000,0.00000,1.93649,3.35410,0.00000,0.00000,2.29129,3.96863,0.00000
					,0.00000,0.86603,3.35410,7.93725,-0.00000,-1.50000,-5.80948,-13.74773,0.00000,1.93649,7.50000,17.74824,-0.00000,-2.29129,-8.87412,-21.00000
					,0.00000,0.86603,3.35410,7.93725,-0.00000,-0.67082,-2.59808,-6.14817,-0.00000,-0.38730,-1.50000,-3.54965,0.00000,1.02470,3.96863,9.39149
					,0.00000,0.86603,3.35410,7.93725,0.00000,0.67082,2.59808,6.14817,-0.00000,-0.38730,-1.50000,-3.54965,-0.00000,-1.02470,-3.96863,-9.39149
					,0.00000,0.86603,3.35410,7.93725,0.00000,1.50000,5.80948,13.74773,0.00000,1.93649,7.50000,17.74824,0.00000,2.29129,8.87412,21.00000;

				V2DsExpected = 0.00000,-0.00000,0.00000,-0.00000,0.86603,-1.50000,1.93649,-2.29129,-3.35410,5.80948,-7.50000,8.87412,7.93725,-13.74773,17.74824,-21.00000
					,0.00000,-0.00000,0.00000,-0.00000,0.86603,-1.50000,1.93649,-2.29129,-1.50000,2.59808,-3.35410,3.96863,0.00000,-0.00000,0.00000,-0.00000
					,0.00000,-0.00000,0.00000,-0.00000,0.86603,-1.50000,1.93649,-2.29129,1.50000,-2.59808,3.35410,-3.96863,0.00000,-0.00000,0.00000,-0.00000
					,0.00000,-0.00000,0.00000,-0.00000,0.86603,-1.50000,1.93649,-2.29129,3.35410,-5.80948,7.50000,-8.87412,7.93725,-13.74773,17.74824,-21.00000
					,0.00000,-0.00000,-0.00000,0.00000,0.86603,-0.67082,-0.38730,1.02470,-3.35410,2.59808,1.50000,-3.96863,7.93725,-6.14817,-3.54965,9.39149
					,0.00000,-0.00000,-0.00000,0.00000,0.86603,-0.67082,-0.38730,1.02470,-1.50000,1.16190,0.67082,-1.77482,0.00000,-0.00000,-0.00000,0.00000
					,0.00000,-0.00000,-0.00000,0.00000,0.86603,-0.67082,-0.38730,1.02470,1.50000,-1.16190,-0.67082,1.77482,0.00000,-0.00000,-0.00000,0.00000
					,0.00000,-0.00000,-0.00000,0.00000,0.86603,-0.67082,-0.38730,1.02470,3.35410,-2.59808,-1.50000,3.96863,7.93725,-6.14817,-3.54965,9.39149
					,0.00000,0.00000,-0.00000,-0.00000,0.86603,0.67082,-0.38730,-1.02470,-3.35410,-2.59808,1.50000,3.96863,7.93725,6.14817,-3.54965,-9.39149
					,0.00000,0.00000,-0.00000,-0.00000,0.86603,0.67082,-0.38730,-1.02470,-1.50000,-1.16190,0.67082,1.77482,0.00000,0.00000,-0.00000,-0.00000
					,0.00000,0.00000,-0.00000,-0.00000,0.86603,0.67082,-0.38730,-1.02470,1.50000,1.16190,-0.67082,-1.77482,0.00000,0.00000,-0.00000,-0.00000
					,0.00000,0.00000,-0.00000,-0.00000,0.86603,0.67082,-0.38730,-1.02470,3.35410,2.59808,-1.50000,-3.96863,7.93725,6.14817,-3.54965,-9.39149
					,0.00000,0.00000,0.00000,0.00000,0.86603,1.50000,1.93649,2.29129,-3.35410,-5.80948,-7.50000,-8.87412,7.93725,13.74773,17.74824,21.00000
					,0.00000,0.00000,0.00000,0.00000,0.86603,1.50000,1.93649,2.29129,-1.50000,-2.59808,-3.35410,-3.96863,0.00000,0.00000,0.00000,0.00000
					,0.00000,0.00000,0.00000,0.00000,0.86603,1.50000,1.93649,2.29129,1.50000,2.59808,3.35410,3.96863,0.00000,0.00000,0.00000,0.00000
					,0.00000,0.00000,0.00000,0.00000,0.86603,1.50000,1.93649,2.29129,3.35410,5.80948,7.50000,8.87412,7.93725,13.74773,17.74824,21.00000;				

				V2DrExpected -= V2Dr;
				V2DsExpected -= V2Ds;

				Assert::That(normMax(V2DrExpected), IsLessThan(1.e-4));
				Assert::That(normMax(V2DsExpected), IsLessThan(1.e-4));
            }

			It(Should_Build_Face_Masks) {
				cout << "Should_Build_Face_Masks" << endl;
                QuadNodesProvisioner & quadNodes = *quadNodesProvisioner;

				const index_matrix_type& Fm = quadNodes.get_Fmask();

				const index_type numFacePoints = (NOrder + 1);
				index_matrix_type FmExpected(numFacePoints, 4);

				FmExpected = 0,12,3,0
					,4,13,7,1
					,8,14,11,2
					,12,15,15,3;

				FmExpected -= Fm;

				Assert::That(normMax(FmExpected), Equals(0.0));
			}

			It(Should_Build_Lift) {
				cout << "Should_Build_Lift" << endl;
                QuadNodesProvisioner & quadNodes = *quadNodesProvisioner;
				const real_matrix_type& lift = quadNodes.get_Lift();
				
				real_matrix_type liftExpected(4*(NOrder+1), (NOrder+1)*(NOrder+1));

				liftExpected = 8.00e+00,2.22e-15,-3.33e-16,4.44e-16,-2.00e+00,-6.66e-16,-3.61e-16,-5.55e-16,-2.00e+00,-6.66e-16,-3.61e-16,-5.55e-16,8.00e+00,2.22e-15,-3.33e-16,4.44e-16
					,-8.94e-01,-1.39e-17,4.16e-17,-5.55e-17,-7.63e-17,-2.00e+00,-1.25e-16,-5.55e-17,8.94e-01,5.55e-17,1.53e-16,0.00e+00,2.78e-16,8.00e+00,-3.33e-16,0.00e+00
					,8.94e-01,5.55e-17,1.53e-16,0.00e+00,-1.04e-16,-1.25e-16,-2.00e+00,-1.11e-16,-8.94e-01,-1.39e-17,4.16e-17,-5.55e-17,3.05e-16,-7.22e-16,8.00e+00,4.44e-16
					,-2.00e+00,-6.66e-16,-3.61e-16,-5.55e-16,-5.55e-16,0.00e+00,-6.66e-16,-2.00e+00,8.00e+00,2.22e-15,-3.33e-16,4.44e-16,2.22e-16,-8.88e-16,2.22e-15,8.00e+00
					,2.78e-17,8.00e+00,1.11e-16,0.00e+00,8.94e-01,2.78e-16,2.08e-16,0.00e+00,-9.02e-17,-2.00e+00,-1.39e-17,-5.55e-17,-8.94e-01,-9.71e-17,-1.25e-16,-5.55e-17
					,-2.43e-17,-8.94e-01,-6.94e-18,6.94e-17,3.47e-17,8.94e-01,7.63e-17,2.78e-17,-2.08e-17,8.94e-01,1.04e-16,1.39e-17,-6.94e-18,-8.94e-01,2.78e-17,-1.39e-17
					,-2.08e-17,8.94e-01,1.04e-16,1.39e-17,4.86e-17,2.08e-17,8.94e-01,5.55e-17,-2.43e-17,-8.94e-01,-6.94e-18,6.94e-17,-5.20e-17,1.39e-17,-8.94e-01,-2.78e-17
					,-9.02e-17,-2.00e+00,-1.39e-17,-5.55e-17,2.78e-17,1.67e-16,2.78e-16,8.94e-01,2.78e-17,8.00e+00,1.11e-16,0.00e+00,-8.33e-17,-1.11e-16,-5.55e-17,-8.94e-01
					,3.05e-16,-2.22e-16,8.00e+00,2.22e-16,-8.94e-01,-9.71e-17,-1.25e-16,-5.55e-17,-1.04e-16,9.71e-17,-2.00e+00,-1.11e-16,8.94e-01,2.78e-16,2.08e-16,0.00e+00
					,3.12e-17,1.39e-17,-8.94e-01,-4.16e-17,-6.94e-18,-8.94e-01,2.78e-17,-1.39e-17,3.82e-17,4.16e-17,8.94e-01,0.00e+00,3.47e-17,8.94e-01,7.63e-17,2.78e-17
					,3.82e-17,4.16e-17,8.94e-01,0.00e+00,-5.20e-17,1.39e-17,-8.94e-01,-2.78e-17,3.12e-17,1.39e-17,-8.94e-01,-4.16e-17,4.86e-17,2.08e-17,8.94e-01,5.55e-17
					,-1.04e-16,9.71e-17,-2.00e+00,-1.11e-16,-8.33e-17,-1.11e-16,-5.55e-17,-8.94e-01,3.05e-16,-2.22e-16,8.00e+00,2.22e-16,2.78e-17,1.67e-16,2.78e-16,8.94e-01
					,2.22e-16,-8.88e-16,1.33e-15,8.00e+00,8.00e+00,2.22e-15,-3.33e-16,4.44e-16,-5.55e-16,0.00e+00,-6.66e-16,-2.00e+00,-2.00e+00,-6.66e-16,-3.61e-16,-5.55e-16
					,-5.55e-17,5.55e-17,0.00e+00,-8.94e-01,2.78e-16,8.00e+00,-3.33e-16,0.00e+00,2.78e-17,1.67e-16,5.55e-17,8.94e-01,-7.63e-17,-2.00e+00,-1.25e-16,-5.55e-17
					,2.78e-17,1.67e-16,5.55e-17,8.94e-01,3.05e-16,-7.22e-16,8.00e+00,4.44e-16,-5.55e-17,5.55e-17,0.00e+00,-8.94e-01,-1.04e-16,-1.25e-16,-2.00e+00,-1.11e-16
					,-5.55e-16,0.00e+00,-6.66e-16,-2.00e+00,2.22e-16,-8.88e-16,2.22e-15,8.00e+00,2.22e-16,-8.88e-16,1.33e-15,8.00e+00,-5.55e-16,0.00e+00,-6.66e-16,-2.00e+00;
				
				liftExpected -= lift;
				Assert::That(normMax(liftExpected), IsLessThan(1.e-3));
			}

		};
   } // namespace QuadNodesProvisionerTests
} // namespace blitzdg
