#include <gtest/gtest.h>
#include "test_common.h"

// This test verifies the homology calculation for the Klein Bottle,
// a non-orientable surface with both a non-zero Betti number and torsion in H1.
// This confirms the correctness of the Betti number formula over integers.
TEST(HomologyIntegrationTest, KleinBottle) {
    LegacyMesh mesh;
    // Standard square representation of a Klein Bottle.
    mesh.vertices = {{0, {}}, {1, {}}, {2, {}}, {3, {}}};
    // a=(0,1), b=(1,2), c=(2,3), d=(3,0)
    mesh.edges = { {0,0,1,{}}, {1,1,2,{}}, {2,2,3,{}}, {3,3,0,{}} };
    mesh.faces = { {0, {0,1,2,3}, {}} };
    
    // We need to build the complex manually for this non-standard identification.
    atopo::CellComplex complex;
    complex.setCellCount(0, 1); // v
    complex.setCellCount(1, 2); // a, b
    complex.setCellCount(2, 1); // f

    // d1(a) = v-v = 0, d1(b) = v-v = 0.
    atopo::IncidenceMatrix<atopo::IncidenceCoeff> d1(1, 2); // All zeros
    complex.setBoundaryOperator(1, std::move(d1));

    // d2(f) = a + b - a + b = 2b
    atopo::IncidenceMatrix<atopo::IncidenceCoeff> d2(2, 1);
    d2.insert(1, 0) = 2; // Coefficient of b is 2
    complex.setBoundaryOperator(2, std::move(d2));

    auto homology_groups = complex.computeHomology();

    // Expected for Klein Bottle: H0=Z, H1=Z x Z_2, H2=0
    
    ASSERT_TRUE(homology_groups.count(0));
    EXPECT_EQ(homology_groups.at(0).rank, 1);
    EXPECT_TRUE(homology_groups.at(0).torsion_coeffs.empty());

    ASSERT_TRUE(homology_groups.count(1));
    EXPECT_EQ(homology_groups.at(1).rank, 1); // β1 = 1
    
    auto torsion1 = homology_groups.at(1).torsion_coeffs;
    ASSERT_EQ(torsion1.size(), 1);
    EXPECT_EQ(torsion1[0], 2); // Torsion is Z_2
}
