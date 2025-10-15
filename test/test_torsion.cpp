#include <gtest/gtest.h>
#include "atopo.h"

TEST(HomologyTorsionTest, RealProjectivePlaneCanonical) {
    atopo::CellComplex complex;

    // Build the canonical cell complex for RP^2 directly.
    // One cell in each dimension 0, 1, 2.
    complex.setCellCount(0, 1); // v
    complex.setCellCount(1, 1); // e
    complex.setCellCount(2, 1); // f

    // Define boundary map d1: C1 -> C0.
    // d1(e) = v - v = 0. Matrix is [0].
    atopo::IncidenceMatrix<atopo::IncidenceCoeff> d1(1, 1);
    complex.setBoundaryOperator(1, std::move(d1));

    // Define boundary map d2: C2 -> C1.
    // d2(f) = 2e. Matrix is [2].
    atopo::IncidenceMatrix<atopo::IncidenceCoeff> d2(1, 1);
    d2.insert(0, 0) = 2;
    complex.setBoundaryOperator(2, std::move(d2));

    // --- Compute Homology ---
    auto homology_groups = complex.computeHomology();
    
    // Expected for RP^2: H0=Z, H1=Z_2, H2=0
    // Betti numbers: β0=1, β1=0, β2=0
    // Torsion: T1={2}
    
    ASSERT_TRUE(homology_groups.count(0));
    EXPECT_EQ(homology_groups.at(0).rank, 1);
    EXPECT_TRUE(homology_groups.at(0).torsion_coeffs.empty());

    ASSERT_TRUE(homology_groups.count(1));
    EXPECT_EQ(homology_groups.at(1).rank, 0) << "Betti-1 should be 0 for RP^2";
    
    auto torsion1 = homology_groups.at(1).torsion_coeffs;
    ASSERT_EQ(torsion1.size(), 1) << "There should be one torsion coefficient in H1";
    EXPECT_EQ(torsion1[0], 2) << "The torsion coefficient of H1 should be 2";
    
    // H2 should be trivial and might not even be in the map if all its values are zero.
    if (homology_groups.count(2)) {
        EXPECT_EQ(homology_groups.at(2).rank, 0);
        EXPECT_TRUE(homology_groups.at(2).torsion_coeffs.empty());
    }
}
