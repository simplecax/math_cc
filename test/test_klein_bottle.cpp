#include <gtest/gtest.h>
#include "test_common.h"

TEST(HomologyIntegrationTest, KleinBottle) {
    // We can't use build(mesh) directly if we need custom boundary entries for a 1-vertex representation
    // But we can just use the manual building logic in the test.
    atopo::CellComplex complex;
    complex.setCellCount(0, 1); // v
    complex.setCellCount(1, 2); // a, b
    complex.setCellCount(2, 1); // f

    atopo::IncidenceMatrix<atopo::IncidenceCoeff> d1(1, 2); 
    complex.setBoundaryOperator(1, std::move(d1));

    atopo::IncidenceMatrix<atopo::IncidenceCoeff> d2(2, 1);
    d2.insert(1, 0) = 2; 
    complex.setBoundaryOperator(2, std::move(d2));

    auto homology_groups = complex.computeHomology();

    ASSERT_TRUE(homology_groups.count(0));
    EXPECT_EQ(homology_groups.at(0).rank, 1);
    ASSERT_TRUE(homology_groups.count(1));
    EXPECT_EQ(homology_groups.at(1).rank, 2); 
    ASSERT_TRUE(homology_groups.count(2));
    EXPECT_EQ(homology_groups.at(2).rank, 1);
}
