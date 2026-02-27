#include <gtest/gtest.h>
#include "test_common.h"
#include <atopo/concepts.h>

using namespace atopo::concepts;
using namespace atopo::test;

TEST(ConceptsUnit, CellConcept) {
    EXPECT_TRUE(Cell<MockVertex>);
    EXPECT_TRUE(Cell<MockEdge>);
    EXPECT_TRUE(Cell<MockFace>);
    EXPECT_TRUE(Cell<MockCell>);
}

TEST(ConceptsUnit, ComplexConcept) {
    EXPECT_TRUE(Complex<MockMesh>);
}

TEST(ConceptsUnit, GenericBuild) {
    auto circle = make_circle();
    // Verify we can build CellComplex from MockMesh
    auto complex = atopo::CellComplex::build(circle);
    
    EXPECT_EQ(complex.getNumberOfCells(0), 3);
    EXPECT_EQ(complex.getNumberOfCells(1), 3);
    EXPECT_EQ(complex.getNumberOfCells(2), 0);
}
