#include <gtest/gtest.h>
#include "atopo.h"
#include <vector>
#include <algorithm>
#include <numeric>

// Test fixture for the RP^2 triangulation integration tests.
class RP2IntegrationTest : public ::testing::TestWithParam<int> {
protected:
    static void SetUpTestSuite() {
        LegacyMesh mesh;
        // The minimal triangulation of the Real Projective Plane (RP^2)
        mesh.vertices = {{0,{}}, {1,{}}, {2,{}}, {3,{}}, {4,{}}, {5,{}}};
        mesh.edges = {
            {0,0,1,{}}, {1,0,2,{}}, {2,0,3,{}}, {3,0,4,{}}, {4,0,5,{}},  // Edges from v0
            {5,1,2,{}}, {6,1,3,{}}, {7,1,4,{}}, {8,1,5,{}},              // Edges from v1
            {9,2,3,{}}, {10,2,4,{}}, {11,2,5,{}},                        // Edges from v2
            {12,3,4,{}}, {13,3,5,{}},                                    // Edges from v3
            {14,4,5,{}}                                                 // Edge from v4
        };
        // --- CORRECTED FACE DEFINITIONS ---
        // A standard, mathematically correct triangulation of RP^2
        mesh.faces = {
            {0, {0, 5, 1}},    // v0,v1,v2
            {1, {1, 9, 2}},    // v0,v2,v3
            {2, {2, 12, 3}},   // v0,v3,v4
            {3, {3, 14, 4}},   // v0,v4,v5
            {4, {4, 8, 0}},    // v0,v5,v1
            {5, {5, 10, 7}},   // v1,v2,v4
            {6, {6, 13, 8}},   // v1,v3,v5
            {7, {9, 13, 11}},  // v2,v3,v5
            {8, {11, 14, 10}}, // v2,v5,v4
            {9, {12, 7, 6}}    // v3,v4,v1
        };

        complex = atopo::create_complex_from_source(mesh);
    }

    static atopo::CellComplex complex;
};

atopo::CellComplex RP2IntegrationTest::complex;

TEST_F(RP2IntegrationTest, ComplexCreation) {
    ASSERT_EQ(complex.getNumberOfCells(0), 6);
    ASSERT_EQ(complex.getNumberOfCells(1), 15);
    ASSERT_EQ(complex.getNumberOfCells(2), 10);
}

TEST_P(RP2IntegrationTest, BoundaryOfBoundaryIsZero_PerFace) {
    int face_idx = GetParam(); 

    Eigen::SparseVector<int> face_vector(complex.getNumberOfCells(2));
    face_vector.insert(face_idx) = 1;
    atopo::Chain<int> single_face_chain(2, std::move(face_vector));

    auto boundary_of_face = atopo::boundary(complex, single_face_chain);
    auto boundary_of_boundary = atopo::boundary(complex, boundary_of_face);
    
    EXPECT_EQ(boundary_of_boundary.data.nonZeros(), 0) 
        << "Failure on face with index " << face_idx;
}

INSTANTIATE_TEST_SUITE_P(
    BoundaryCheck,
    RP2IntegrationTest,
    ::testing::Range(0, 10)
);

TEST_F(RP2IntegrationTest, HomologyCalculation) {
    auto homology_groups = complex.computeHomology();
    
    ASSERT_TRUE(homology_groups.count(0));
    EXPECT_EQ(homology_groups.at(0).rank, 1);
    EXPECT_TRUE(homology_groups.at(0).torsion_coeffs.empty());

    ASSERT_TRUE(homology_groups.count(1));
    EXPECT_EQ(homology_groups.at(1).rank, 0);
    
    auto torsion1 = homology_groups.at(1).torsion_coeffs;
    std::sort(torsion1.begin(), torsion1.end());
    ASSERT_EQ(torsion1.size(), 1);
    EXPECT_EQ(torsion1[0], 2);
    
    if (homology_groups.count(2)) {
        EXPECT_EQ(homology_groups.at(2).rank, 0);
        EXPECT_TRUE(homology_groups.at(2).torsion_coeffs.empty());
    }
}
