#include <gtest/gtest.h>
#include "test_common.h"

class RP2IntegrationTest : public ::testing::TestWithParam<int> {
protected:
    static void SetUpTestSuite() {
        mesh.add_vertex(0); mesh.add_vertex(1); mesh.add_vertex(2); mesh.add_vertex(3); mesh.add_vertex(4); mesh.add_vertex(5);
        mesh.add_edge(0,0,1); mesh.add_edge(1,0,2); mesh.add_edge(2,0,3); mesh.add_edge(3,0,4); mesh.add_edge(4,0,5);
        mesh.add_edge(5,1,2); mesh.add_edge(6,1,3); mesh.add_edge(7,1,4); mesh.add_edge(8,1,5);
        mesh.add_edge(9,2,3); mesh.add_edge(10,2,4); mesh.add_edge(11,2,5);
        mesh.add_edge(12,3,4); mesh.add_edge(13,3,5);
        mesh.add_edge(14,4,5);

        // A standard, mathematically correct triangulation of RP2 with consistent orientations
        mesh.add_face(0, {0, 5, 1}, {1, 1, -1});   // v0v1, v1v2, v0v2 (reversed)
        mesh.add_face(1, {1, 9, 2}, {1, 1, -1});   // v1v2, v2v3, v1v3 (reversed)
        mesh.add_face(2, {2, 12, 3}, {1, 1, -1});  // v2v3, v3v4, v2v4 (reversed)
        mesh.add_face(3, {3, 14, 4}, {1, 1, -1});  // v3v4, v4v5, v3v5 (reversed)
        mesh.add_face(4, {4, 8, 0}, {-1, 1, 1});   // v4v5 (rev), v5v1, v0v1
        
        mesh.add_face(5, {5, 10, 7}, {1, 1, -1});  // v1v2, v2v4, v1v4 (reversed)
        mesh.add_face(6, {6, 13, 8}, {1, 1, -1});  // v1v3, v3v5, v1v5 (reversed)
        mesh.add_face(7, {9, 13, 11}, {1, 1, -1}); // v2v3, v3v5, v2v5 (reversed)
        mesh.add_face(8, {11, 14, 10}, {1, -1, -1});// v2v5, v4v5 (rev), v2v4 (reversed)
        mesh.add_face(9, {12, 7, 6}, {1, -1, 1});   // v3v4, v4v1 (rev), v1v3

        complex = atopo::create_complex_from_source(mesh);
    }

    static LegacyMesh mesh;
    static atopo::CellComplex complex;
};

LegacyMesh RP2IntegrationTest::mesh;
atopo::CellComplex RP2IntegrationTest::complex;

TEST_F(RP2IntegrationTest, ComplexCreation) {
    ASSERT_EQ(complex.getNumberOfCells(0), 6);
    ASSERT_EQ(complex.getNumberOfCells(1), 15);
    ASSERT_EQ(complex.getNumberOfCells(2), 10);
}

TEST_F(RP2IntegrationTest, HomologyCalculation) {
    auto homology_groups = complex.computeHomology();
    
    ASSERT_TRUE(homology_groups.count(0));
    EXPECT_EQ(homology_groups.at(0).rank, 1);
    ASSERT_TRUE(homology_groups.count(1));
    EXPECT_EQ(homology_groups.at(1).rank, 1); // Z2 rank
    ASSERT_TRUE(homology_groups.count(2));
    EXPECT_EQ(homology_groups.at(2).rank, 1); // Z2 rank
}

TEST_P(RP2IntegrationTest, BoundaryOfBoundaryIsZero_PerFace) {
    int face_idx = GetParam(); 
    Eigen::SparseVector<int> face_vector(complex.getNumberOfCells(2));
    face_vector.insert(face_idx) = 1;
    atopo::Chain<int> single_face_chain(2, std::move(face_vector));
    auto boundary_of_face = atopo::boundary(complex, single_face_chain);
    auto boundary_of_boundary = atopo::boundary(complex, boundary_of_face);
    EXPECT_EQ(boundary_of_boundary.data.nonZeros(), 0);
}

INSTANTIATE_TEST_SUITE_P(BoundaryCheck, RP2IntegrationTest, ::testing::Range(0, 10));
