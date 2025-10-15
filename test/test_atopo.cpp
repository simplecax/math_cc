#include <gtest/gtest.h>
#include "test_common.h"
#include <gtest/gtest.h>

// A test fixture class to set up a common environment for all tests.
// This avoids duplicating the mesh and complex creation code.
class AtopoTest : public ::testing::Test {
protected:
    // This function is called before each test is run
    void SetUp() override {
        // 1. Define a simple mesh (a single square face)
        my_mesh.vertices = {
            {{0}, {}}, {{1}, {}}, {{2}, {}}, {{3}, {}}
        };
        my_mesh.edges = {
            {0, 0, 1, {}}, {1, 1, 2, {}}, {2, 2, 3, {}}, {3, 3, 0, {}}
        };
        my_mesh.faces = {
            {0, {0, 1, 2, 3}, {}}
        };

        // 2. Build the CellComplex from the mesh
        complex = atopo::create_complex_from_source(my_mesh);
    }

    LegacyMesh my_mesh;
    atopo::CellComplex complex;
    using Coeff = int;
};

// Test case to verify that the complex is built correctly
TEST_F(AtopoTest, ComplexCreation) {
    ASSERT_EQ(complex.getNumberOfCells(0), 4); // 4 vertices
    ASSERT_EQ(complex.getNumberOfCells(1), 4); // 4 edges
    ASSERT_EQ(complex.getNumberOfCells(2), 1); // 1 face
}

// Test case to verify the fundamental property: the boundary of a boundary is zero (∂∂=0)
TEST_F(AtopoTest, BoundaryOfBoundaryIsZero) {
    // Create a 2-chain representing the single face
    Eigen::SparseVector<Coeff> face_vector(complex.getNumberOfCells(2));
    face_vector.insert(0) = 1;
    atopo::Chain<Coeff> face_chain(2, std::move(face_vector));

    // First boundary: ∂(face) -> edges
    auto boundary_of_face = atopo::boundary(complex, face_chain);
    
    // The boundary should be the sum of 4 edges.
    EXPECT_EQ(boundary_of_face.dimension, 1);
    EXPECT_EQ(boundary_of_face.data.nonZeros(), 4);

    // Second boundary: ∂(∂(face)) -> vertices
    auto boundary_of_boundary = atopo::boundary(complex, boundary_of_face);
    
    // The boundary of a boundary must be zero.
    EXPECT_EQ(boundary_of_boundary.dimension, 0);
    EXPECT_EQ(boundary_of_boundary.data.nonZeros(), 0);
}

// Test case to verify the fundamental property: the coboundary of a coboundary is zero (δδ=0)
TEST_F(AtopoTest, CoboundaryOfCoboundaryIsZero) {
    // Create a 0-cochain representing a potential difference (e.g., +1 at v0, -1 at v2)
    Eigen::SparseVector<Coeff> v_cochain_vector(complex.getNumberOfCells(0));
    v_cochain_vector.insert(0) = 1;
    v_cochain_vector.insert(2) = -1;
    atopo::Cochain<Coeff> vertex_cochain(0, std::move(v_cochain_vector));

    // First coboundary: δ(potential) -> "flow" on edges
    auto coboundary_of_v = atopo::coboundary(complex, vertex_cochain);
    
    // Expecting non-zero values on edges connected to v0 and v2.
    // Edge 0 (0->1), Edge 1 (1->2), Edge 2 (2->3), Edge 3 (3->0)
    // δ(f)(e0) = f(v1)-f(v0) = 0-1 = -1
    // δ(f)(e1) = f(v2)-f(v1) = -1-0 = -1
    // δ(f)(e2) = f(v3)-f(v2) = 0-(-1) = 1
    // δ(f)(e3) = f(v0)-f(v3) = 1-0 = 1
    EXPECT_EQ(coboundary_of_v.dimension, 1);
    EXPECT_EQ(coboundary_of_v.data.nonZeros(), 4);

    // Second coboundary: δ(δ(potential)) -> "circulation" on faces
    auto coboundary_of_coboundary = atopo::coboundary(complex, coboundary_of_v);

    // The coboundary of a coboundary must be zero.
    EXPECT_EQ(coboundary_of_coboundary.dimension, 2);
    EXPECT_EQ(coboundary_of_coboundary.data.nonZeros(), 0);
}