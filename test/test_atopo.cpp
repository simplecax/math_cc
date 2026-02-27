#include <gtest/gtest.h>
#include "test_common.h"

// A test fixture class to set up a common environment for all tests.
class AtopoTest : public ::testing::Test {
protected:
    void SetUp() override {
        my_mesh.add_vertex(0); my_mesh.add_vertex(1); my_mesh.add_vertex(2); my_mesh.add_vertex(3);
        my_mesh.add_edge(0, 0, 1); my_mesh.add_edge(1, 1, 2); my_mesh.add_edge(2, 2, 3); my_mesh.add_edge(3, 3, 0);
        my_mesh.add_face(0, {0, 1, 2, 3});

        complex = atopo::create_complex_from_source(my_mesh);
    }

    LegacyMesh my_mesh;
    atopo::CellComplex complex;
    using Coeff = int;
};

TEST_F(AtopoTest, ComplexCreation) {
    ASSERT_EQ(complex.getNumberOfCells(0), 4);
    ASSERT_EQ(complex.getNumberOfCells(1), 4);
    ASSERT_EQ(complex.getNumberOfCells(2), 1);
}

TEST_F(AtopoTest, BoundaryOfBoundaryIsZero) {
    Eigen::SparseVector<Coeff> face_vector(complex.getNumberOfCells(2));
    face_vector.insert(0) = 1;
    atopo::Chain<Coeff> face_chain(2, std::move(face_vector));

    auto boundary_of_face = atopo::boundary(complex, face_chain);
    EXPECT_EQ(boundary_of_face.dimension, 1);
    EXPECT_EQ(boundary_of_face.data.nonZeros(), 4);

    auto boundary_of_boundary = atopo::boundary(complex, boundary_of_face);
    EXPECT_EQ(boundary_of_boundary.dimension, 0);
    EXPECT_EQ(boundary_of_boundary.data.nonZeros(), 0);
}

TEST_F(AtopoTest, CoboundaryOfCoboundaryIsZero) {
    Eigen::SparseVector<Coeff> v_cochain_vector(complex.getNumberOfCells(0));
    v_cochain_vector.insert(0) = 1;
    v_cochain_vector.insert(2) = -1;
    atopo::Cochain<Coeff> vertex_cochain(0, std::move(v_cochain_vector));

    auto coboundary_of_v = atopo::coboundary(complex, vertex_cochain);
    EXPECT_EQ(coboundary_of_v.dimension, 1);
    EXPECT_EQ(coboundary_of_v.data.nonZeros(), 4);

    auto coboundary_of_coboundary = atopo::coboundary(complex, coboundary_of_v);
    EXPECT_EQ(coboundary_of_coboundary.dimension, 2);
    EXPECT_EQ(coboundary_of_coboundary.data.nonZeros(), 0);
}
