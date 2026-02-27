#include <gtest/gtest.h>
#include "test_common.h"

class OperatorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Simple 1x1 grid (4 vertices, 4 edges, 1 face)
        mesh.add_vertex(0); mesh.add_vertex(1); mesh.add_vertex(2); mesh.add_vertex(3);
        mesh.add_edge(0, 0, 1); mesh.add_edge(1, 1, 2); mesh.add_edge(2, 2, 3); mesh.add_edge(3, 3, 0);
        mesh.add_face(0, {0, 1, 2, 3}, {1, 1, 1, 1}); // Consistent orientation

        complex = atopo::create_complex_from_source(mesh);
    }

    LegacyMesh mesh;
    atopo::CellComplex complex;
};

TEST_F(OperatorTest, LaplacianMatrixSize) {
    auto L0 = atopo::laplacian_matrix<double>(complex, 0);
    EXPECT_EQ(L0.rows(), 4);
    EXPECT_EQ(L0.cols(), 4);

    auto L1 = atopo::laplacian_matrix<double>(complex, 1);
    EXPECT_EQ(L1.rows(), 4);
    EXPECT_EQ(L1.cols(), 4);
}

TEST_F(OperatorTest, ZeroFormLaplacian) {
    // 0-form laplacian on a graph is the standard graph Laplacian: L = delta_{-1}*d0 + d1*delta0 = d1*d1^T
    // Wait, for 0-forms, L0 = delta_{-1}*d0 + d1*delta0. Since d0 is usually considered 0, L0 = d1^T*d1.
    // Let's check our definition: L_p = d_{p+1} * delta_p + delta_{p-1} * d_p
    // For p=0: L0 = d1 * (d1)^T + delta_{-1} * d0. 
    // In our code: L = d1 * d1.transpose() if p=0. 
    // This is the vertex Laplacian.
    
    auto L0 = atopo::laplacian_matrix<double>(complex, 0);
    
    // Diagonal elements of Graph Laplacian d1*d1^T are the degrees of vertices
    // In our square, each vertex has degree 2.
    for (int i=0; i<4; ++i) {
        EXPECT_EQ(L0.coeff(i, i), 2.0);
    }
}

TEST_F(OperatorTest, LaplacianApplication) {
    // Apply L0 to a 0-cochain [1, 0, 0, 0]
    Eigen::SparseVector<double> v(4);
    v.insert(0) = 1.0;
    atopo::Cochain<double> c0(0, std::move(v));

    auto result = atopo::laplacian(complex, c0);
    EXPECT_EQ(result.dimension, 0);
    // Result should be the first column of L0: [2, -1, 0, -1] (depending on edge directions)
    EXPECT_EQ(result.data.coeff(0), 2.0);
}
