#include <gtest/gtest.h>
#include "test_common.h"
#include <atopo/operators.h>

using namespace atopo;
using namespace atopo::test;

class OperatorsUnit : public ::testing::Test {
protected:
    void SetUp() override {
        mesh = make_sphere();
        complex = CellComplex::build(mesh);
    }
    MockMesh mesh;
    CellComplex complex;
};

TEST_F(OperatorsUnit, BoundaryOfBoundaryIsZero) {
    // Face 0 chain
    Eigen::SparseVector<int> f0(4);
    f0.insert(0) = 1;
    Chain<int> c2(2, std::move(f0));

    auto d2_c2 = boundary(complex, c2);
    auto d1_d2_c2 = boundary(complex, d2_c2);
    
    EXPECT_EQ(d1_d2_c2.data.nonZeros(), 0);
}

TEST_F(OperatorsUnit, CoboundaryOfCoboundaryIsZero) {
    // Vertex 0 cochain
    Eigen::SparseVector<int> v0(4);
    v0.insert(0) = 1;
    Cochain<int> co0(0, std::move(v0));

    auto delta0_co0 = coboundary(complex, co0);
    auto delta1_delta0_co0 = coboundary(complex, delta0_co0);
    
    EXPECT_EQ(delta1_delta0_co0.data.nonZeros(), 0);
}

TEST_F(OperatorsUnit, LaplacianMatrix) {
    auto L0 = laplacian_matrix<double>(complex, 0);
    EXPECT_EQ(L0.rows(), 4);
    // On sphere, each vertex has degree 3
    for (int i=0; i<4; ++i) {
        EXPECT_EQ(L0.coeff(i, i), 3.0);
    }
}
