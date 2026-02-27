#include <gtest/gtest.h>
#include "test_common.h"
#include <atopo/query.h>

using namespace atopo;
using namespace atopo::test;

TEST(HomologyIntegration, SphereTopologicalFeatures) {
    auto sphere = make_sphere();
    auto complex = CellComplex::build(sphere);
    
    EXPECT_EQ(query::betti_number(complex, 0), 1);
    EXPECT_EQ(query::betti_number(complex, 1), 0);
    EXPECT_EQ(query::betti_number(complex, 2), 1);
    
    EXPECT_TRUE(query::is_connected(complex));
    EXPECT_TRUE(query::has_void(complex));
    EXPECT_FALSE(query::has_tunnel(complex));
}

TEST(HomologyIntegration, MultipleComponents) {
    auto c1 = make_circle();
    auto c2 = make_circle();
    
    MockMesh m;
    // Component 1: Vertices [0,1,2], Edges [0,1,2]
    for(auto v : c1.v) m.add_v(v.id);
    for(auto e : c1.e) m.add_e(e.id, e.v0, e.v1);
    
    // Component 2: Vertices [3,4,5], Edges [3,4,5]
    // Indices in boundary() must match the cumulative index in the vertex list.
    for(size_t i=0; i<3; ++i) m.add_v(i + 3);
    m.add_e(3, 3, 4); m.add_e(4, 4, 5); m.add_e(5, 5, 3);
    
    auto complex = CellComplex::build(m);
    EXPECT_EQ(query::count_connected_components(complex), 2);
    EXPECT_EQ(query::betti_number(complex, 1), 2);
}

TEST(HomologyIntegration, KleinBottleZ2) {
    CellComplex complex;
    complex.setCellCount(0, 1);
    complex.setCellCount(1, 2);
    complex.setCellCount(2, 1);

    atopo::IncidenceMatrix<atopo::IncidenceCoeff> d1(1, 2); 
    complex.setBoundaryOperator(1, std::move(d1));

    atopo::IncidenceMatrix<atopo::IncidenceCoeff> d2(2, 1);
    d2.insert(1, 0) = 2; // d2(f) = 2b
    complex.setBoundaryOperator(2, std::move(d2));

    // Over Z2, 2 becomes 0.
    // H0=Z2, H1=Z2+Z2, H2=Z2
    EXPECT_EQ(query::betti_number(complex, 0), 1);
    EXPECT_EQ(query::betti_number(complex, 1), 2);
    EXPECT_EQ(query::betti_number(complex, 2), 1);
}
