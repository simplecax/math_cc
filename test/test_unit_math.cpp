#include <gtest/gtest.h>
#include <atopo/math.h>

using namespace atopo::math;

TEST(MathUnit, Z2VectorXOR) {
    Z2Vector a{{0, 1, 2}};
    Z2Vector b{{1, 2, 3}};
    
    a.add(b);
    // XOR results in symmetric difference: {0, 3}
    EXPECT_EQ(a.indices.size(), 2);
    EXPECT_EQ(a.indices[0], 0);
    EXPECT_EQ(a.indices[1], 3);
    
    a.add(b);
    // {0, 3} XOR {1, 2, 3} = {0, 1, 2}
    EXPECT_EQ(a.indices.size(), 3);
    EXPECT_EQ(a.indices[0], 0);
    EXPECT_EQ(a.indices[1], 1);
    EXPECT_EQ(a.indices[2], 2);
}

TEST(MathUnit, DSUConnectivity) {
    DSU dsu(5);
    EXPECT_TRUE(dsu.unite(0, 1));
    EXPECT_TRUE(dsu.unite(1, 2));
    EXPECT_FALSE(dsu.unite(0, 2)); // Already connected
    
    EXPECT_EQ(dsu.find(0), dsu.find(2));
    EXPECT_NE(dsu.find(0), dsu.find(3));
}

TEST(MathUnit, Z2View) {
    std::vector<int> vals = {0, 1, 2, 3, 4, 5};
    auto filtered = vals | views::z2;
    
    std::vector<int> result;
    for (int v : filtered) result.push_back(v);
    
    ASSERT_EQ(result.size(), 3);
    EXPECT_EQ(result[0], 1);
    EXPECT_EQ(result[1], 3);
    EXPECT_EQ(result[2], 5);
}
