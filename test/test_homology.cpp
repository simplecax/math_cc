#include <gtest/gtest.h>
#include "atopo.h"

// Test fixture for creating complexes from user-defined mesh data
class HomologyTest : public ::testing::Test {
protected:
    using Coeff = int;
    
    // Helper to build and test a complex
    void RunTest(const LegacyMesh& mesh, const std::map<int, int>& expected_betti) {
        auto complex = atopo::create_complex_from_source(mesh);
        auto betti_numbers = complex.computeBettiNumbers();

        // Check that all expected Betti numbers are present and correct
        for(const auto& pair : expected_betti) {
            ASSERT_TRUE(betti_numbers.count(pair.first)) << "Betti number for dimension " << pair.first << " is missing.";
            EXPECT_EQ(betti_numbers.at(pair.first), pair.second) << "Betti number for dimension " << pair.first << " is incorrect.";
        }
        
        // Check that no unexpected Betti numbers are present
        for(const auto& pair : betti_numbers) {
            if(pair.second != 0) { // Only check for non-zero unexpected values
                 ASSERT_TRUE(expected_betti.count(pair.first)) << "Unexpected non-zero Betti number at dimension " << pair.first;
            }
        }
    }
};

TEST_F(HomologyTest, SingleEdge) {
    LegacyMesh mesh;
    mesh.vertices = {{0, {}}, {1, {}}};
    mesh.edges = {{0, 0, 1, {}}};
    // Expected: 1 connected component (β₀=1), 0 loops (β₁=0)
    RunTest(mesh, {{0, 1}, {1, 0}});
}

TEST_F(HomologyTest, TwoDisconnectedEdges) {
    LegacyMesh mesh;
    mesh.vertices = {{0, {}}, {1, {}}, {2, {}}, {3, {}}};
    mesh.edges = {{0, 0, 1, {}}, {1, 2, 3, {}}};
    // Expected: 2 connected components (β₀=2), 0 loops (β₁=0)
    RunTest(mesh, {{0, 2}, {1, 0}});
}

TEST_F(HomologyTest, CircleLoop) {
    LegacyMesh mesh;
    mesh.vertices = {{0, {}}, {1, {}}, {2, {}}};
    mesh.edges = {{0, 0, 1, {}}, {1, 1, 2, {}}, {2, 2, 0, {}}};
    // Expected: 1 connected component (β₀=1), 1 loop (β₁=1)
    RunTest(mesh, {{0, 1}, {1, 1}});
}

TEST_F(HomologyTest, Sphere) {
    // A sphere represented as a tetrahedron
    LegacyMesh mesh;
    mesh.vertices = {{0,{}}, {1,{}}, {2,{}}, {3,{}}};
    // Edges of a complete graph K4
    mesh.edges = {
        {0,0,1,{}}, {1,0,2,{}}, {2,0,3,{}},
        {3,1,2,{}}, {4,1,3,{}}, {5,2,3,{}}
    };
    // --- CORRECTED FACE DEFINITIONS ---
    // The four valid triangular faces of the tetrahedron
    mesh.faces = {
        {0, {0, 3, 1}, {}}, // Face on vertices {0,1,2}
        {1, {0, 4, 2}, {}}, // Face on vertices {0,1,3}
        {2, {1, 5, 2}, {}}, // Face on vertices {0,2,3}
        {3, {3, 5, 4}, {}}  // Face on vertices {1,2,3}
    };
    // Expected: 1 component (β₀=1), 0 loops (β₁=0), 1 void (β₂=1)
    RunTest(mesh, {{0, 1}, {1, 0}, {2, 1}});
}
