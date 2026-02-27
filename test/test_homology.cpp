#include <gtest/gtest.h>
#include "test_common.h"

// Test fixture for creating complexes from user-defined mesh data
class HomologyTest : public ::testing::Test {
protected:
    using Coeff = int;
    
    void RunTest(const LegacyMesh& mesh, const std::map<int, atopo::HomologyGroup>& expected_homology) {
        auto complex = atopo::create_complex_from_source(mesh);
        auto homology_groups = complex.computeHomology();

        for(const auto& pair : expected_homology) {
            int dim = pair.first;
            const auto& expected = pair.second;

            ASSERT_TRUE(homology_groups.count(dim)) << "Homology group for dimension " << dim << " is missing.";
            const auto& actual = homology_groups.at(dim);

            EXPECT_EQ(actual.rank, expected.rank) << "Betti number for dimension " << dim << " is incorrect.";
            
            // Sort both vectors to ensure comparison is order-independent
            auto actual_torsion = actual.torsion_coeffs;
            auto expected_torsion = expected.torsion_coeffs;
            std::sort(actual_torsion.begin(), actual_torsion.end());
            std::sort(expected_torsion.begin(), expected_torsion.end());

            EXPECT_EQ(actual_torsion, expected_torsion) << "Torsion coefficients for dimension " << dim << " are incorrect.";
        }
    }
};

TEST_F(HomologyTest, SingleEdge) {
    LegacyMesh mesh;
    mesh.add_vertex(0); mesh.add_vertex(1);
    mesh.add_edge(0, 0, 1);
    RunTest(mesh, {{0, {1, {}}}, {1, {0, {}}}});
}

TEST_F(HomologyTest, TwoDisconnectedEdges) {
    LegacyMesh mesh;
    mesh.add_vertex(0); mesh.add_vertex(1); mesh.add_vertex(2); mesh.add_vertex(3);
    mesh.add_edge(0, 0, 1); mesh.add_edge(1, 2, 3);
    RunTest(mesh, {{0, {2, {}}}, {1, {0, {}}}});
}

TEST_F(HomologyTest, CircleLoop) {
    LegacyMesh mesh;
    mesh.add_vertex(0); mesh.add_vertex(1); mesh.add_vertex(2);
    mesh.add_edge(0, 0, 1); mesh.add_edge(1, 1, 2); mesh.add_edge(2, 2, 0);
    RunTest(mesh, {{0, {1, {}}}, {1, {1, {}}}});
}

TEST_F(HomologyTest, Sphere) {
    LegacyMesh mesh;
    mesh.add_vertex(0); mesh.add_vertex(1); mesh.add_vertex(2); mesh.add_vertex(3);
    mesh.add_edge(0, 0, 1); mesh.add_edge(1, 0, 2); mesh.add_edge(2, 0, 3);
    mesh.add_edge(3, 1, 2); mesh.add_edge(4, 1, 3); mesh.add_edge(5, 2, 3);
    mesh.add_face(0, {0, 3, 1}); 
    mesh.add_face(1, {0, 4, 2}); 
    mesh.add_face(2, {1, 5, 2}); 
    mesh.add_face(3, {3, 5, 4}); 
    RunTest(mesh, {{0, {1, {}}}, {1, {0, {}}}, {2, {1, {}}}});
}
