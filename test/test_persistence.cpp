#include <gtest/gtest.h>
#include "test_common.h"
#include <atopo/persistence.h>

struct SimpleFiltration {
    auto ordered_cells(const LegacyMesh& mesh) const {
        std::vector<std::pair<int, size_t>> order;
        // Vertices first, then edges, then faces
        for (int d = 0; d <= 2; ++d) {
            for (size_t i = 0; i < mesh.cell_count(d); ++i) {
                order.push_back({d, i});
            }
        }
        return order;
    }
};

TEST(PersistenceTest, CirclePersistence) {
    LegacyMesh mesh;
    mesh.add_vertex(0); mesh.add_vertex(1); mesh.add_vertex(2);
    mesh.add_edge(0, 0, 1); mesh.add_edge(1, 1, 2); mesh.add_edge(2, 2, 0);

    SimpleFiltration filtration;
    auto diagram = atopo::compute_persistence(mesh, filtration);

    // Expected for a circle:
    // One persistent 0-cycle (connected component)
    // One persistent 1-cycle (the loop)
    
    int persistent_0 = 0;
    int persistent_1 = 0;
    for (const auto& p : diagram) {
        if (p.death_idx == -1) {
            if (p.dimension == 0) persistent_0++;
            if (p.dimension == 1) persistent_1++;
        }
    }

    EXPECT_EQ(persistent_0, 1);
    EXPECT_EQ(persistent_1, 1);
}
