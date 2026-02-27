#include <gtest/gtest.h>
#include "test_common.h"
#include <atopo/persistence.h>

using namespace atopo;
using namespace atopo::test;

struct StepFiltration {
    auto ordered_cells(const MockMesh& mesh) const {
        std::vector<std::pair<int, size_t>> order;
        // Vertices first
        for (size_t i = 0; i < mesh.cell_count(0); ++i) order.push_back({0, i});
        // Then edges
        for (size_t i = 0; i < mesh.cell_count(1); ++i) order.push_back({1, i});
        // Then faces
        for (size_t i = 0; i < mesh.cell_count(2); ++i) order.push_back({2, i});
        return order;
    }
};

TEST(PersistenceIntegration, CircleFiltration) {
    auto circle = make_circle();
    StepFiltration filtration;
    auto diagram = compute_persistence(circle, filtration);
    
    // Summary check
    std::string summary = summarize_persistence(diagram);
    EXPECT_TRUE(summary.find("Persistence Summary") != std::string::npos);

    int essential_0 = 0;
    int essential_1 = 0;
    for (const auto& p : diagram) {
        if (p.death_idx == -1) {
            if (p.dimension == 0) essential_0++;
            if (p.dimension == 1) essential_1++;
        }
    }
    
    // One connected component, one loop
    EXPECT_EQ(essential_0, 1);
    EXPECT_EQ(essential_1, 1);
    
    // Check birth/death details:
    // First 3 cells (vertices) should give birth to 3 0-cycles.
    // Next 2 edges should kill 2 of them.
    // Last edge should give birth to a 1-cycle.
    int killed_0 = 0;
    for (const auto& p : diagram) {
        if (p.dimension == 0 && p.death_idx != -1) killed_0++;
    }
    EXPECT_EQ(killed_0, 2);
}

TEST(PersistenceIntegration, EmptyComplex) {
    MockMesh empty;
    StepFiltration f;
    auto diagram = compute_persistence(empty, f);
    EXPECT_TRUE(diagram.empty());
}
