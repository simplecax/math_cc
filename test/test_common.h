#ifndef ATOPO_TEST_COMMON_H
#define ATOPO_TEST_COMMON_H

#include <atopo.h>
#include <vector>
#include <ranges>
#include <variant>
#include <memory>

namespace atopo::test {

    // --- Mock User Structures ---
    struct MockBoundaryEntry {
        std::size_t index;
        int coefficient;
    };

    struct MockVertex {
        std::size_t id;
        int dimension() const { return 0; }
        std::vector<MockBoundaryEntry> boundary() const { return {}; }
    };

    struct MockEdge {
        std::size_t id;
        std::size_t v0, v1;
        int dimension() const { return 1; }
        std::vector<MockBoundaryEntry> boundary() const {
            return {{v0, -1}, {v1, 1}};
        }
    };

    struct MockFace {
        std::size_t id;
        std::vector<std::size_t> edges;
        std::vector<int> orientations;
        int dimension() const { return 2; }
        std::vector<MockBoundaryEntry> boundary() const {
            std::vector<MockBoundaryEntry> result;
            for (size_t i = 0; i < edges.size(); ++i) {
                result.push_back({edges[i], (i < orientations.size() ? orientations[i] : 1)});
            }
            return result;
        }
    };

    using MockCellVariant = std::variant<MockVertex, MockEdge, MockFace>;

    struct MockCell {
        MockCellVariant data;
        int dimension() const { return std::visit([](const auto& c) { return c.dimension(); }, data); }
        auto boundary() const { return std::visit([](const auto& c) { return c.boundary(); }, data); }
    };

    struct MockMesh {
        std::vector<MockVertex> v;
        std::vector<MockEdge> e;
        std::vector<MockFace> f;

        int max_dim() const { return 2; }
        std::size_t cell_count(int dim) const {
            if (dim == 0) return v.size();
            if (dim == 1) return e.size();
            if (dim == 2) return f.size();
            return 0;
        }
        std::vector<MockCell> cells(int dim) const {
            std::vector<MockCell> result;
            if (dim == 0) for(const auto& item : v) result.push_back({item});
            else if (dim == 1) for(const auto& item : e) result.push_back({item});
            else if (dim == 2) for(const auto& item : f) result.push_back({item});
            return result;
        }

        // Helpers
        void add_v(size_t id) { v.push_back({id}); }
        void add_e(size_t id, size_t v0, size_t v1) { e.push_back({id, v0, v1}); }
        void add_f(size_t id, std::vector<size_t> es, std::vector<int> o = {}) { f.push_back({id, es, o}); }
    };

    // --- Standard Topological Factories ---
    inline MockMesh make_circle() {
        MockMesh m;
        m.add_v(0); m.add_v(1); m.add_v(2);
        m.add_e(0, 0, 1); m.add_e(1, 1, 2); m.add_e(2, 2, 0);
        return m;
    }

    inline MockMesh make_sphere() {
        MockMesh m;
        m.add_v(0); m.add_v(1); m.add_v(2); m.add_v(3);
        m.add_e(0, 0, 1); m.add_e(1, 0, 2); m.add_e(2, 0, 3);
        m.add_e(3, 1, 2); m.add_e(4, 1, 3); m.add_e(5, 2, 3);
        m.add_f(0, {0, 3, 1}, {1, 1, -1});
        m.add_f(1, {0, 4, 2}, {1, 1, -1});
        m.add_f(2, {1, 5, 2}, {1, 1, -1});
        m.add_f(3, {3, 5, 4}, {1, 1, -1});
        return m;
    }

} // namespace atopo::test

#endif
