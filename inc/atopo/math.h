#ifndef ATOPO_MATH_H
#define ATOPO_MATH_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <ranges>

namespace atopo::math {

    /**
     * @brief A modern, C++20-based XOR vector engine for Z2 calculations.
     * Replaces libstick's sorted_boolean_vector.
     */
    struct Z2Vector {
        std::vector<int> indices;

        void add(const Z2Vector& other) {
            std::vector<int> next;
            next.reserve(indices.size() + other.indices.size());
            std::set_symmetric_difference(
                indices.begin(), indices.end(),
                other.indices.begin(), other.indices.end(),
                std::back_inserter(next)
            );
            indices = std::move(next);
        }

        bool empty() const { return indices.empty(); }
        int pivot() const { return indices.empty() ? -1 : indices.back(); }
    };

    /**
     * @brief Optimized Disjoint Set Union for fast graph rank calculation.
     */
    class DSU {
        std::vector<int> parent;
    public:
        DSU(std::size_t n) : parent(n) {
            std::iota(parent.begin(), parent.end(), 0);
        }

        int find(int i) {
            if (parent[i] == i) return i;
            return parent[i] = find(parent[i]);
        }

        bool unite(int i, int j) {
            int root_i = find(i);
            int root_j = find(j);
            if (root_i != root_j) {
                parent[root_i] = root_j;
                return true;
            }
            return false;
        }
    };

    /**
     * @brief Z2 View projection: transforms signed integers to Z2 (parity).
     */
    namespace views {
        inline auto z2 = std::views::filter([](int v) { return v % 2 != 0; });
    }

} // namespace atopo::math

#endif // ATOPO_MATH_H
