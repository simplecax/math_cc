#include "atopo.h"
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <numeric>

namespace atopo::detail {

    /**
     * @brief Union-Find (Disjoint Set Union) for fast rank calculation.
     */
    struct DSU {
        std::vector<int> parent;
        int components;
        DSU(int n) : parent(n), components(n) {
            std::iota(parent.begin(), parent.end(), 0);
        }
        int find(int i) {
            if (parent[i] == i) return i;
            return parent[i] = find(parent[i]);
        }
        void unite(int i, int j) {
            int root_i = find(i);
            int root_j = find(j);
            if (root_i != root_j) {
                parent[root_i] = root_j;
                components--;
            }
        }
    };

    /**
     * @brief Computes rank using Union-Find for matrices with <= 2 ones per column.
     * In Z2, the rank of a graph incidence matrix is V - (number of components).
     * If a component has an "open" edge (column with only one 1), it doesn't reduce rank.
     */
    size_t computeGraphRankZ2(const IncidenceMatrix<IncidenceCoeff>& mat) {
        int numRows = static_cast<int>(mat.rows());
        DSU dsu(numRows);
        std::vector<bool> has_boundary_edge(numRows, false);

        for (int j = 0; j < mat.cols(); ++j) {
            std::vector<int> rows;
            for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(mat, j); it; ++it) {
                if (it.value() % 2 != 0) rows.push_back(it.row());
            }

            if (rows.size() == 2) {
                dsu.unite(rows[0], rows[1]);
            } else if (rows.size() == 1) {
                has_boundary_edge[rows[0]] = true;
            }
        }

        // A component is "closed" if none of its vertices connect to a single-entry column.
        std::unordered_map<int, bool> component_is_open;
        for (int i = 0; i < numRows; ++i) {
            int root = dsu.find(i);
            if (has_boundary_edge[i]) component_is_open[root] = true;
        }

        int closed_components = 0;
        std::vector<bool> root_seen(numRows, false);
        for (int i = 0; i < numRows; ++i) {
            int root = dsu.find(i);
            if (!root_seen[root]) {
                root_seen[root] = true;
                if (!component_is_open[root]) closed_components++;
            }
        }

        return numRows - closed_components;
    }

    /**
     * @brief Enhanced Z2 solver with fast-paths for CAD/Manifold structures.
     */
    HomologyGroup compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat) {
        if (sparse_mat.nonZeros() == 0) return {};

        int numRows = static_cast<int>(sparse_mat.rows());
        int numCols = static_cast<int>(sparse_mat.cols());

        // --- Fast Path Check ---
        bool graph_like_cols = true;
        for (int j = 0; j < numCols; ++j) {
            int count = 0;
            for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(sparse_mat, j); it; ++it) {
                if (it.value() % 2 != 0) count++;
            }
            if (count > 2) { graph_like_cols = false; break; }
        }

        if (graph_like_cols) {
            HomologyGroup res;
            res.rank = computeGraphRankZ2(sparse_mat);
            return res;
        }

        // --- Fallback: Optimized Gaussian Elimination ---
        using Z2Vector = std::vector<int>;
        std::unordered_map<int, Z2Vector> pivot_cols; // pivot_row -> reduced_column

        for (int j = 0; j < numCols; ++j) {
            Z2Vector curr;
            for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(sparse_mat, j); it; ++it) {
                if (it.value() % 2 != 0) curr.push_back(it.row());
            }
            std::sort(curr.begin(), curr.end());

            while (!curr.empty()) {
                int pivot = curr.back();
                if (pivot_cols.find(pivot) == pivot_cols.end()) {
                    pivot_cols[pivot] = std::move(curr);
                    break;
                } else {
                    const Z2Vector& other = pivot_cols[pivot];
                    Z2Vector next_curr;
                    std::set_symmetric_difference(curr.begin(), curr.end(),
                                                other.begin(), other.end(),
                                                std::back_inserter(next_curr));
                    curr = std::move(next_curr);
                }
            }
        }

        HomologyGroup result;
        result.rank = pivot_cols.size();
        return result;
    }

} // namespace atopo::detail
