#include "atopo.h"
#include <atopo/math.h>
#include <unordered_map>

namespace atopo::detail {

    /**
     * @brief Computes rank using DSU for matrices with <= 2 ones per column.
     * In Z2, the rank of a graph incidence matrix is V - (number of components).
     */
    size_t computeGraphRankZ2(const IncidenceMatrix<IncidenceCoeff>& mat) {
        int numRows = static_cast<int>(mat.rows());
        math::DSU dsu(numRows);
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

        return static_cast<size_t>(numRows - closed_components);
    }

    /**
     * @brief Enhanced Z2 solver with modernized C++20 DSA core.
     */
    HomologyGroup compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat) {
        if (sparse_mat.nonZeros() == 0) return {};

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

        // --- Unified Gaussian Elimination Core ---
        std::unordered_map<int, math::Z2Vector> pivot_cols;

        for (int j = 0; j < numCols; ++j) {
            math::Z2Vector curr;
            for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(sparse_mat, j); it; ++it) {
                if (it.value() % 2 != 0) curr.indices.push_back(it.row());
            }
            std::sort(curr.indices.begin(), curr.indices.end());

            while (!curr.empty()) {
                int pivot = curr.pivot();
                auto it_pivot = pivot_cols.find(pivot);
                if (it_pivot == pivot_cols.end()) {
                    pivot_cols[pivot] = std::move(curr);
                    break;
                } else {
                    curr.add(it_pivot->second);
                }
            }
        }

        HomologyGroup result;
        result.rank = pivot_cols.size();
        return result;
    }

} // namespace atopo::detail
