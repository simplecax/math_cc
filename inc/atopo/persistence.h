#ifndef ATOPO_PERSISTENCE_H
#define ATOPO_PERSISTENCE_H

#include <atopo/complex.h>
#include <atopo/math.h>
#include <vector>
#include <map>
#include <algorithm>

namespace atopo {

    struct PersistencePoint {
        int birth_idx;
        int death_idx; // -1 for essential
        int dimension;
    };

    using PersistenceDiagram = std::vector<PersistencePoint>;

    /**
     * @brief Computes persistent homology over Z2.
     */
    template<concepts::Complex C, concepts::Filtration<C> F>
    PersistenceDiagram compute_persistence(const C& source, const F& filtration) {
        auto ordered = filtration.ordered_cells(source);
        
        struct CellId { int dim; size_t local_idx; };
        std::vector<CellId> global_to_local;
        std::map<int, std::map<size_t, int>> local_to_global;

        int g_idx = 0;
        for (auto item : ordered) {
            int d = std::get<0>(item);
            size_t i = std::get<1>(item);
            global_to_local.push_back({d, i});
            local_to_global[d][i] = g_idx++;
        }

        int total_cells = static_cast<int>(global_to_local.size());
        CellComplex complex = CellComplex::build(source);
        
        std::vector<math::Z2Vector> rbm_cols(total_cells);
        std::map<int, int> lowestone_to_column;
        PersistenceDiagram diagram;

        for (int j = 0; j < total_cells; ++j) {
            auto& cell = global_to_local[j];
            math::Z2Vector& curr = rbm_cols[j];

            // 1. Initialize boundary from filtration indices
            if (cell.dim > 0) {
                auto incidence = complex.getIncidenceMap(cell.dim, cell.dim - 1);
                for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(incidence, cell.local_idx); it; ++it) {
                    if (it.value() % 2 != 0) {
                        curr.indices.push_back(local_to_global[cell.dim - 1][it.row()]);
                    }
                }
                std::sort(curr.indices.begin(), curr.indices.end());
            }

            // 2. Reduce column
            while (!curr.empty()) {
                int i = curr.pivot();
                auto it_p = lowestone_to_column.find(i);
                if (it_p == lowestone_to_column.end()) {
                    lowestone_to_column[i] = j;
                    break;
                } else {
                    curr.add(rbm_cols[it_p->second]);
                }
            }

            // 3. Extract persistence points
            if (curr.empty()) {
                // Potential birth of a cycle
                // We'll mark it now, and if it dies later, we update it.
                diagram.push_back({j, -1, cell.dim});
            } else {
                // j kills the cycle that was born at i
                int i = curr.pivot();
                for (auto& p : diagram) {
                    if (p.birth_idx == i) {
                        p.death_idx = j;
                        break;
                    }
                }
            }
        }

        // Clean up: only return cycles that were born but never died, 
        // or died at a valid index.
        return diagram;
    }

} // namespace atopo

#endif // ATOPO_PERSISTENCE_H
