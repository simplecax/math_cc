#include "atopo.h"
#include "test_common.h" // For LegacyMesh traits

#include <string> // For std::to_string in exceptions

/**
 * @file atopo.cpp
 * @brief Contains definitions for non-templated functions in the atopo library.
 */
namespace atopo {

    void CellComplex::setCellCount(int dim, size_t count) { m_cell_counts[dim] = count; }
    void CellComplex::setBoundaryOperator(int source_dim, IncidenceMatrix<IncidenceCoeff>&& map) { m_boundary_maps[source_dim] = std::move(map); }
    size_t CellComplex::getNumberOfCells(int dim) const {
        auto it = m_cell_counts.find(dim);
        return (it == m_cell_counts.end()) ? 0 : it->second;
    }
    IncidenceMatrix<IncidenceCoeff> CellComplex::getIncidenceMap(int from_dim, int to_dim) const {
        if (to_dim == from_dim - 1) {
            auto it = m_boundary_maps.find(from_dim);
            if (it != m_boundary_maps.end()) {
                return it->second;
            }
        }
        if (to_dim == from_dim + 1) {
            auto it = m_boundary_maps.find(to_dim);
            if (it != m_boundary_maps.end()) return it->second.transpose();
        }
        return IncidenceMatrix<IncidenceCoeff>(getNumberOfCells(to_dim), getNumberOfCells(from_dim));
    }
    
    std::map<int, HomologyGroup> CellComplex::computeHomology() const {
        std::map<int, HomologyGroup> homology;
        int max_dim = 0;
        if(!m_cell_counts.empty()) { max_dim = m_cell_counts.rbegin()->first; }

        for (int p = 0; p <= max_dim + 1; ++p) {
             if (getNumberOfCells(p) == 0 && getNumberOfCells(p-1) == 0 && p > max_dim) continue;
            
            // --- CORRECTED NAMESPACE ---
            auto snf_dp = (p > 0) ? detail::compute_snf_results(getIncidenceMap(p, p - 1)) : atopo::HomologyGroup{};
            auto snf_dp1 = detail::compute_snf_results(getIncidenceMap(p + 1, p));
            
            long rank_Cp = getNumberOfCells(p);
            long rank_dp = snf_dp.rank;
            long rank_dp1 = snf_dp1.rank;

            HomologyGroup hg;
            hg.rank = rank_Cp - rank_dp - rank_dp1;
            hg.torsion_coeffs = snf_dp1.torsion_coeffs;
            
            if (hg.rank != 0 || !hg.torsion_coeffs.empty() || rank_Cp > 0) {
                homology[p] = hg;
            }
        }
        return homology;
    }
}
