#!/bin/bash
#
# This script applies the definitive fix to the project structure, resolving
# the redefinition and namespace errors from the previous script.
#
# USAGE: ./final_fix.sh from the project root.

set -e

echo "🚀 Correcting project structure to resolve redefinition errors..."

# ==============================================================================
# 1. DELETE the problematic .inc file
# ==============================================================================
rm -f src/atopo_templates.inc
echo " -> Removed src/atopo_templates.inc"

# ==============================================================================
# 2. UPDATE inc/atopo/operators.h with the full template definitions
# ==============================================================================
cat << 'EOF' > inc/atopo/operators.h
#ifndef ATOPO_OPERATORS_H
#define ATOPO_OPERATORS_H

#include <atopo/complex.h>
#include <atopo/chain.h>
#include <atopo/traits.h>

namespace atopo {
    // --- TEMPLATE DEFINITIONS MOVED HERE ---
    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex& complex, const Chain<T>& chain) {
        if (chain.dimension <= 0) return Chain<T>(-1, Eigen::SparseVector<T>());
        const IncidenceMatrix<IncidenceCoeff>& d_p = complex.getIncidenceMap(chain.dimension, chain.dimension - 1);
        Eigen::SparseVector<T> boundary_vector = d_p.template cast<T>() * chain.data;
        boundary_vector.prune(0, 0);
        return Chain<T>(chain.dimension - 1, std::move(boundary_vector));
    }

    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex& complex, const Cochain<T>& cochain) {
        int p = cochain.dimension;
        const IncidenceMatrix<IncidenceCoeff> delta_p = complex.getIncidenceMap(p + 1, p).transpose();
        Eigen::SparseVector<T> coboundary_vector = delta_p.template cast<T>() * cochain.data;
        coboundary_vector.prune(0, 0);
        return Cochain<T>(p + 1, std::move(coboundary_vector));
    }

    template<typename TopologySource>
    [[nodiscard]] CellComplex create_complex_from_source(const TopologySource& source) {
        return TopologySourceTraits<TopologySource>::build(source);
    }
}
#endif // ATOPO_OPERATORS_H
EOF
echo " -> Moved template definitions to inc/atopo/operators.h"

# ==============================================================================
# 3. UPDATE src/atopo.cpp to be lean and correct
# ==============================================================================
cat << 'EOF' > src/atopo.cpp
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
            if (it != m_boundary_maps.end()) return it->second;
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
    
    CellComplex TopologySourceTraits<LegacyMesh>::build(const LegacyMesh& mesh) {
        CellComplex complex;
        complex.setCellCount(0, mesh.vertices.size());
        complex.setCellCount(1, mesh.edges.size());
        complex.setCellCount(2, mesh.faces.size());
        
        IncidenceMatrix<IncidenceCoeff> d1(mesh.vertices.size(), mesh.edges.size());
        for (const auto& edge : mesh.edges) {
            d1.insert(edge.v_start, edge.id) = -1;
            d1.insert(edge.v_end, edge.id) = 1;
        }
        complex.setBoundaryOperator(1, std::move(d1));

        IncidenceMatrix<IncidenceCoeff> d2(mesh.edges.size(), mesh.faces.size());
        for (const auto& face : mesh.faces) {
            if (face.edge_loop.size() < 2) continue;

            const auto& e0 = mesh.edges.at(face.edge_loop[0]);
            const auto& e1 = mesh.edges.at(face.edge_loop[1]);

            size_t v_start;
            if (e0.v_end == e1.v_start || e0.v_end == e1.v_end) { v_start = e0.v_start; }
            else { v_start = e0.v_end; }
            
            size_t current_v = v_start;
            for (size_t edge_id : face.edge_loop) {
                const auto& edge = mesh.edges.at(edge_id);
                if (edge.v_start == current_v) {
                    d2.insert(edge.id, face.id) = 1;
                    current_v = edge.v_end;
                } else if (edge.v_end == current_v) {
                    d2.insert(edge.id, face.id) = -1;
                    current_v = edge.v_start;
                } else {
                    throw std::runtime_error("Face " + std::to_string(face.id) + " has a discontinuous edge loop.");
                }
            }
            if (current_v != v_start) {
                throw std::runtime_error("Face " + std::to_string(face.id) + " has an open edge loop.");
            }
        }
        complex.setBoundaryOperator(2, std::move(d2));
        return complex;
    }
}
EOF
echo " -> Corrected and simplified src/atopo.cpp"

echo "✅ Project structure has been corrected."
echo "   Please clean and rebuild."
echo "   Run: rm -rf build && export PKG_CONFIG_PATH=/opt/lib/pkgconfig:\$PKG_CONFIG_PATH && cmake -S . -B build && cmake --build build && ctest --test-dir build"