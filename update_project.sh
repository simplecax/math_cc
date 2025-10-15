#!/bin/bash
#
# This script applies the definitive fix to inc/atopo.h, restoring the
# full definition of the ChainBase struct that was corrupted by a previous script.
#
# USAGE: ./final_fix.sh from the project root.

set -e

echo "🚀 Correcting the corrupted atopo.h file..."

# Overwrite the header file with a completely corrected version.
cat << 'EOF' > inc/atopo.h
#ifndef ATOPO_H
#define ATOPO_H

#include <iostream>
#include <vector>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <map>
#include <utility>
#include <algorithm>
#include <numeric>

#include <Eigen/Sparse>

namespace atopo {

    // --- Core Data Types ---
    template<typename T> using Coefficient = T;
    using DimensionPair = std::pair<int, int>;
    using IncidenceCoeff = signed short;
    template<typename T> using IncidenceMatrix = Eigen::SparseMatrix<T>;

    struct HomologyGroup {
        int rank = 0;
        std::vector<long> torsion_coeffs;
    };

    // --- TRAITS ---
    template<typename CellType>
    struct CellTraits {
        static constexpr int dimension = -1;
        using GeometryType = void;
        static GeometryType get_geometry(const CellType& cell);
    };
    template<typename TopologySource>
    struct TopologySourceTraits { /* User specializes */ };

    // --- Implementation Detail for Homology ---
    namespace detail {
        std::vector<long> compute_invariant_factors(const IncidenceMatrix<IncidenceCoeff>& sparse_mat);
    }

    // --- CORRECTED DEFINITIONS ---
    template<typename T>
    struct ChainBase {
        int dimension;
        Eigen::SparseVector<Coefficient<T>> data;
        ChainBase(int p, const Eigen::SparseVector<Coefficient<T>>& vec) : dimension(p), data(vec) {}
        ChainBase(int p, Eigen::SparseVector<Coefficient<T>>&& vec) : dimension(p), data(std::move(vec)) {}
    };
    template<typename T> struct Chain : public ChainBase<T> { using ChainBase<T>::ChainBase; };
    template<typename T> struct Cochain : public ChainBase<T> { using ChainBase<T>::ChainBase; };

    class CellComplex {
    private:
        std::map<int, size_t> m_cell_counts;
        std::map<int, IncidenceMatrix<IncidenceCoeff>> m_boundary_maps;
        std::map<DimensionPair, IncidenceMatrix<IncidenceCoeff>> m_general_incidence_maps;
    public:
        void setCellCount(int dim, size_t count) { m_cell_counts[dim] = count; }
        void setBoundaryOperator(int source_dim, IncidenceMatrix<IncidenceCoeff>&& map) { m_boundary_maps[source_dim] = std::move(map); }
        void setGeneralIncidenceMap(int from_dim, int to_dim, IncidenceMatrix<IncidenceCoeff>&& map);
        [[nodiscard]] size_t getNumberOfCells(int dim) const;
        [[nodiscard]] IncidenceMatrix<IncidenceCoeff> getIncidenceMap(int from_dim, int to_dim) const;
        
        [[nodiscard]] std::map<int, HomologyGroup> computeHomology() const;
    };
    
    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex& complex, const Chain<T>& chain);
    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex& complex, const Cochain<T>& cochain);
    template<typename TopologySource>
    [[nodiscard]] CellComplex create_complex_from_source(const TopologySource& source);
}

// Glue code and example implementations
struct MyPoint3D { double x, y, z; }; struct MyCurve   { /* ... */ }; struct MySurface { /* ... */ };
struct LegacyVertex { size_t id; MyPoint3D geometry; }; struct LegacyEdge   { size_t id; size_t v_start, v_end; MyCurve geometry; };
struct LegacyFace   { size_t id; std::vector<size_t> edge_loop; MySurface geometry; };
struct LegacyMesh { std::vector<LegacyVertex> vertices; std::vector<LegacyEdge> edges; std::vector<LegacyFace> faces; };
namespace atopo {
    template<> struct CellTraits<LegacyVertex> { static constexpr int dimension = 0; using GeometryType = MyPoint3D; static GeometryType get_geometry(const LegacyVertex& cell) { return cell.geometry; } };
    template<> struct CellTraits<LegacyEdge> { static constexpr int dimension = 1; using GeometryType = MyCurve; static GeometryType get_geometry(const LegacyEdge& cell) { return cell.geometry; } };
    template<> struct CellTraits<LegacyFace> { static constexpr int dimension = 2; using GeometryType = MySurface; static GeometryType get_geometry(const LegacyFace& cell) { return cell.geometry; } };
    template<> struct TopologySourceTraits<LegacyMesh> {
        static CellComplex build(const LegacyMesh& mesh);
    };
}

#endif // ATOPO_H
EOF

# Overwrite the atopo.cpp file to match the new header declarations
cat << 'EOF' > src/atopo.cpp
#include "atopo.h"

namespace atopo {
    void CellComplex::setGeneralIncidenceMap(int from_dim, int to_dim, IncidenceMatrix<IncidenceCoeff>&& map) {
        int low_dim = std::min(from_dim, to_dim);
        int high_dim = std::max(from_dim, to_dim);
        DimensionPair key = {low_dim, high_dim};
        if (from_dim < to_dim) { m_general_incidence_maps[key] = std::move(map); }
        else { m_general_incidence_maps[key] = map.transpose(); }
    }

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
        int low_dim = std::min(from_dim, to_dim);
        int high_dim = std::max(from_dim, to_dim);
        DimensionPair key = {low_dim, high_dim};
        auto it = m_general_incidence_maps.find(key);
        if (it == m_general_incidence_maps.end()) {
             return IncidenceMatrix<IncidenceCoeff>(getNumberOfCells(to_dim), getNumberOfCells(from_dim));
        }
        return (from_dim < to_dim) ? it->second : it->second.transpose();
    }

    std::map<int, HomologyGroup> CellComplex::computeHomology() const {
        std::map<int, HomologyGroup> homology;
        int max_dim = 0;
        if(!m_cell_counts.empty()) { max_dim = m_cell_counts.rbegin()->first; }

        for (int p = 0; p <= max_dim + 1; ++p) {
            if (getNumberOfCells(p) == 0 && getNumberOfCells(p-1) == 0 && p > max_dim) continue;

            auto factors_dp = (p > 0) ? detail::compute_invariant_factors(getIncidenceMap(p, p - 1)) : std::vector<long>();
            auto factors_dp1 = detail::compute_invariant_factors(getIncidenceMap(p + 1, p));
            
            long rank_Cp = getNumberOfCells(p);
            long rank_dp = factors_dp.size();
            long rank_dp1 = factors_dp1.size();

            HomologyGroup hg;
            hg.rank = rank_Cp - rank_dp - rank_dp1;
            hg.torsion_coeffs = factors_dp;
            
            if (hg.rank != 0 || !hg.torsion_coeffs.empty() || rank_Cp > 0) {
                homology[p] = hg;
            }
        }
        return homology;
    }

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
    
    CellComplex TopologySourceTraits<LegacyMesh>::build(const LegacyMesh& mesh) {
        CellComplex complex;
        complex.setCellCount(0, mesh.vertices.size()); complex.setCellCount(1, mesh.edges.size()); complex.setCellCount(2, mesh.faces.size());
        IncidenceMatrix<IncidenceCoeff> d1(mesh.vertices.size(), mesh.edges.size());
        for (const auto& edge : mesh.edges) { d1.insert(edge.v_start, edge.id) = -1; d1.insert(edge.v_end,   edge.id) = 1; }
        complex.setBoundaryOperator(1, std::move(d1));
        IncidenceMatrix<IncidenceCoeff> d2(mesh.edges.size(), mesh.faces.size());
        for (const auto& face : mesh.faces) {
            if (face.edge_loop.empty()) continue;
            std::vector<size_t> vertex_loop;
            const auto& first_edge = mesh.edges.at(face.edge_loop.front());
            vertex_loop.push_back(first_edge.v_start);
            size_t last_vertex = first_edge.v_start;
            for (size_t edge_id : face.edge_loop) {
                const auto& edge = mesh.edges.at(edge_id);
                if (edge.v_start == last_vertex) { last_vertex = edge.v_end; } else { last_vertex = edge.v_start; }
                if (last_vertex != vertex_loop.front()) { vertex_loop.push_back(last_vertex); }
            }
            for (size_t i = 0; i < face.edge_loop.size(); ++i) {
                size_t edge_id = face.edge_loop[i];
                const auto& edge = mesh.edges.at(edge_id);
                size_t v_start_of_traversal = vertex_loop[i];
                size_t v_end_of_traversal = vertex_loop[(i + 1) % vertex_loop.size()];
                IncidenceCoeff orientation = (edge.v_start == v_start_of_traversal && edge.v_end == v_end_of_traversal) ? 1 : -1;
                d2.insert(edge_id, face.id) = orientation;
            }
        }
        complex.setBoundaryOperator(2, std::move(d2));
        return complex;
    }

    // Explicit template instantiations to ensure functions are compiled
    template Chain<int> boundary<int>(const CellComplex&, const Chain<int>&);
    template Cochain<int> coboundary<int>(const CellComplex&, const Cochain<int>&);
    template CellComplex create_complex_from_source<LegacyMesh>(const LegacyMesh&);

} // namespace atopo
EOF

echo "✅ File inc/atopo.h has been restored."
echo "   Please clean and rebuild."
echo "   Run: rm -rf build && export PKG_CONFIG_PATH=/opt/lib/pkgconfig:\$PKG_CONFIG_PATH && cmake -S . -B build && cmake --build build"