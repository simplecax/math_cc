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

// --- Third-Party Includes (CORRECTED ORDER AND PATHS) ---
// LinBox/Givaro must come BEFORE Eigen to resolve macro conflicts.
#include <gmp++/gmp++.h>
#include <gmp++/gmp++_int.h>
#include <givaro/zring.h> // Correct header for ZRing
#include <linbox/matrix/dense-matrix.h>
#include <linbox/solutions/smith-form.h>

#include <Eigen/Sparse>

namespace atopo {

    // --- Core Data Types ---
    template<typename T> using Coefficient = T;
    using DimensionPair = std::pair<int, int>;
    using IncidenceCoeff = signed short;
    template<typename T> using IncidenceMatrix = Eigen::SparseMatrix<T>;

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

        inline size_t compute_rank_with_linbox(const IncidenceMatrix<IncidenceCoeff>& sparse_mat) {
            if (sparse_mat.nonZeros() == 0) return 0;

            using Integer = Givaro::Integer;
            using IntRing = Givaro::ZRing<Integer>;
            using LinBoxMatrix = LinBox::BlasMatrix<IntRing>;

            IntRing Z;
            LinBoxMatrix A(Z, sparse_mat.rows(), sparse_mat.cols());

            // Copy data from Eigen::SparseMatrix to LinBox::BlasMatrix.
            // This is more robust than a custom wrapper.
            for (int k=0; k < sparse_mat.outerSize(); ++k) {
                for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(sparse_mat,k); it; ++it) {
                    A.setEntry(it.row(), it.col(), Integer(it.value()));
                }
            }
            
            LinBox::SmithList<IntRing> invariant_factors;
            LinBox::smithForm(invariant_factors, A);
            
            // The rank is the number of non-one invariant factors.
            size_t rank = 0;
            Integer one(1); // Create an Integer object for comparison.
            for(const auto& factor_pair : invariant_factors) {
                if (factor_pair.first != one) {
                    rank += factor_pair.second; // Add multiplicity
                }
            }
            return rank;
        }
    } // namespace detail

    // The rest of the file is unchanged, but included for completeness
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
        void setGeneralIncidenceMap(int from_dim, int to_dim, IncidenceMatrix<IncidenceCoeff>&& map) {
            int low_dim = std::min(from_dim, to_dim);
            int high_dim = std::max(from_dim, to_dim);
            DimensionPair key = {low_dim, high_dim};
            if (from_dim < to_dim) { m_general_incidence_maps[key] = std::move(map); }
            else { m_general_incidence_maps[key] = map.transpose(); }
        }
        [[nodiscard]] size_t getNumberOfCells(int dim) const {
            auto it = m_cell_counts.find(dim);
            return (it == m_cell_counts.end()) ? 0 : it->second;
        }
        [[nodiscard]] IncidenceMatrix<IncidenceCoeff> getIncidenceMap(int from_dim, int to_dim) const {
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
        [[nodiscard]] std::map<int, int> computeBettiNumbers() const {
            std::map<int, int> betti_numbers;
            int max_dim = 0;
            if(!m_cell_counts.empty()) { max_dim = m_cell_counts.rbegin()->first; }
            for (int p = 0; p <= max_dim; ++p) {
                long rank_Cp = getNumberOfCells(p);
                long rank_dp = (p > 0) ? detail::compute_rank_with_linbox(getIncidenceMap(p, p - 1)) : 0;
                long rank_dp1 = detail::compute_rank_with_linbox(getIncidenceMap(p + 1, p));
                int betti_p = rank_Cp - rank_dp - rank_dp1;
                if (betti_p != 0 || rank_Cp > 0) { betti_numbers[p] = betti_p; }
            }
            return betti_numbers;
        }
    };
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
struct MyPoint3D { double x, y, z; }; struct MyCurve   { /* ... */ }; struct MySurface { /* ... */ };
struct LegacyVertex { size_t id; MyPoint3D geometry; }; struct LegacyEdge   { size_t id; size_t v_start, v_end; MyCurve geometry; };
struct LegacyFace   { size_t id; std::vector<size_t> edge_loop; MySurface geometry; };
struct LegacyMesh { std::vector<LegacyVertex> vertices; std::vector<LegacyEdge> edges; std::vector<LegacyFace> faces; };
namespace atopo {
    template<> struct CellTraits<LegacyVertex> { static constexpr int dimension = 0; using GeometryType = MyPoint3D; static GeometryType get_geometry(const LegacyVertex& cell) { return cell.geometry; } };
    template<> struct CellTraits<LegacyEdge> { static constexpr int dimension = 1; using GeometryType = MyCurve; static GeometryType get_geometry(const LegacyEdge& cell) { return cell.geometry; } };
    template<> struct CellTraits<LegacyFace> { static constexpr int dimension = 2; using GeometryType = MySurface; static GeometryType get_geometry(const LegacyFace& cell) { return cell.geometry; } };
    template<> struct TopologySourceTraits<LegacyMesh> {
        static CellComplex build(const LegacyMesh& mesh) {
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
    };
}
#endif // ATOPO_H
