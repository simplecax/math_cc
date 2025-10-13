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

// Eigen is a header-only library for linear algebra, used for our sparse matrix and vector operations.
#include <Eigen/Sparse>

/**
 * @file atopo.h
 * @brief A single-header, trait-based C++ framework for algebraic topology.
 *
 * This file provides a set of tools to represent and manipulate topological spaces
 * (cell complexes) using the principles of algebraic topology. It is designed to be
 * generic and adaptable to existing user-defined geometry and mesh data structures
 * through a trait-based specialization system.
 */

/**
 * @namespace atopo
 * @brief The core namespace for the algebraic topology framework.
 */
namespace atopo {

    // --- Core Data Types ---
    template<typename T>
    using Coefficient = T;
    template<typename T>
    using IncidenceMatrix = Eigen::SparseMatrix<Coefficient<T>>;
    using DimensionPair = std::pair<int, int>;


    // --- TRAITS (Compile-Time Polymorphism via Specialization) ---

    /**
     * @brief Trait class to adapt a user-defined type to be a "cell".
     * @details The user must specialize this template for their own geometry types.
     * The specialization must provide the cell's dimension, its associated geometry
     * C++ type, and a static function to access that geometry.
     */
    template<typename CellType>
    struct CellTraits {
        // The dimension of the cell (e.g., 0 for a vertex, 1 for an edge).
        static constexpr int dimension = -1; 

        // The C++ type of the associated geometric entity (e.g., Point3D, NurbsCurve).
        using GeometryType = void;

        // A static function to retrieve the geometric entity from a cell instance.
        static GeometryType get_geometry(const CellType& cell);
    };

    /**
     * @brief Trait class to adapt a user's data source (e.g., a mesh) to build a CellComplex.
     * @details The user must specialize this template for their own mesh or topology data structure.
     * The specialization must provide a static `build` function that constructs and returns a
     * `CellComplex` from the user's data source.
     */
    template<typename TopologySource>
    struct TopologySourceTraits { /* User specializes: static CellComplex<T> build(const Source&); */ };


    // --- CORE FRAMEWORK CLASSES & FUNCTIONS --- (Unchanged from previous version)

    template<typename T>
    struct ChainBase {
        int dimension;
        Eigen::SparseVector<Coefficient<T>> data;
        ChainBase(int p, const Eigen::SparseVector<Coefficient<T>>& vec) : dimension(p), data(vec) {}
        ChainBase(int p, Eigen::SparseVector<Coefficient<T>>&& vec) : dimension(p), data(std::move(vec)) {}
    };

    template<typename T>
    struct Chain : public ChainBase<T> { using ChainBase<T>::ChainBase; };
    
    template<typename T>
    struct Cochain : public ChainBase<T> { using ChainBase<T>::ChainBase; };

    template<typename T>
    class CellComplex {
    private:
        std::map<int, size_t> m_cell_counts;
        std::map<int, IncidenceMatrix<T>> m_boundary_maps;
        std::map<DimensionPair, IncidenceMatrix<T>> m_general_incidence_maps;

    public:
        void setCellCount(int dim, size_t count) { m_cell_counts[dim] = count; }

        void setBoundaryOperator(int source_dim, IncidenceMatrix<T>&& map) {
            m_boundary_maps[source_dim] = std::move(map);
        }

        void setGeneralIncidenceMap(int from_dim, int to_dim, IncidenceMatrix<T>&& map) {
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
        
        [[nodiscard]] IncidenceMatrix<T> getIncidenceMap(int from_dim, int to_dim) const {
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
                 return IncidenceMatrix<T>(getNumberOfCells(to_dim), getNumberOfCells(from_dim));
            }
            return (from_dim < to_dim) ? it->second : it->second.transpose();
        }
    };

    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex<T>& complex, const Chain<T>& chain) {
        if (chain.dimension <= 0) return Chain<T>(-1, Eigen::SparseVector<T>());
        const IncidenceMatrix<T>& d_p = complex.getIncidenceMap(chain.dimension, chain.dimension - 1);
        Eigen::SparseVector<T> boundary_vector = d_p * chain.data;
        boundary_vector.prune(0, 0);
        return Chain<T>(chain.dimension - 1, std::move(boundary_vector));
    }

    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex<T>& complex, const Cochain<T>& cochain) {
        int p = cochain.dimension;
        const IncidenceMatrix<T> delta_p = complex.getIncidenceMap(p, p + 1);
        Eigen::SparseVector<T> coboundary_vector = delta_p * cochain.data;
        coboundary_vector.prune(0, 0);
        return Cochain<T>(p + 1, std::move(coboundary_vector));
    }

    template<typename T, typename TopologySource>
    [[nodiscard]] CellComplex<T> create_complex_from_source(const TopologySource& source) {
        return TopologySourceTraits<TopologySource>::template build<T>(source);
    }

} // namespace atopo


// --- EXAMPLE: USER'S LEGACY CODE ---
// This section demonstrates how a user's geometry types might be defined.
// Note they now include associated geometry (e.g., MyPoint3D).

// Dummy geometry types for the example
struct MyPoint3D { double x, y, z; };
struct MyCurve   { /* e.g., NURBS curve data */ };
struct MySurface { /* e.g., NURBS surface data */ };

struct LegacyVertex { size_t id; MyPoint3D geometry; };
struct LegacyEdge   { size_t id; size_t v_start, v_end; MyCurve geometry; };
struct LegacyFace   { size_t id; std::vector<size_t> edge_loop; MySurface geometry; };

struct LegacyMesh {
    std::vector<LegacyVertex> vertices;
    std::vector<LegacyEdge> edges;
    std::vector<LegacyFace> faces;
};


// --- EXAMPLE: GLUE CODE ---
// This section shows how the user "glues" their legacy code to the `atopo` framework.
// Note how the CellTraits specialization now provides all three required pieces of information.

namespace atopo {

    template<> struct CellTraits<LegacyVertex> {
        static constexpr int dimension = 0;
        using GeometryType = MyPoint3D;
        static GeometryType get_geometry(const LegacyVertex& cell) { return cell.geometry; }
    };
    template<> struct CellTraits<LegacyEdge> {
        static constexpr int dimension = 1;
        using GeometryType = MyCurve;
        static GeometryType get_geometry(const LegacyEdge& cell) { return cell.geometry; }
    };
    template<> struct CellTraits<LegacyFace> {
        static constexpr int dimension = 2;
        using GeometryType = MySurface;
        static GeometryType get_geometry(const LegacyFace& cell) { return cell.geometry; }
    };

    // The TopologySourceTraits specialization remains the same, as it only deals with topology.
    template<> struct TopologySourceTraits<LegacyMesh> {
        template<typename T>
        static CellComplex<T> build(const LegacyMesh& mesh) {
            CellComplex<T> complex;
            complex.setCellCount(0, mesh.vertices.size());
            complex.setCellCount(1, mesh.edges.size());
            complex.setCellCount(2, mesh.faces.size());

            IncidenceMatrix<T> d1(mesh.vertices.size(), mesh.edges.size());
            for (const auto& edge : mesh.edges) {
                d1.insert(edge.v_start, edge.id) = -1;
                d1.insert(edge.v_end,   edge.id) = 1;
            }
            complex.setBoundaryOperator(1, std::move(d1));

            IncidenceMatrix<T> d2(mesh.edges.size(), mesh.faces.size());
            for (const auto& face : mesh.faces) {
                for (size_t edge_id : face.edge_loop) {
                    d2.insert(edge_id, face.id) = 1; 
                }
            }
            complex.setBoundaryOperator(2, std::move(d2));

            return complex;
        }
    };
} // namespace atopo

#endif // ATOPO_H