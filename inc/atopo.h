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

    /// @brief A template alias for the coefficient type used in chains and matrices.
    /// This can be an integer type (like int for homology over integers) or a floating-point type.
    template<typename T>
    using Coefficient = T;

    /// @brief The primary data structure for representing boundary maps (∂).
    /// It's a sparse matrix from the Eigen library for memory and computational efficiency.
    template<typename T>
    using IncidenceMatrix = Eigen::SparseMatrix<Coefficient<T>>;

    /// @brief A pair of integers representing the source and target dimension for an incidence map.
    using DimensionPair = std::pair<int, int>;


    // --- TRAITS (Compile-Time Polymorphism via Specialization) ---

    /**
     * @brief Trait class to adapt a user-defined type to be a "cell".
     * @details The user must specialize this template for their own geometry types
     * (e.g., MyVertex, MyEdge) to inform the framework of the dimension of each cell type.
     */
    template<typename CellType>
    struct CellTraits { /* User specializes: static constexpr int dimension = D; */ };

    /**
     * @brief Trait class to adapt a user's data source (e.g., a mesh) to build a CellComplex.
     * @details The user must specialize this template for their own mesh or topology data structure.
     * The specialization must provide a static `build` function that constructs and returns a
     * `CellComplex` from the user's data source.
     */
    template<typename TopologySource>
    struct TopologySourceTraits { /* User specializes: static CellComplex<T> build(const Source&); */ };


    // --- CORE FRAMEWORK CLASSES & FUNCTIONS ---

    /**
     * @struct ChainBase
     * @brief A non-templated base struct for `Chain` and `Cochain` to reduce code duplication.
     * @details It holds the essential data: the dimension of the cells and a sparse vector
     * representing a linear combination of those cells.
     */
    template<typename T>
    struct ChainBase {
        int dimension;
        Eigen::SparseVector<Coefficient<T>> data;

        ChainBase(int p, const Eigen::SparseVector<Coefficient<T>>& vec) : dimension(p), data(vec) {}
        ChainBase(int p, Eigen::SparseVector<Coefficient<T>>&& vec) : dimension(p), data(std::move(vec)) {}
    };

    /**
     * @struct Chain
     * @brief Represents a p-chain, a formal linear combination of p-cells.
     * @details For example, a 1-chain could be `2*edge1 - 1*edge2`, representing a path.
     * A 2-chain could be `1*face1 + 1*face2`, representing a surface region.
     */
    template<typename T>
    struct Chain : public ChainBase<T> {
        using ChainBase<T>::ChainBase;
    };
    
    /**
     * @struct Cochain
     * @brief Represents a p-cochain, a linear function that maps p-chains to a coefficient value.
     * @details Cochains are the "dual" concept to chains. For example, a 1-cochain could represent
     * the voltage potential on each 1-cell (edge). Applying it to a 1-chain (a path) would
     * calculate the total voltage drop along that path.
     */
    template<typename T>
    struct Cochain : public ChainBase<T> {
        using ChainBase<T>::ChainBase;
    };

    /**
     * @class CellComplex
     * @brief The central data structure representing a topological space.
     * @details It stores the topology of a space by maintaining a count of cells for each dimension
     * and the incidence matrices that define how they connect. It uses a hybrid storage model:
     * an optimized map for standard boundary operators (p -> p-1) and a general map for all
     * other incidence relations, employing a canonical key to save memory and ensure consistency.
     */
    template<typename T>
    class CellComplex {
    private:
        /// @brief Master list of cell counts, keyed by dimension.
        std::map<int, size_t> m_cell_counts;
        
        /// @brief Optimized storage for the common case: standard boundary operators, ∂_p: C_p -> C_{p-1}.
        /// The map is keyed by the source dimension `p`.
        std::map<int, IncidenceMatrix<T>> m_boundary_maps;
        
        /// @brief General storage for uncommon, non-adjacent incidence relations (e.g., ∂_{3,0}).
        /// Uses a canonical key (min_dim, max_dim) to store each matrix only once.
        std::map<DimensionPair, IncidenceMatrix<T>> m_general_incidence_maps;

    public:
        /// @brief Sets the number of cells for a given dimension.
        void setCellCount(int dim, size_t count) {
            m_cell_counts[dim] = count;
        }

        /// @brief Sets the standard boundary operator matrix (∂_p) for a given source dimension.
        void setBoundaryOperator(int source_dim, IncidenceMatrix<T>&& map) {
            m_boundary_maps[source_dim] = std::move(map);
        }

        /**
         * @brief Sets a general (potentially non-adjacent) incidence map.
         * @details This method automatically handles the canonical keying strategy. It stores
         * the matrix under the key `(min(from, to), max(from, to))`, transposing the input
         * matrix if necessary to maintain consistency. This ensures ∂_{i,j} = (∂_{j,i})^T.
         */
        void setGeneralIncidenceMap(int from_dim, int to_dim, IncidenceMatrix<T>&& map) {
            int low_dim = std::min(from_dim, to_dim);
            int high_dim = std::max(from_dim, to_dim);
            DimensionPair key = {low_dim, high_dim};

            if (from_dim < to_dim) {
                m_general_incidence_maps[key] = std::move(map);
            } else {
                m_general_incidence_maps[key] = map.transpose();
            }
        }

        /// @brief Retrieves the number of cells for a given dimension. Returns 0 if dimension is not present.
        [[nodiscard]] size_t getNumberOfCells(int dim) const {
            auto it = m_cell_counts.find(dim);
            return (it == m_cell_counts.end()) ? 0 : it->second;
        }
        
        /**
         * @brief A unified accessor to get any incidence map between two dimensions.
         * @details It intelligently searches for the requested map.
         * 1. Checks the optimized storage for standard boundary (p -> p-1) or coboundary (p -> p+1) maps.
         * 2. If not found, it checks the general storage using the canonical key.
         * 3. It automatically returns the transpose of a stored matrix if the request direction is reversed.
         * @return The requested incidence matrix. Returns an empty matrix if no relationship is defined.
         */
        [[nodiscard]] IncidenceMatrix<T> getIncidenceMap(int from_dim, int to_dim) const {
            if (to_dim == from_dim - 1) { // Standard boundary operator (∂_p)
                auto it = m_boundary_maps.find(from_dim);
                if (it != m_boundary_maps.end()) return it->second;
            }
            if (to_dim == from_dim + 1) { // Standard coboundary (transpose of next boundary, δ^p = (∂_{p+1})^T)
                auto it = m_boundary_maps.find(to_dim);
                if (it != m_boundary_maps.end()) return it->second.transpose();
            }

            // Fallback to general incidence storage
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

    /**
     * @brief Computes the boundary (∂) of a p-chain.
     * @details The boundary of a p-chain is a (p-1)-chain. A key property in topology
     * is that the "boundary of a boundary is zero" (∂ ∘ ∂ = 0).
     */
    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex<T>& complex, const Chain<T>& chain) {
        if (chain.dimension <= 0) return Chain<T>(-1, Eigen::SparseVector<T>());
        const IncidenceMatrix<T>& d_p = complex.getIncidenceMap(chain.dimension, chain.dimension - 1);
        Eigen::SparseVector<T> boundary_vector = d_p * chain.data;
        boundary_vector.prune(0, 0);
        return Chain<T>(chain.dimension - 1, std::move(boundary_vector));
    }

    /**
     * @brief Computes the coboundary (δ) of a p-cochain.
     * @details The coboundary of a p-cochain is a (p+1)-cochain. The coboundary operator
     * is the dual of the boundary operator. It also satisfies "coboundary of a coboundary is zero" (δ ∘ δ = 0).
     */
    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex<T>& complex, const Cochain<T>& cochain) {
        int p = cochain.dimension;
        // The coboundary δ^p is the transpose of the boundary ∂_{p+1}.
        // Our getIncidenceMap(p, p+1) correctly returns (∂_{p+1})^T.
        const IncidenceMatrix<T> delta_p = complex.getIncidenceMap(p, p + 1);
        Eigen::SparseVector<T> coboundary_vector = delta_p * cochain.data;
        coboundary_vector.prune(0, 0);
        return Cochain<T>(p + 1, std::move(coboundary_vector));
    }

    /// @brief A factory function to build a CellComplex from a user-defined source using the trait system.
    template<typename T, typename TopologySource>
    [[nodiscard]] CellComplex<T> create_complex_from_source(const TopologySource& source) {
        return TopologySourceTraits<TopologySource>::template build<T>(source);
    }

} // namespace atopo


// --- EXAMPLE: USER'S LEGACY CODE ---
// This section demonstrates how a user might have their own existing data structures for geometry.
// The framework is designed to adapt to these without requiring them to be rewritten.

struct LegacyVertex { size_t id; };
struct LegacyEdge { size_t id; size_t v_start, v_end; };
struct LegacyFace { size_t id; std::vector<size_t> edge_loop; };
struct LegacyMesh {
    std::vector<LegacyVertex> vertices;
    std::vector<LegacyEdge> edges;
    std::vector<LegacyFace> faces;
};


// --- EXAMPLE: GLUE CODE ---
// This section shows how the user "glues" their legacy code to the `atopo` framework
// by providing specializations for the required traits.

namespace atopo {
    // Specialization of CellTraits to define the dimension of each legacy type.
    template<> struct CellTraits<LegacyVertex> { static constexpr int dimension = 0; };
    template<> struct CellTraits<LegacyEdge>   { static constexpr int dimension = 1; };
    template<> struct CellTraits<LegacyFace>   { static constexpr int dimension = 2; };

    // Specialization of TopologySourceTraits to define how to build a CellComplex from a LegacyMesh.
    template<> struct TopologySourceTraits<LegacyMesh> {
        template<typename T>
        static CellComplex<T> build(const LegacyMesh& mesh) {
            CellComplex<T> complex;
            complex.setCellCount(0, mesh.vertices.size());
            complex.setCellCount(1, mesh.edges.size());
            complex.setCellCount(2, mesh.faces.size());

            // Define the boundary operator ∂_1: C_1 -> C_0 (edges to vertices)
            IncidenceMatrix<T> d1(mesh.vertices.size(), mesh.edges.size());
            for (const auto& edge : mesh.edges) {
                d1.insert(edge.v_start, edge.id) = -1; // From start vertex
                d1.insert(edge.v_end,   edge.id) = 1;  // To end vertex
            }
            complex.setBoundaryOperator(1, std::move(d1));

            // Define the boundary operator ∂_2: C_2 -> C_1 (faces to edges)
            IncidenceMatrix<T> d2(mesh.edges.size(), mesh.faces.size());
            for (const auto& face : mesh.faces) {
                for (size_t edge_id : face.edge_loop) {
                    // NOTE: A robust implementation would need orientation info from the mesh
                    // to set the correct sign (+1 or -1). We assume +1 for this example.
                    d2.insert(edge_id, face.id) = 1; 
                }
            }
            complex.setBoundaryOperator(2, std::move(d2));

            return complex;
        }
    };
} // namespace atopo

#endif // ATOPO_H