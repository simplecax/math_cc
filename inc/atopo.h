#ifndef AT_OPO_H
#define AT_OPO_H

#include <iostream>
#include <vector>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <map>
#include <utility>
#include <algorithm>

// Eigen is a header-only library for linear algebra.
#include <Eigen/Sparse>

/**
 * @namespace atopo (Algebraic Topology)
 * @brief The core namespace for the generic topology framework.
 */
namespace atopo {

    // --- Core Data Types ---
    template<typename T>
    using Coefficient = T;
    template<typename T>
    using IncidenceMatrix = Eigen::SparseMatrix<Coefficient<T>>;
    using DimensionPair = std::pair<int, int>;

    // --- TRAITS (Interfaces for User Specialization) ---
    template<typename CellType>
    struct CellTraits { /* User specializes */ };
    template<typename TopologySource>
    struct TopologySourceTraits { /* User specializes */ };

    // --- CORE FRAMEWORK CLASSES & FUNCTIONS ---

    /**
     * @struct ChainBase
     * @brief A common base struct for Chains and Cochains to avoid code duplication.
     */
    template<typename T>
    struct ChainBase {
        int dimension;
        Eigen::SparseVector<Coefficient<T>> data;

        // Constructor
        ChainBase(int p, const Eigen::SparseVector<Coefficient<T>>& vec) : dimension(p), data(vec) {}
        ChainBase(int p, Eigen::SparseVector<Coefficient<T>>&& vec) : dimension(p), data(std::move(vec)) {}
    };

    /**
     * @struct Chain
     * @brief Represents a p-chain. Inherits its structure from ChainBase.
     */
    template<typename T>
    struct Chain : public ChainBase<T> {
        using ChainBase<T>::ChainBase;
    };
    
    /**
     * @struct Cochain
     * @brief Represents a p-cochain. Inherits its structure from ChainBase.
     */
    template<typename T>
    struct Cochain : public ChainBase<T> {
        using ChainBase<T>::ChainBase;
    };

    /**
     * @class CellComplex
     * @brief A robust and efficient data structure for a general cell complex.
     * Uses a hybrid storage model optimized for standard boundary operators
     * while supporting arbitrary incidence relations via canonical keying.
     */
    template<typename T>
    class CellComplex {
    private:
        std::map<int, size_t> m_cell_counts;
        std::map<int, IncidenceMatrix<T>> m_boundary_maps; // Optimized for d_p: C_p -> C_{p-1}
        std::map<DimensionPair, IncidenceMatrix<T>> m_general_incidence_maps; // For all other d_{i,j}

    public:
        void setCellCount(int dim, size_t count) {
            m_cell_counts[dim] = count;
        }

        void setBoundaryOperator(int source_dim, IncidenceMatrix<T>&& map) {
            m_boundary_maps[source_dim] = std::move(map);
        }

        void setGeneralIncidenceMap(int from_dim, int to_dim, IncidenceMatrix<T>&& map) {
            int low_dim = std::min(from_dim, to_dim);
            int high_dim = std::max(from_dim, to_dim);
            DimensionPair key = {low_dim, high_dim};

            if (from_dim < to_dim) {
                m_general_incidence_maps[key] = std::move(map);
            } else {
                // Store the transpose to match the canonical key order
                m_general_incidence_maps[key] = map.transpose();
            }
        }

        [[nodiscard]] size_t getNumberOfCells(int dim) const {
            auto it = m_cell_counts.find(dim);
            if (it == m_cell_counts.end()) {
                return 0;
            }
            return it->second;
        }
        
        [[nodiscard]] IncidenceMatrix<T> getIncidenceMap(int from_dim, int to_dim) const {
            // Case 1: Standard boundary operator (p -> p-1)
            if (to_dim == from_dim - 1) {
                auto it = m_boundary_maps.find(from_dim);
                if (it != m_boundary_maps.end()) return it->second;
            }
            
            // Case 2: Standard coboundary operator (p -> p+1)
            if (to_dim == from_dim + 1) {
                auto it = m_boundary_maps.find(to_dim); // Coboundary is transpose of next boundary
                if (it != m_boundary_maps.end()) return it->second.transpose();
            }

            // Case 3: General incidence, using canonical key
            int low_dim = std::min(from_dim, to_dim);
            int high_dim = std::max(from_dim, to_dim);
            DimensionPair key = {low_dim, high_dim};

            auto it = m_general_incidence_maps.find(key);
            if (it == m_general_incidence_maps.end()) {
                 // Return empty matrix if no map exists
                 return IncidenceMatrix<T>(getNumberOfCells(to_dim), getNumberOfCells(from_dim));
            }
            
            if (from_dim < to_dim) {
                return it->second; // Stored in correct order
            } else {
                return it->second.transpose(); // Stored in reverse order
            }
        }
    };

    /**
     * @brief Computes the boundary (∂) of a p-chain.
     */
    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex<T>& complex, const Chain<T>& chain) {
        if (chain.dimension <= 0) return Chain<T>(-1, Eigen::SparseVector<T>());

        const IncidenceMatrix<T> d_p = complex.getIncidenceMap(chain.dimension, chain.dimension - 1);
        Eigen::SparseVector<T> boundary_vector = d_p * chain.data;
        boundary_vector.prune(0, 0);
        
        return Chain<T>(chain.dimension - 1, std::move(boundary_vector));
    }

    /**
     * @brief Computes the coboundary (δ) of a p-cochain.
     */
    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex<T>& complex, const Cochain<T>& cochain) {
        int p = cochain.dimension;
        const IncidenceMatrix<T> d_p_plus_1 = complex.getIncidenceMap(p + 1, p);
        Eigen::SparseVector<T> coboundary_vector = d_p_plus_1.transpose() * cochain.data;
        coboundary_vector.prune(0, 0);
        
        return Cochain<T>(p + 1, std::move(coboundary_vector));
    }

    template<typename T, typename TopologySource>
    [[nodiscard]] CellComplex<T> create_complex_from_source(const TopologySource& source) {
        return TopologySourceTraits<TopologySource>::template build<T>(source);
    }

} // namespace atopo


// --- USER'S LEGACY CODE ---
struct LegacyVertex { size_t id; };
struct LegacyEdge { size_t id; size_t v_start, v_end; };
struct LegacyFace { size_t id; std::vector<size_t> edge_loop; };
struct LegacyMesh {
    std::vector<LegacyVertex> vertices;
    std::vector<LegacyEdge> edges;
    std::vector<LegacyFace> faces;
};


// --- GLUE CODE ---
namespace atopo {
    template<> struct CellTraits<LegacyVertex> { static constexpr int dimension = 0; };
    template<> struct CellTraits<LegacyEdge> { static constexpr int dimension = 1; };
    template<> struct CellTraits<LegacyFace> { static constexpr int dimension = 2; };

    template<> struct TopologySourceTraits<LegacyMesh> {
        template<typename T>
        static CellComplex<T> build(const LegacyMesh& mesh) {
            CellComplex<T> complex;
            complex.setCellCount(0, mesh.vertices.size());
            complex.setCellCount(1, mesh.edges.size());
            complex.setCellCount(2, mesh.faces.size());

            // Boundary d_1: C_1 -> C_0
            IncidenceMatrix<T> d1(mesh.vertices.size(), mesh.edges.size());
            for (const auto& edge : mesh.edges) {
                d1.insert(edge.v_start, edge.id) = -1;
                d1.insert(edge.v_end,   edge.id) = 1;
            }
            complex.setBoundaryOperator(1, std::move(d1));

            // Boundary d_2: C_2 -> C_1
            IncidenceMatrix<T> d2(mesh.edges.size(), mesh.faces.size());
            for (const auto& face : mesh.faces) {
                for (size_t edge_id : face.edge_loop) {
                    // This is a simplification; correct orientation requires more info
                    d2.insert(edge_id, face.id) = 1; 
                }
            }
            complex.setBoundaryOperator(2, std::move(d2));

            return complex;
        }
    };
} // namespace atopo

#endif // AT_OPO_H