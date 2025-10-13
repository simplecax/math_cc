#ifndef AT_OPO_H
#define AT_OPO_H

#include <iostream>
#include <vector>
#include <cstddef>
#include <memory>
#include <stdexcept>

// Eigen is a header-only library for linear algebra.
#include <Eigen/Sparse>

/**
 * @namespace atopo (Algebraic Topology)
 * @brief The core namespace for the generic topology framework.
 */
namespace atopo {

    // --- Core Data Types ---
    template<typename T>
    using Coefficient = T; /// The numeric type for chain coefficients (e.g., int, double).

    template<typename T>
    using IncidenceMatrix = Eigen::SparseMatrix<Coefficient<T>>; /// The boundary map matrix.

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
        // Inherit constructors from the base class
        using ChainBase<T>::ChainBase;
    };
    
    /**
     * @struct Cochain
     * @brief Represents a p-cochain. Inherits its structure from ChainBase.
     * It is a distinct type from Chain, ensuring type safety.
     */
    template<typename T>
    struct Cochain : public ChainBase<T> {
        // Inherit constructors from the base class
        using ChainBase<T>::ChainBase;
    };

    /**
     * @class CellComplex
     * @brief A concrete, efficient data structure holding the topological space.
     */
    template<typename T>
    class CellComplex {
    private:
        std::vector<IncidenceMatrix<T>> marrBoundaryMaps;
        std::vector<size_t> marrCellCounts;

    public:
        CellComplex(size_t max_dim) : marrBoundaryMaps(max_dim), marrCellCounts(max_dim + 1, 0) {}
        
        void setBoundaryMap(size_t dim_of_target_chain, IncidenceMatrix<T>&& boundary_map) {
            if (dim_of_target_chain >= marrBoundaryMaps.size()) throw std::out_of_range("Dimension index is out of range.");
            marrBoundaryMaps[dim_of_target_chain] = std::move(boundary_map);
        }

        void setCellCount(size_t dim, size_t count) {
            if (dim >= marrCellCounts.size()) throw std::out_of_range("Dimension index is out of range.");
            marrCellCounts[dim] = count;
        }

        [[nodiscard]] const IncidenceMatrix<T>& getBoundaryMap(size_t dim_of_source_chain) const {
            if (dim_of_source_chain == 0 || dim_of_source_chain > marrBoundaryMaps.size()) throw std::out_of_range("Dimension is out of range for boundary map.");
            return marrBoundaryMaps[dim_of_source_chain - 1];
        }
        
        [[nodiscard]] size_t getNumberOfCells(size_t dim) const {
            if (dim >= marrCellCounts.size()) throw std::out_of_range("Dimension index out of range for cell count.");
            return marrCellCounts[dim];
        }

        [[nodiscard]] size_t getMaxDimension() const {
            return marrBoundaryMaps.size();
        }
    };

    /**
     * @brief Computes the boundary (∂) of a p-chain.
     * @return A (p-1)-chain representing the boundary.
     */
    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex<T>& complex, const Chain<T>& chain) {
        if (chain.dimension <= 0) return Chain<T>(-1, Eigen::SparseVector<T>());

        const IncidenceMatrix<T>& d_p = complex.getBoundaryMap(chain.dimension);
        Eigen::SparseVector<T> boundary_vector = d_p * chain.data;
        boundary_vector.prune(0, 0);
        
        return Chain<T>(chain.dimension - 1, std::move(boundary_vector));
    }

    /**
     * @brief Computes the coboundary (δ) of a p-cochain.
     * @return A (p+1)-cochain.
     */
    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex<T>& complex, const Cochain<T>& cochain) {
        int p = cochain.dimension;
        if (p < 0 || static_cast<size_t>(p) >= complex.getMaxDimension()) {
            return Cochain<T>(p + 1, Eigen::SparseVector<T>());
        }

        const IncidenceMatrix<T>& d_p_plus_1 = complex.getBoundaryMap(p + 1);
        Eigen::SparseVector<T> coboundary_vector = d_p_plus_1.transpose() * cochain.data;
        coboundary_vector.prune(0, 0);
        
        return Cochain<T>(p + 1, std::move(coboundary_vector));
    }

    /**
     * @brief Builder function to construct a `CellComplex` from a legacy source.
     */
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
            CellComplex<T> complex(2);
            complex.setCellCount(0, mesh.vertices.size());
            complex.setCellCount(1, mesh.edges.size());
            complex.setCellCount(2, mesh.faces.size());

            IncidenceMatrix<T> d1(mesh.vertices.size(), mesh.edges.size());
            for (const auto& edge : mesh.edges) {
                d1.insert(edge.v_start, edge.id) = -1; d1.insert(edge.v_end,   edge.id) = 1;
            }
            complex.setBoundaryMap(0, std::move(d1));

            IncidenceMatrix<T> d2(mesh.edges.size(), mesh.faces.size());
            for (const auto& face : mesh.faces) {
                for (size_t edge_id : face.edge_loop) { d2.insert(edge_id, face.id) = 1; }
            }
            complex.setBoundaryMap(1, std::move(d2));

            return complex;
        }
    };
} // namespace atopo


#endif // AT_OPO_H
