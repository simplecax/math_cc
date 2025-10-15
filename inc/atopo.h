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

/**
 * @file atopo.h
 * @brief A trait-based C++ framework for algebraic topology computations.
 *
 * This header provides the core public interface for creating cell complexes,
 * adapting user data structures, and computing homology groups.
 */
namespace atopo {

    // --- Core Data Types ---
    template<typename T> using Coefficient = T;
    using DimensionPair = std::pair<int, int>;
    using IncidenceCoeff = signed short;
    template<typename T> using IncidenceMatrix = Eigen::SparseMatrix<T>;

    /**
     * @brief Represents a homology group, containing its rank and torsion components.
     */
    struct HomologyGroup {
        int rank = 0; //!< The rank of the free part (the Betti number, β).
        std::vector<long> torsion_coeffs; //!< The torsion coefficients (e.g., {2} for Z₂).
    };

    // --- TRAITS ---
    /**
     * @brief Trait class to adapt a user-defined type to be a "cell".
     * Specialize this struct for your own cell types (e.g., Vertex, Edge).
     */
    template<typename CellType>
    struct CellTraits {
        static constexpr int dimension = -1;
    };

    /**
     * @brief Trait class to adapt a user's data source to build a CellComplex.
     * Specialize this struct for your own mesh or topology data structures.
     */
    template<typename TopologySource>
    struct TopologySourceTraits { /* User must specialize the `build` static method */ };

    // --- Implementation Detail ---
    namespace detail {
        struct SNFResult {
            size_t rank = 0;
            std::vector<long> torsion_coeffs;
        };
        // The definition is in atopo_linbox.cpp to isolate the LinBox dependency.
        SNFResult compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat);
    }

    /**
     * @brief Base class for chains and cochains.
     * @tparam T The coefficient type (e.g., int, double).
     */
    template<typename T>
    struct ChainBase {
        int dimension;
        Eigen::SparseVector<Coefficient<T>> data;
        ChainBase(int p, const Eigen::SparseVector<Coefficient<T>>& vec) : dimension(p), data(vec) {}
        ChainBase(int p, Eigen::SparseVector<Coefficient<T>>&& vec) : dimension(p), data(std::move(vec)) {}
    };

    template<typename T> struct Chain : public ChainBase<T> { using ChainBase<T>::ChainBase; };
    template<typename T> struct Cochain : public ChainBase<T> { using ChainBase<T>::ChainBase; };

    /**
     * @brief Represents a CW complex using incidence matrices.
     */
    class CellComplex {
    private:
        std::map<int, size_t> m_cell_counts;
        std::map<int, IncidenceMatrix<IncidenceCoeff>> m_boundary_maps;
        std::map<DimensionPair, IncidenceMatrix<IncidenceCoeff>> m_general_incidence_maps;
    public:
        void setCellCount(int dim, size_t count);
        void setBoundaryOperator(int source_dim, IncidenceMatrix<IncidenceCoeff>&& map);
        [[nodiscard]] size_t getNumberOfCells(int dim) const;
        [[nodiscard]] IncidenceMatrix<IncidenceCoeff> getIncidenceMap(int from_dim, int to_dim) const;
        [[nodiscard]] std::map<int, HomologyGroup> computeHomology() const;
    };
    
    /**
     * @brief Computes the boundary of a chain.
     * @param complex The cell complex defining the boundary operators.
     * @param chain The input p-chain.
     * @return The resulting (p-1)-chain.
     */
    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex& complex, const Chain<T>& chain);

    /**
     * @brief Computes the coboundary of a cochain.
     * @param complex The cell complex defining the boundary operators.
     * @param cochain The input p-cochain.
     * @return The resulting (p+1)-cochain.
     */
    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex& complex, const Cochain<T>& cochain);

    /**
     * @brief Factory function to create a CellComplex from a user data source.
     * @tparam TopologySource The type of the user's data source.
     * @param source The user's data source object.
     * @return A fully constructed CellComplex.
     */
    template<typename TopologySource>
    [[nodiscard]] CellComplex create_complex_from_source(const TopologySource& source);
}

// --- Example Glue Code ---
struct MyPoint3D { double x, y, z; }; struct MyCurve { /*...*/ }; struct MySurface { /*...*/ };
struct LegacyVertex { size_t id; MyPoint3D geometry; }; struct LegacyEdge { size_t id; size_t v_start, v_end; MyCurve geometry; };
struct LegacyFace { size_t id; std::vector<size_t> edge_loop; MySurface geometry; };
struct LegacyMesh { std::vector<LegacyVertex> vertices; std::vector<LegacyEdge> edges; std::vector<LegacyFace> faces; };

namespace atopo {
    template<> struct CellTraits<LegacyVertex> { static constexpr int dimension = 0; };
    template<> struct CellTraits<LegacyEdge> { static constexpr int dimension = 1; };
    template<> struct CellTraits<LegacyFace> { static constexpr int dimension = 2; };
    template<> struct TopologySourceTraits<LegacyMesh> { static CellComplex build(const LegacyMesh& mesh); };
}

#endif // ATOPO_H
