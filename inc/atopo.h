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
    template<typename T> using Coefficient = T;
    using DimensionPair = std::pair<int, int>;
    using IncidenceCoeff = signed short;
    template<typename T> using IncidenceMatrix = Eigen::SparseMatrix<T>;
    struct HomologyGroup {
        int rank = 0;
        std::vector<long> torsion_coeffs;
    };
    namespace detail {
        struct SNFResult {
            size_t rank = 0;
            std::vector<long> torsion_coeffs;
        };
        SNFResult compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat);
    }
    template<typename CellType>
    struct CellTraits {
        static constexpr int dimension = -1;
        using GeometryType = void;
        static GeometryType get_geometry(const CellType& cell);
    };
    template<typename TopologySource>
    struct TopologySourceTraits { /* User specializes */ };
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
        void setBoundaryOperator(int source_dim, IncidenceMatrix<IncidenceCoeff>&& map);
        void setGeneralIncidenceMap(int from_dim, int to_dim, IncidenceMatrix<IncidenceCoeff>&& map);
        [[nodiscard]] size_t getNumberOfCells(int dim) const;
        [[nodiscard]] IncidenceMatrix<IncidenceCoeff> getIncidenceMap(int from_dim, int to_dim) const;
        [[nodiscard]] std::map<int, HomologyGroup> computeHomology() const;
    };
    template<typename T> [[nodiscard]] Chain<T> boundary(const CellComplex& complex, const Chain<T>& chain);
    template<typename T> [[nodiscard]] Cochain<T> coboundary(const CellComplex& complex, const Cochain<T>& cochain);
    template<typename TopologySource> [[nodiscard]] CellComplex create_complex_from_source(const TopologySource& source);
}
struct MyPoint3D { double x, y, z; }; struct MyCurve   { /* ... */ }; struct MySurface { /* ... */ };
struct LegacyVertex { size_t id; MyPoint3D geometry; }; struct LegacyEdge   { size_t id; size_t v_start, v_end; MyCurve geometry; };
struct LegacyFace   { size_t id; std::vector<size_t> edge_loop; MySurface geometry; };
struct LegacyMesh { std::vector<LegacyVertex> vertices; std::vector<LegacyEdge> edges; std::vector<LegacyFace> faces; };
namespace atopo {
    template<> struct CellTraits<LegacyVertex> { static constexpr int dimension = 0; using GeometryType = MyPoint3D; static GeometryType get_geometry(const LegacyVertex& cell) { return cell.geometry; } };
    template<> struct CellTraits<LegacyEdge> { static constexpr int dimension = 1; using GeometryType = MyCurve; static GeometryType get_geometry(const LegacyEdge& cell) { return cell.geometry; } };
    template<> struct CellTraits<LegacyFace> { static constexpr int dimension = 2; using GeometryType = MySurface; static GeometryType get_geometry(const LegacyFace& cell) { return cell.geometry; } };
    template<> struct TopologySourceTraits<LegacyMesh> { static CellComplex build(const LegacyMesh& mesh); };
}
#endif // ATOPO_H
