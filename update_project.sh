#!/bin/bash
#
# This script refactors the single atopo.h header into a modular, multi-file
# structure for improved maintainability and compile-time performance.
#
# USAGE: ./refactor_header.sh from the project root.

set -e

echo "🚀 Refactoring atopo.h into a modular header structure..."

# ==============================================================================
# 1. CREATE NEW DIRECTORY STRUCTURE AND FILES
# ==============================================================================
mkdir -p inc/atopo
echo " -> Created directory inc/atopo/"

# --- Create atopo/core.h ---
cat << 'EOF' > inc/atopo/core.h
#ifndef ATOPO_CORE_H
#define ATOPO_CORE_H

#include <vector>
#include <map>
#include <utility>

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
}
#endif // ATOPO_CORE_H
EOF

# --- Create atopo/traits.h ---
cat << 'EOF' > inc/atopo/traits.h
#ifndef ATOPO_TRAITS_H
#define ATOPO_TRAITS_H

namespace atopo {
    template<typename CellType>
    struct CellTraits {
        static constexpr int dimension = -1;
    };

    template<typename TopologySource>
    struct TopologySourceTraits { /* User must specialize */ };
}
#endif // ATOPO_TRAITS_H
EOF

# --- Create atopo/chain.h ---
cat << 'EOF' > inc/atopo/chain.h
#ifndef ATOPO_CHAIN_H
#define ATOPO_CHAIN_H

#include <atopo/core.h>

namespace atopo {
    template<typename T>
    struct ChainBase {
        int dimension;
        Eigen::SparseVector<Coefficient<T>> data;
        ChainBase(int p, const Eigen::SparseVector<Coefficient<T>>& vec) : dimension(p), data(vec) {}
        ChainBase(int p, Eigen::SparseVector<Coefficient<T>>&& vec) : dimension(p), data(std::move(vec)) {}
    };

    template<typename T> struct Chain : public ChainBase<T> { using ChainBase<T>::ChainBase; };
    template<typename T> struct Cochain : public ChainBase<T> { using ChainBase<T>::ChainBase; };
}
#endif // ATOPO_CHAIN_H
EOF

# --- Create atopo/complex.h ---
cat << 'EOF' > inc/atopo/complex.h
#ifndef ATOPO_COMPLEX_H
#define ATOPO_COMPLEX_H

#include <atopo/core.h>

namespace atopo {
    namespace detail {
        struct SNFResult {
            size_t rank = 0;
            std::vector<long> torsion_coeffs;
        };
        SNFResult compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat);
    }

    class CellComplex {
    private:
        std::map<int, size_t> m_cell_counts;
        std::map<int, IncidenceMatrix<IncidenceCoeff>> m_boundary_maps;
    public:
        void setCellCount(int dim, size_t count);
        void setBoundaryOperator(int source_dim, IncidenceMatrix<IncidenceCoeff>&& map);
        [[nodiscard]] size_t getNumberOfCells(int dim) const;
        [[nodiscard]] IncidenceMatrix<IncidenceCoeff> getIncidenceMap(int from_dim, int to_dim) const;
        [[nodiscard]] std::map<int, HomologyGroup> computeHomology() const;
    };
}
#endif // ATOPO_COMPLEX_H
EOF

# --- Create atopo/operators.h ---
cat << 'EOF' > inc/atopo/operators.h
#ifndef ATOPO_OPERATORS_H
#define ATOPO_OPERATORS_H

#include <atopo/complex.h>
#include <atopo/chain.h>
#include <atopo/traits.h>

namespace atopo {
    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex& complex, const Chain<T>& chain);

    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex& complex, const Cochain<T>& cochain);

    template<typename TopologySource>
    [[nodiscard]] CellComplex create_complex_from_source(const TopologySource& source);
}
#endif // ATOPO_OPERATORS_H
EOF

# --- Create test/test_common.h ---
cat << 'EOF' > test/test_common.h
#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include "atopo.h"

// --- Example User Code ---
// This code is an example of how a user would adapt their own data
// structures. It belongs with the tests, not in the library headers.
struct MyPoint3D { double x, y, z; };
struct MyCurve   { /* ... */ };
struct MySurface { /* ... */ };

struct LegacyVertex { size_t id; MyPoint3D geometry; };
struct LegacyEdge   { size_t id; size_t v_start, v_end; MyCurve geometry; };
struct LegacyFace   { size_t id; std::vector<size_t> edge_loop; MySurface geometry; };

struct LegacyMesh {
    std::vector<LegacyVertex> vertices;
    std::vector<LegacyEdge> edges;
    std::vector<LegacyFace> faces;
};

// --- Glue Code for the Example ---
namespace atopo {
    template<> struct CellTraits<LegacyVertex> { static constexpr int dimension = 0; };
    template<> struct CellTraits<LegacyEdge>   { static constexpr int dimension = 1; };
    template<> struct CellTraits<LegacyFace>   { static constexpr int dimension = 2; };

    template<> struct TopologySourceTraits<LegacyMesh> {
        static CellComplex build(const LegacyMesh& mesh);
    };
}

#endif // TEST_COMMON_H
EOF

echo " -> Created new modular header files and test_common.h"

# ==============================================================================
# 2. UPDATE MAIN HEADER AND IMPLEMENTATION FILES
# ==============================================================================

# --- Update main atopo.h to be a convenience header ---
cat << 'EOF' > inc/atopo.h
#ifndef ATOPO_H
#define ATOPO_H

// This is a convenience header that includes the entire atopo framework.
// For faster compile times, users can optionally include only the specific
// headers they need from the 'atopo/' subdirectory.

#include <atopo/core.h>
#include <atopo/traits.h>
#include <atopo/chain.h>
#include <atopo/complex.h>
#include <atopo/operators.h>

#endif // ATOPO_H
EOF

# --- Update src/atopo.cpp to use the main header and move the build logic ---
cat << 'EOF' > src/atopo.cpp
#include "atopo.h"
#include "test_common.h" // For LegacyMesh traits

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
            auto snf_dp = (p > 0) ? detail::compute_snf_results(getIncidenceMap(p, p - 1)) : detail::SNFResult{};
            auto snf_dp1 = detail::compute_snf_results(getIncidenceMap(p + 1, p));
            long rank_Cp = getNumberOfCells(p);
            long rank_dp = snf_dp.rank;
            long rank_dp1 = snf_dp1.rank;
            HomologyGroup hg;
            hg.rank = rank_Cp - rank_dp - rank_dp1;
            hg.torsion_coeffs = snf_dp1.torsion_coeffs;
            if (hg.rank != 0 || !hg.torsion_coeffs.empty() || rank_Cp > 0) { homology[p] = hg; }
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
                } else {
                    d2.insert(edge.id, face.id) = -1;
                    current_v = edge.v_start;
                }
            }
        }
        complex.setBoundaryOperator(2, std::move(d2));
        return complex;
    }

    template Chain<int> boundary<int>(const CellComplex&, const Chain<int>&);
    template Cochain<int> coboundary<int>(const CellComplex&, const Cochain<int>&);
    template CellComplex create_complex_from_source<LegacyMesh>(const LegacyMesh&);
}
EOF
echo " -> Updated src/atopo.cpp"

# ==============================================================================
# 3. UPDATE TEST FILES
# ==============================================================================
# All test files now just need to include gtest and test_common.h
sed -i '1s/^/#include <gtest\/gtest.h>\n#include "test_common.h"\n/' test/test_atopo.cpp
sed -i '/#include "atopo.h"/d' test/test_atopo.cpp

sed -i '1s/^/#include <gtest\/gtest.h>\n#include "test_common.h"\n/' test/test_homology.cpp
sed -i '/#include "atopo.h"/d' test/test_homology.cpp
sed -i '/#include <vector>/d' test/test_homology.cpp

sed -i '1s/^/#include <gtest\/gtest.h>\n#include "test_common.h"\n/' test/test_rp2_integration.cpp
sed -i '/#include "atopo.h"/d' test/test_rp2_integration.cpp
sed -i '/#include <vector>/d' test/test_rp2_integration.cpp
sed -i '/#include <algorithm>/d' test/test_rp2_integration.cpp
sed -i '/#include <numeric>/d' test/test_rp2_integration.cpp
echo " -> Updated test files to use test_common.h"

# ==============================================================================
# 4. UPDATE CMAKE
# ==============================================================================
# The CMakeLists.txt needs to know about the new test include path
sed -i '/target_link_directories(atopo_tests PRIVATE ${LINBOX_LIBRARY_DIRS})/a \
# The test target needs to find test_common.h\
target_include_directories(atopo_tests PRIVATE test)' CMakeLists.txt
sed -i '/target_include_directories(atopo PRIVATE/a \    test # For test_common.h needed by atopo.cpp' CMakeLists.txt
echo " -> Updated CMakeLists.txt"

echo "✅ Header refactoring complete."
echo "   Please clean and rebuild."
echo "   Run: rm -rf build && export PKG_CONFIG_PATH=/opt/lib/pkgconfig:\$PKG_CONFIG_PATH && cmake -S . -B build && cmake --build build && ctest --test-dir build"