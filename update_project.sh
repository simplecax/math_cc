#!/bin/bash
#
# This script refactors the project to unify the SNFResult and HomologyGroup
# structs, with a fix for the sed command.
#
# USAGE: ./unify_structs_fixed.sh from the project root.

set -e

echo "🚀 Unifying HomologyGroup and SNFResult structs..."

# ==============================================================================
# 1. UPDATE inc/atopo/core.h
# ==============================================================================
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

    /**
     * @brief Represents a homology group.
     * @note Its properties (rank and torsion) are derived from the metadata
     * of the Smith Normal Form (SNF) of the boundary matrices.
     */
    struct HomologyGroup {
        size_t rank = 0; //!< The rank of the free part (the Betti number, β).
        std::vector<long> torsion_coeffs; //!< The torsion coefficients.
    };
}
#endif // ATOPO_CORE_H
EOF
echo " -> Updated inc/atopo/core.h"

# ==============================================================================
# 2. UPDATE inc/atopo/complex.h
# ==============================================================================
# Remove the entire detail::SNFResult struct definition.
sed -i '/namespace detail {/,/};/d' inc/atopo/complex.h
# Re-insert the namespace with the updated function declaration.
sed -i '/class CellComplex {/i \
namespace atopo { \
    namespace detail { \
        HomologyGroup compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat); \
    } \
' inc/atopo/complex.h
echo " -> Updated inc/atopo/complex.h"

# ==============================================================================
# 3. UPDATE src/atopo_linbox.cpp
# ==============================================================================
sed -i 's/SNFResult compute_snf_results/HomologyGroup compute_snf_results/' src/atopo_linbox.cpp
sed -i 's/SNFResult result;/HomologyGroup result;/' src/atopo_linbox.cpp
echo " -> Updated src/atopo_linbox.cpp"

# ==============================================================================
# 4. UPDATE src/atopo.cpp (CORRECTED COMMAND)
# ==============================================================================
# This command is now a single line, fixing the "unterminated 's' command" error.
sed -i 's/detail::SNFResult{}/atopo::HomologyGroup{}/' src/atopo.cpp
echo " -> Updated src/atopo.cpp"

echo "✅ Struct unification complete."
echo "   Please clean and rebuild."
echo "   Run: rm -rf build && export PKG_CONFIG_PATH=/opt/lib/pkgconfig:\$PKG_CONFIG_PATH && cmake -S . -B build && cmake --build build && ctest --test-dir build"