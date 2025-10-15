#!/bin/bash
#
# This script applies the definitive fix to the LinBox integration, based on
# the ground truth from the compiler error messages. It uses the correct
# 3-argument constructor for the default LinBox sparse matrix.
#
# USAGE: ./final_fix.sh from the project root.

set -e

echo "🚀 Applying final, evidence-based API corrections for LinBox::SparseMatrix..."

# Overwrite the implementation file with the final, correct version.
cat << 'EOF' > src/atopo_linbox.cpp
/**
 * @file atopo_linbox.cpp
 * @brief Isolates the LinBox dependency and implementation details.
 */

// --- Third-Party Includes (CORRECT ORDER) ---
#include <gmp++/gmp++.h>
#include <gmp++/gmp++_int.h>
#include <givaro/zring.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/smith-form.h>

#include "atopo.h" // Our header, which includes Eigen, comes last.

namespace atopo::detail {

    /**
     * @brief Computes the Smith Normal Form using a native LinBox sparse matrix.
     */
    HomologyGroup compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat) {
        if (sparse_mat.nonZeros() == 0) return {};

        using Integer = Givaro::Integer;
        using IntRing = Givaro::ZRing<Integer>;
        using LinBoxSparseMatrix = LinBox::SparseMatrix<IntRing>;

        IntRing Z;
        
        // --- CORRECTED CONSTRUCTOR ---
        // The compiler error showed that the correct constructor for the default
        // sparse matrix type takes exactly 3 arguments.
        LinBoxSparseMatrix A(Z, sparse_mat.rows(), sparse_mat.cols());

        // Copy non-zero entries using the standard 'setEntry' method.
        for (int k=0; k < sparse_mat.outerSize(); ++k) {
            for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(sparse_mat,k); it; ++it) {
                A.setEntry(it.row(), it.col(), Integer(it.value()));
            }
        }
        
        // Call the high-performance Smith Normal Form solver.
        LinBox::SmithList<IntRing> snf_result;
        smithForm(snf_result, A);
        
        // Process the results.
        HomologyGroup result;
        Integer zero(0);
        Integer one(1);

        for(const auto& factor_pair : snf_result) {
            if (factor_pair.first != zero) { result.rank += factor_pair.second; }
            if (factor_pair.first != zero && factor_pair.first != one) {
                for(size_t i = 0; i < factor_pair.second; ++i) {
                    result.torsion_coeffs.push_back(static_cast<long>(factor_pair.first));
                }
            }
        }
        return result;
    }

} // namespace atopo::detail
EOF

echo "✅ Implementation corrected based on compiler feedback."
echo "   This is the final fix. Please clean, rebuild, and run the tests."
echo "   Run: rm -rf build && export PKG_CONFIG_PATH=/opt/lib/pkgconfig:\$PKG_CONFIG_PATH && cmake -S . -B build && cmake --build build && ctest --test-dir build"