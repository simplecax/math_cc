/**
 * @file atopo_linbox.cpp
 * @brief Isolates the LinBox dependency and implementation details.
 */

// --- Third-Party Includes (CORRECT ORDER) ---
// LinBox/Givaro headers must come BEFORE atopo.h to resolve macro conflicts.
#include <gmp++/gmp++.h>
#include <gmp++/gmp++_int.h>
#include <givaro/zring.h>
#include <linbox/matrix/dense-matrix.h>
#include <linbox/solutions/smith-form.h>

#include "atopo.h" // Our header, which includes Eigen, comes last.

namespace atopo::detail {

    /**
     * @brief Computes the Smith Normal Form of a matrix using LinBox.
     * @param sparse_mat The input sparse incidence matrix from Eigen.
     * @return A struct containing the rank and torsion coefficients.
     */
    HomologyGroup compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat) {
        if (sparse_mat.nonZeros() == 0) return {};

        using Integer = Givaro::Integer;
        using IntRing = Givaro::ZRing<Integer>;
        using LinBoxMatrix = LinBox::BlasMatrix<IntRing>;

        // 1. Create a LinBox matrix and its associated integer ring.
        IntRing Z;
        LinBoxMatrix A(Z, sparse_mat.rows(), sparse_mat.cols());

        // 2. Copy data from Eigen to LinBox.
        for (int k=0; k < sparse_mat.outerSize(); ++k) {
            for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(sparse_mat,k); it; ++it) {
                A.setEntry(it.row(), it.col(), Integer(it.value()));
            }
        }
        
        // 3. Call the high-performance Smith Normal Form solver.
        LinBox::SmithList<IntRing> snf_result;
        LinBox::smithForm(snf_result, A);
        
        // 4. Process the results.
        HomologyGroup result;
        Integer zero(0);
        Integer one(1);

        for(const auto& factor_pair : snf_result) {
            // Rank is the number of non-zero factors.
            if (factor_pair.first != zero) {
                result.rank += factor_pair.second;
            }
            // Torsion coefficients are factors that are not 0 or 1.
            if (factor_pair.first != zero && factor_pair.first != one) {
                for(size_t i = 0; i < factor_pair.second; ++i) {
                    result.torsion_coeffs.push_back(static_cast<long>(factor_pair.first));
                }
            }
        }
        return result;
    }

} // namespace atopo::detail
