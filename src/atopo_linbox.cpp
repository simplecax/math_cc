// --- Third-Party Includes (CORRECTED ORDER) ---
// LinBox/Givaro must come BEFORE our own header to resolve macro conflicts.
#include <gmp++/gmp++.h>
#include <gmp++/gmp++_int.h>
#include <givaro/zring.h>
#include <linbox/matrix/dense-matrix.h>
#include <linbox/solutions/smith-form.h>

#include "atopo.h" // Our header, which includes Eigen, comes last.

namespace atopo::detail {

    size_t compute_rank_with_linbox(const IncidenceMatrix<IncidenceCoeff>& sparse_mat) {
        if (sparse_mat.nonZeros() == 0) return 0;

        using Integer = Givaro::Integer;
        using IntRing = Givaro::ZRing<Integer>;
        using LinBoxMatrix = LinBox::BlasMatrix<IntRing>;

        IntRing Z;
        LinBoxMatrix A(Z, sparse_mat.rows(), sparse_mat.cols());

        // Copy data from Eigen::SparseMatrix to LinBox::BlasMatrix.
        for (int k=0; k < sparse_mat.outerSize(); ++k) {
            for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(sparse_mat,k); it; ++it) {
                A.setEntry(it.row(), it.col(), Integer(it.value()));
            }
        }
        
        LinBox::SmithList<IntRing> invariant_factors;
        LinBox::smithForm(invariant_factors, A);
        
        size_t rank = 0;
        for(const auto& factor_pair : invariant_factors) {
            Integer zero(0);
        if (factor_pair.first != zero) {
                rank += factor_pair.second;
            }
        }
        return rank;
    }

} // namespace atopo::detail
