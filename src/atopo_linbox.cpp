#include <gmp++/gmp++.h>
#include <gmp++/gmp++_int.h>
#include <givaro/zring.h>
#include <linbox/matrix/dense-matrix.h>
#include <linbox/solutions/smith-form.h>

#include "atopo.h"

namespace atopo::detail {

    std::vector<long> compute_invariant_factors(const IncidenceMatrix<IncidenceCoeff>& sparse_mat) {
        if (sparse_mat.nonZeros() == 0) return {};

        using Integer = Givaro::Integer;
        using IntRing = Givaro::ZRing<Integer>;
        using LinBoxMatrix = LinBox::BlasMatrix<IntRing>;

        IntRing Z;
        LinBoxMatrix A(Z, sparse_mat.rows(), sparse_mat.cols());

        for (int k=0; k < sparse_mat.outerSize(); ++k) {
            for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(sparse_mat,k); it; ++it) {
                A.setEntry(it.row(), it.col(), Integer(it.value()));
            }
        }
        
        LinBox::SmithList<IntRing> snf_result;
        LinBox::smithForm(snf_result, A);
        
        std::vector<long> factors;
        Integer one(1);
        for(const auto& factor_pair : snf_result) {
            if (factor_pair.first != one && factor_pair.first != Integer(0)) {
                for(size_t i = 0; i < factor_pair.second; ++i) {
                    factors.push_back(static_cast<long>(factor_pair.first));
                }
            }
        }
        return factors;
    }

} // namespace atopo::detail
