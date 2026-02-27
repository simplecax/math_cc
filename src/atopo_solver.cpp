#include "atopo.h"
#include <vector>
#include <algorithm>
#include <unordered_map>

namespace atopo::detail {

    /**
     * @brief Computes the rank of a sparse matrix over Z2 (binary field).
     * This replaces the LinBox Smith Normal Form solver.
     * Note: Torsion coefficients will be empty as Z2 does not capture them.
     */
    HomologyGroup compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat) {
        if (sparse_mat.nonZeros() == 0) return {};

        int numCols = static_cast<int>(sparse_mat.cols());
        
        // pivot_to_col[row_index] = col_index
        // Records which column (index) is currently providing the pivot at a given row.
        std::unordered_map<int, int> pivot_to_col;
        
        // Use a vector of sorted row indices to represent each column in Z2.
        // For Z2, addition is XOR, which corresponds to the symmetric difference of index sets.
        using Z2Vector = std::vector<int>;
        std::vector<Z2Vector> reduced_columns(numCols);

        for (int j = 0; j < numCols; ++j) {
            Z2Vector& curr = reduced_columns[j];
            
            // 1. Initialize current column from sparse matrix (project to Z2).
            // We use the parity of the incidence coefficients.
            for (Eigen::SparseMatrix<IncidenceCoeff>::InnerIterator it(sparse_mat, j); it; ++it) {
                if (it.value() % 2 != 0) {
                    curr.push_back(it.row());
                }
            }
            // Eigen's SparseMatrix InnerIterator typically visits in increasing row order for ColMajor.
            // But we sort just to be absolutely safe for set_symmetric_difference.
            std::sort(curr.begin(), curr.end());

            // 2. Reduce the column using Gaussian Elimination (XOR-based).
            while (!curr.empty()) {
                // The pivot is the highest row index (last element in our sorted vector).
                int pivot = curr.back();
                auto it_pivot = pivot_to_col.find(pivot);
                
                if (it_pivot == pivot_to_col.end()) {
                    // This row hasn't been used as a pivot yet.
                    pivot_to_col[pivot] = j;
                    break;
                } else {
                    // Collision: XOR with the existing column that already has this pivot.
                    int other_idx = it_pivot->second;
                    const Z2Vector& other = reduced_columns[other_idx];
                    
                    Z2Vector next_curr;
                    next_curr.reserve(curr.size() + other.size());
                    
                    // Symmetric difference is equivalent to XOR for sets of non-zero indices.
                    std::set_symmetric_difference(
                        curr.begin(), curr.end(),
                        other.begin(), other.end(),
                        std::back_inserter(next_curr)
                    );
                    curr = std::move(next_curr);
                }
            }
        }

        // The rank of the matrix over Z2 is the number of pivots found.
        HomologyGroup result;
        result.rank = pivot_to_col.size();
        result.torsion_coeffs = {}; // Z2 coefficients do not capture torsion (Z coefficients required).
        return result;
    }

} // namespace atopo::detail
