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
