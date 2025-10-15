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
