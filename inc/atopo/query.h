#ifndef ATOPO_QUERY_H
#define ATOPO_QUERY_H

#include <atopo/complex.h>
#include <atopo/persistence.h>
#include <algorithm>

namespace atopo::query {

    /**
     * @brief Returns the p-th Betti number of the complex.
     */
    inline size_t betti_number(const CellComplex& complex, int p) {
        auto homology = complex.computeHomology();
        if (homology.count(p)) {
            return homology.at(p).rank;
        }
        return 0;
    }

    /**
     * @brief Checks if the complex is path-connected (Betti 0 == 1).
     */
    inline bool is_connected(const CellComplex& complex) {
        return betti_number(complex, 0) == 1;
    }

    /**
     * @brief Checks if the complex has any "voids" (Betti 2 > 0).
     */
    inline bool has_void(const CellComplex& complex) {
        return betti_number(complex, 2) > 0;
    }

    /**
     * @brief Returns the number of connected components.
     */
    inline size_t count_connected_components(const CellComplex& complex) {
        return betti_number(complex, 0);
    }

    /**
     * @brief Checks if the complex has any "tunnels" (Betti 1 > 0).
     */
    inline bool has_tunnel(const CellComplex& complex) {
        return betti_number(complex, 1) > 0;
    }

} // namespace atopo::query

#endif // ATOPO_QUERY_H
