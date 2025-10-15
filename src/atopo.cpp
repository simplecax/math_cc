#include "atopo.h"

namespace atopo {
    IncidenceMatrix<IncidenceCoeff> CellComplex::getIncidenceMap(int from_dim, int to_dim) const {
        if (to_dim == from_dim - 1) {
            auto it = m_boundary_maps.find(from_dim);
            if (it != m_boundary_maps.end()) return it->second;
        }
        if (to_dim == from_dim + 1) {
            auto it = m_boundary_maps.find(to_dim);
            if (it != m_boundary_maps.end()) return it->second.transpose();
        }
        int low_dim = std::min(from_dim, to_dim);
        int high_dim = std::max(from_dim, to_dim);
        DimensionPair key = {low_dim, high_dim};
        auto it = m_general_incidence_maps.find(key);
        if (it == m_general_incidence_maps.end()) {
             return IncidenceMatrix<IncidenceCoeff>(getNumberOfCells(to_dim), getNumberOfCells(from_dim));
        }
        return (from_dim < to_dim) ? it->second : it->second.transpose();
    }
}
