#ifndef ATOPO_COMPLEX_H
#define ATOPO_COMPLEX_H

#include <atopo/core.h>

namespace atopo {
    namespace detail {
        struct SNFResult {
            size_t rank = 0;
            std::vector<long> torsion_coeffs;
        };
        SNFResult compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat);
    }

    class CellComplex {
    private:
        std::map<int, size_t> m_cell_counts;
        std::map<int, IncidenceMatrix<IncidenceCoeff>> m_boundary_maps;
    public:
        void setCellCount(int dim, size_t count);
        void setBoundaryOperator(int source_dim, IncidenceMatrix<IncidenceCoeff>&& map);
        [[nodiscard]] size_t getNumberOfCells(int dim) const;
        [[nodiscard]] IncidenceMatrix<IncidenceCoeff> getIncidenceMap(int from_dim, int to_dim) const;
        [[nodiscard]] std::map<int, HomologyGroup> computeHomology() const;
    };
}
#endif // ATOPO_COMPLEX_H
