#ifndef ATOPO_COMPLEX_H
#define ATOPO_COMPLEX_H

#include <atopo/core.h>
#include <atopo/concepts.h>
#include <map>

namespace atopo {
    namespace detail { 
        HomologyGroup compute_snf_results(const IncidenceMatrix<IncidenceCoeff>& sparse_mat); 
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

        /**
         * @brief Generic builder that works with any type satisfying the Complex concept.
         */
        template<concepts::Complex T>
        static CellComplex build(const T& source) {
            CellComplex complex;
            int max_d = source.max_dim();
            
            for (int dim = 0; dim <= max_d; ++dim) {
                size_t count = source.cell_count(dim);
                complex.setCellCount(dim, count);

                if (dim > 0) {
                    IncidenceMatrix<IncidenceCoeff> d(source.cell_count(dim-1), count);
                    std::vector<Eigen::Triplet<IncidenceCoeff>> triplets;
                    // Heuristic reserve: assume an average of 4 boundary elements per cell
                    triplets.reserve(count * 4);
                    
                    size_t cell_idx = 0;
                    for (const auto& cell : source.cells(dim)) {
                        for (const auto& entry : cell.boundary()) {
                            triplets.emplace_back(entry.index, cell_idx, static_cast<IncidenceCoeff>(entry.coefficient));
                        }
                        cell_idx++;
                    }
                    d.setFromTriplets(triplets.begin(), triplets.end());
                    if (d.nonZeros() > 0) {
                        complex.setBoundaryOperator(dim, std::move(d));
                    }
                }
            }
            return complex;
        }
    };
    
    // Alias for the dual-layer architecture
    using StaticComplex = CellComplex;
}
#endif // ATOPO_COMPLEX_H
