#ifndef ATOPO_DYNAMIC_COMPLEX_H
#define ATOPO_DYNAMIC_COMPLEX_H

#include <vector>
#include <map>
#include <atopo/core.h>
#include <atopo/concepts.h>
#include <atopo/complex.h>

namespace atopo {

    /**
     * @brief Represents a single boundary entry in a dynamic cell.
     */
    struct DynamicBoundaryEntry {
        std::size_t index;
        int coefficient;
    };

    /**
     * @brief A topology-only cell for the dynamic layer.
     * Caller manages business payloads (e.g. geometric coordinates) externally, 
     * mapping them 1:1 via the cell's index.
     */
    class DynamicCell {
    private:
        int m_dimension;
        std::vector<DynamicBoundaryEntry> m_boundary;
    public:

        explicit DynamicCell(int dim) : m_dimension(dim) {}
        
        [[nodiscard]] int dimension() const { return m_dimension; }
        [[nodiscard]] const std::vector<DynamicBoundaryEntry>& boundary() const { return m_boundary; }
        
        void add_boundary(std::size_t idx, int coeff) {
            m_boundary.push_back({idx, coeff});
        }
    };

    /**
     * @brief The Dynamic Layer for topology construction and high-frequency updates.
     * Satisfies atopo::concepts::Complex, enabling zero-copy traversal or compilation
     * into the Static Layer (CellComplex).
     */
    class DynamicComplex {
    private:
        int m_max_dim = -1;
        std::map<int, std::vector<DynamicCell>> m_cells;

    public:
        DynamicComplex() = default;

        [[nodiscard]] int max_dim() const { return m_max_dim; }
        
        [[nodiscard]] std::size_t cell_count(int dim) const {
            auto it = m_cells.find(dim);
            if (it != m_cells.end()) {
                return it->second.size();
            }
            return 0;
        }

        [[nodiscard]] const std::vector<DynamicCell>& cells(int dim) const {
            static const std::vector<DynamicCell> empty_cells;
            auto it = m_cells.find(dim);
            if (it != m_cells.end()) {
                return it->second;
            }
            return empty_cells;
        }

        // --- Dynamic Operations ---
        
        /**
         * @brief Adds a new cell of the given dimension.
         * @return The topological index of the new cell (useful for external payload mapping).
         */
        std::size_t add_cell(int dim) {
            if (dim > m_max_dim) m_max_dim = dim;
            std::size_t idx = m_cells[dim].size();
            m_cells[dim].emplace_back(dim);
            return idx;
        }

        /**
         * @brief Adds a boundary relation (e.g. a 0-cell bounding a 1-cell).
         */
        void add_boundary_relation(int dim, std::size_t cell_idx, std::size_t boundary_cell_idx, int coeff) {
            if (dim <= 0) return;
            if (cell_idx < m_cells[dim].size()) {
                m_cells[dim][cell_idx].add_boundary(boundary_cell_idx, coeff);
            }
        }

        /**
         * @brief Compiles the dynamic complex into a highly optimized StaticComplex.
         */
        CellComplex compile() const {
            return CellComplex::build(*this);
        }
    };

} // namespace atopo

#endif // ATOPO_DYNAMIC_COMPLEX_H
