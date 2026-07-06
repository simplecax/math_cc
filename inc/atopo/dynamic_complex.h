#ifndef ATOPO_DYNAMIC_COMPLEX_H
#define ATOPO_DYNAMIC_COMPLEX_H

#include <vector>
#include <stdexcept>
#include <atopo/core.h>
#include <atopo/concepts.h>
#include <atopo/complex.h>

namespace atopo {

    /**
     * @brief Represents a single boundary entry in a dynamic cell.
     */
    struct DynamicBoundaryEntry {
        std::size_t index;
        IncidenceCoeff coefficient; // Optimized: Use smaller width math type instead of int
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
        
        void add_boundary(std::size_t idx, IncidenceCoeff coeff) {
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
        // Optimized: O(1) contiguous vector mapping instead of O(log N) tree mapping
        std::vector<std::vector<DynamicCell>> m_cells;

    public:
        DynamicComplex() = default;

        [[nodiscard]] int max_dim() const { return m_max_dim; }
        
        [[nodiscard]] std::size_t cell_count(int dim) const {
            if (dim >= 0 && dim < static_cast<int>(m_cells.size())) {
                return m_cells[dim].size();
            }
            return 0;
        }

        [[nodiscard]] const std::vector<DynamicCell>& cells(int dim) const {
            static const std::vector<DynamicCell> empty_cells;
            if (dim >= 0 && dim < static_cast<int>(m_cells.size())) {
                return m_cells[dim];
            }
            return empty_cells;
        }

        // --- Dynamic Operations ---
        
        /**
         * @brief Adds a new cell of the given dimension.
         * @return The topological index of the new cell (useful for external payload mapping).
         */
        std::size_t add_cell(int dim) {
            if (dim < 0) throw std::invalid_argument("Dimension cannot be negative.");
            if (dim > m_max_dim) m_max_dim = dim;
            if (dim >= static_cast<int>(m_cells.size())) {
                m_cells.resize(dim + 1);
            }
            std::size_t idx = m_cells[dim].size();
            m_cells[dim].emplace_back(dim);
            return idx;
        }

        /**
         * @brief Adds a boundary relation (e.g. a 0-cell bounding a 1-cell).
         */
        void add_boundary_relation(int dim, std::size_t cell_idx, std::size_t boundary_cell_idx, IncidenceCoeff coeff) {
            if (dim <= 0) return;
            // Defensive Programming: Fail-fast out of bounds check
            if (dim >= static_cast<int>(m_cells.size())) {
                throw std::out_of_range("Dimension out of range in add_boundary_relation");
            }
            if (cell_idx >= m_cells[dim].size()) {
                throw std::out_of_range("cell_idx out of bounds in add_boundary_relation");
            }
            m_cells[dim][cell_idx].add_boundary(boundary_cell_idx, coeff);
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
