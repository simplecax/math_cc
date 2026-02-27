#ifndef ATOPO_CONCEPTS_H
#define ATOPO_CONCEPTS_H

#include <concepts>
#include <iterator>
#include <type_traits>

namespace atopo::concepts {

    /**
     * @brief Represents a single incidence entry (index and coefficient).
     */
    template<typename T>
    concept BoundaryEntry = requires(T e) {
        { e.index } -> std::convertible_to<std::size_t>;
        { e.coefficient } -> std::convertible_to<int>;
    };

    /**
     * @brief A Cell is a general topological building block.
     */
    template<typename T>
    concept Cell = requires(const T c) {
        { c.dimension() } -> std::convertible_to<int>;
        { c.boundary() } -> std::ranges::input_range;
    };

    /**
     * @brief A Complex is a collection of cells/simplices.
     */
    template<typename T>
    concept Complex = requires(const T c, int dim) {
        { c.cell_count(dim) } -> std::convertible_to<std::size_t>;
        // Should be able to iterate over all cells of a given dimension.
        { c.cells(dim) } -> std::ranges::input_range;
        requires Cell<std::ranges::range_value_t<decltype(std::declval<T>().cells(dim))>>;
    };

} // namespace atopo::concepts

#endif // ATOPO_CONCEPTS_H
