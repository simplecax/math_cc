#ifndef ATOPO_OPERATORS_H
#define ATOPO_OPERATORS_H

#include <atopo/complex.h>
#include <atopo/chain.h>
#include <atopo/traits.h>

namespace atopo {
    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex& complex, const Chain<T>& chain);

    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex& complex, const Cochain<T>& cochain);

    template<typename TopologySource>
    [[nodiscard]] CellComplex create_complex_from_source(const TopologySource& source);
}
#endif // ATOPO_OPERATORS_H
