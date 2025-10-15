#ifndef ATOPO_TRAITS_H
#define ATOPO_TRAITS_H

namespace atopo {
    template<typename CellType>
    struct CellTraits {
        static constexpr int dimension = -1;
    };

    template<typename TopologySource>
    struct TopologySourceTraits { /* User must specialize */ };
}
#endif // ATOPO_TRAITS_H
