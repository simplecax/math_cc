#ifndef ATOPO_CHAIN_H
#define ATOPO_CHAIN_H

#include <atopo/core.h>

namespace atopo {
    template<typename T>
    struct ChainBase {
        int dimension;
        Eigen::SparseVector<Coefficient<T>> data;
        ChainBase(int p, const Eigen::SparseVector<Coefficient<T>>& vec) : dimension(p), data(vec) {}
        ChainBase(int p, Eigen::SparseVector<Coefficient<T>>&& vec) : dimension(p), data(std::move(vec)) {}
    };

    template<typename T> struct Chain : public ChainBase<T> { using ChainBase<T>::ChainBase; };
    template<typename T> struct Cochain : public ChainBase<T> { using ChainBase<T>::ChainBase; };
}
#endif // ATOPO_CHAIN_H
