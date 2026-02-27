#ifndef ATOPO_OPERATORS_H
#define ATOPO_OPERATORS_H

#include <atopo/complex.h>
#include <atopo/chain.h>

namespace atopo {
    // --- TEMPLATE DEFINITIONS MOVED HERE ---
    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex& complex, const Chain<T>& chain) {
        if (chain.dimension <= 0) return Chain<T>(-1, Eigen::SparseVector<T>());
        const IncidenceMatrix<IncidenceCoeff>& d_p = complex.getIncidenceMap(chain.dimension, chain.dimension - 1);
        Eigen::SparseVector<T> boundary_vector = d_p.template cast<T>() * chain.data;
        boundary_vector.prune(0, 0);
        return Chain<T>(chain.dimension - 1, std::move(boundary_vector));
    }

    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex& complex, const Cochain<T>& cochain) {
        int p = cochain.dimension;
        // In DEC, the coboundary operator delta^p is the transpose of the boundary operator d_{p+1}
        const IncidenceMatrix<IncidenceCoeff> delta_p = complex.getIncidenceMap(p + 1, p).transpose();
        Eigen::SparseVector<T> coboundary_vector = delta_p.template cast<T>() * cochain.data;
        coboundary_vector.prune(0, 0);
        return Cochain<T>(p + 1, std::move(coboundary_vector));
    }

    /**
     * @brief Computes the combinatorial Laplacian matrix for p-forms.
     * L_p = d_{p+1} * delta_p + delta_{p-1} * d_p
     */
    template<typename T = double>
    [[nodiscard]] Eigen::SparseMatrix<T> laplacian_matrix(const CellComplex& complex, int p) {
        // d_{p+1} maps (p+1) -> p
        const Eigen::SparseMatrix<T> d_p1 = complex.getIncidenceMap(p + 1, p).template cast<T>();
        // d_p maps p -> (p-1)
        const Eigen::SparseMatrix<T> d_p = complex.getIncidenceMap(p, p - 1).template cast<T>();
        
        // L_p = d_{p+1} * (d_{p+1})^T + (d_p)^T * d_p
        Eigen::SparseMatrix<T> L(complex.getNumberOfCells(p), complex.getNumberOfCells(p));
        L.setZero();

        L = d_p1 * d_p1.transpose(); 
        if (p > 0) {
            L += Eigen::SparseMatrix<T>(d_p.transpose() * d_p);
        }
        return L;
    }

    /**
     * @brief Applies the p-form Laplacian to a Cochain.
     */
    template<typename T>
    [[nodiscard]] Cochain<T> laplacian(const CellComplex& complex, const Cochain<T>& cochain) {
        int p = cochain.dimension;
        auto L = laplacian_matrix<T>(complex, p);
        Eigen::SparseVector<T> result_vec = L * cochain.data;
        result_vec.prune(0, 0);
        return Cochain<T>(p, std::move(result_vec));
    }

    template<concepts::Complex T>
    [[nodiscard]] CellComplex create_complex_from_source(const T& source) {
        return CellComplex::build(source);
    }
}
#endif // ATOPO_OPERATORS_H
