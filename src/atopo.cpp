#include "atopo.h"

namespace atopo {
    // Member function definitions
    void CellComplex::setBoundaryOperator(int source_dim, IncidenceMatrix<IncidenceCoeff>&& map) { m_boundary_maps[source_dim] = std::move(map); }
    void CellComplex::setGeneralIncidenceMap(int from_dim, int to_dim, IncidenceMatrix<IncidenceCoeff>&& map) {
        int low_dim = std::min(from_dim, to_dim);
        int high_dim = std::max(from_dim, to_dim);
        DimensionPair key = {low_dim, high_dim};
        if (from_dim < to_dim) { m_general_incidence_maps[key] = std::move(map); } else { m_general_incidence_maps[key] = map.transpose(); }
    }
    size_t CellComplex::getNumberOfCells(int dim) const {
        auto it = m_cell_counts.find(dim);
        return (it == m_cell_counts.end()) ? 0 : it->second;
    }
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
        if (it == m_general_incidence_maps.end()) { return IncidenceMatrix<IncidenceCoeff>(getNumberOfCells(to_dim), getNumberOfCells(from_dim)); }
        return (from_dim < to_dim) ? it->second : it->second.transpose();
    }
    
    std::map<int, HomologyGroup> CellComplex::computeHomology() const {
        std::map<int, HomologyGroup> homology;
        int max_dim = 0;
        if(!m_cell_counts.empty()) { max_dim = m_cell_counts.rbegin()->first; }

        for (int p = 0; p <= max_dim + 1; ++p) {
             if (getNumberOfCells(p) == 0 && getNumberOfCells(p-1) == 0 && p > max_dim) continue;

            auto snf_dp = (p > 0) ? detail::compute_snf_results(getIncidenceMap(p, p - 1)) : detail::SNFResult{};
            auto snf_dp1 = detail::compute_snf_results(getIncidenceMap(p + 1, p));
            
            long rank_Cp = getNumberOfCells(p);
            long rank_dp = snf_dp.rank;
            long rank_dp1 = snf_dp1.rank;

            HomologyGroup hg;
            hg.rank = rank_Cp - rank_dp - rank_dp1;
            hg.torsion_coeffs = snf_dp1.torsion_coeffs; // Torsion of Hp is from ∂(p+1)
            
            if (hg.rank != 0 || !hg.torsion_coeffs.empty() || rank_Cp > 0) {
                homology[p] = hg;
            }
        }
        return homology;
    }
    
    // Templated free functions and glue code
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
        const IncidenceMatrix<IncidenceCoeff> delta_p = complex.getIncidenceMap(p + 1, p).transpose();
        Eigen::SparseVector<T> coboundary_vector = delta_p.template cast<T>() * cochain.data;
        coboundary_vector.prune(0, 0);
        return Cochain<T>(p + 1, std::move(coboundary_vector));
    }
    template<typename TopologySource>
    [[nodiscard]] CellComplex create_complex_from_source(const TopologySource& source) {
        return TopologySourceTraits<TopologySource>::build(source);
    }
    
    CellComplex TopologySourceTraits<LegacyMesh>::build(const LegacyMesh& mesh) {
        CellComplex complex;
        complex.setCellCount(0, mesh.vertices.size());
        complex.setCellCount(1, mesh.edges.size());
        complex.setCellCount(2, mesh.faces.size());
        
        IncidenceMatrix<IncidenceCoeff> d1(mesh.vertices.size(), mesh.edges.size());
        for (const auto& edge : mesh.edges) {
            d1.insert(edge.v_start, edge.id) = -1;
            d1.insert(edge.v_end, edge.id) = 1;
        }
        complex.setBoundaryOperator(1, std::move(d1));

        IncidenceMatrix<IncidenceCoeff> d2(mesh.edges.size(), mesh.faces.size());
        for (const auto& face : mesh.faces) {
            if (face.edge_loop.size() < 2) continue;

            // --- FINAL, CORRECT ORIENTATION LOGIC ---
            const auto& e0 = mesh.edges.at(face.edge_loop[0]);
            const auto& e1 = mesh.edges.at(face.edge_loop[1]);

            // Find the common vertex between the first two edges to establish traversal order.
            size_t v_start, v_mid, v_next;
            if (e0.v_end == e1.v_start) { v_start = e0.v_start; v_mid = e0.v_end; v_next = e1.v_end; }
            else if (e0.v_end == e1.v_end) { v_start = e0.v_start; v_mid = e0.v_end; v_next = e1.v_start; }
            else if (e0.v_start == e1.v_start) { v_start = e0.v_end; v_mid = e0.v_start; v_next = e1.v_end; }
            else { v_start = e0.v_end; v_mid = e0.v_start; v_next = e1.v_start; }
            
            size_t current_v = v_start;
            for (size_t edge_id : face.edge_loop) {
                const auto& edge = mesh.edges.at(edge_id);
                if (edge.v_start == current_v) {
                    d2.insert(edge.id, face.id) = 1; // Agrees with traversal
                    current_v = edge.v_end;
                } else {
                    d2.insert(edge.id, face.id) = -1; // Opposes traversal
                    current_v = edge.v_start;
                }
            }
        }
        complex.setBoundaryOperator(2, std::move(d2));
        return complex;
    }

    // Explicit template instantiations to ensure functions are compiled
    template Chain<int> boundary<int>(const CellComplex&, const Chain<int>&);
    template Cochain<int> coboundary<int>(const CellComplex&, const Cochain<int>&);
    template CellComplex create_complex_from_source<LegacyMesh>(const LegacyMesh&);
}
