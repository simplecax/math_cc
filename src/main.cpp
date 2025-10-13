#include "atopo.h"
#include <iostream>
#include <string>

// Helper to print a chain
template<typename T>
void print_chain(const std::string& name, const atopo::Chain<T>& chain) {
    std::cout << name << " (a " << chain.dimension << "-chain):\n" << chain.data << "\n" << std::endl;
}

// Helper to print a cochain
template<typename T>
void print_cochain(const std::string& name, const atopo::Cochain<T>& cochain) {
    std::cout << name << " (a " << cochain.dimension << "-cochain):\n" << cochain.data << "\n" << std::endl;
}


int main() {
    // --- 1. SETUP & BUILD ---
    LegacyMesh my_mesh;
    my_mesh.vertices = {{0}, {1}, {2}, {3}};
    my_mesh.edges = {{0, 0, 1}, {1, 1, 2}, {2, 2, 3}, {3, 3, 0}};
    my_mesh.faces = {{0, {0, 1, 2, 3}}};
    
    using Coeff = int;

    auto complex = atopo::create_complex_from_source<Coeff>(my_mesh);
    std::cout << "--- `atopo::CellComplex` built successfully from `LegacyMesh`! ---\n\n";

    // --- 2. BOUNDARY EXAMPLE (∂² = 0) ---
    std::cout << "## Boundary Operator Example ##\n";
    Eigen::SparseVector<Coeff> face_vector(complex.getNumberOfCells(2));
    face_vector.insert(0) = 1;
    atopo::Chain<Coeff> face_chain(2, std::move(face_vector));
    print_chain("Face f0", face_chain);

    auto boundary_of_face = atopo::boundary(complex, face_chain);
    print_chain("Boundary of f0 (∂f0)", boundary_of_face);

    auto boundary_of_boundary = atopo::boundary(complex, boundary_of_face);
    print_chain("Boundary of boundary (∂²f0)", boundary_of_boundary);

    std::cout << "Is ∂²f0 empty? " << (boundary_of_boundary.data.nonZeros() == 0 ? "Yes" : "No") << "\n";
    std::cout << "--------------------------------\n\n";
    
    // --- 3. COBOUNDARY EXAMPLE (δ² = 0) ---
    std::cout << "## Coboundary Operator Example ##\n";
    Eigen::SparseVector<Coeff> v_cochain_vector(complex.getNumberOfCells(0));
    v_cochain_vector.insert(0) = 1;
    v_cochain_vector.insert(2) = -1;
    atopo::Cochain<Coeff> vertex_cochain(0, std::move(v_cochain_vector));
    print_cochain("Cochain f on Vertices", vertex_cochain);
    
    auto coboundary_of_v = atopo::coboundary(complex, vertex_cochain);
    print_cochain("Coboundary of f (δf)", coboundary_of_v);

    auto coboundary_of_coboundary = atopo::coboundary(complex, coboundary_of_v);
    print_cochain("Coboundary of coboundary (δ²f)", coboundary_of_coboundary);

    std::cout << "Is δ²f empty? " << (coboundary_of_coboundary.data.nonZeros() == 0 ? "Yes" : "No") << std::endl;

    return 0;
}
