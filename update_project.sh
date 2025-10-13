#!/bin/bash

# This script refactors the single-file Algebraic Topology project into
# a structured format with a separate header (inc/atopo.h), a main source
# file (src/main.cpp), and an updated CMakeLists.txt.

# Exit immediately if any command fails.
set -e

echo "⚙️  Starting project refactor..."
echo "-----------------------------------"

# Step 1: Create the 'inc' directory if it doesn't exist.
echo "1. Creating directory 'inc/'..."
mkdir -p inc
echo "✅ Directory created."
echo ""

# Step 2: Create the new header file 'inc/atopo.h'.
# The 'cat << EOF' syntax creates a file from the text that follows.
echo "2. Creating new header file 'inc/atopo.h'..."
cat << 'EOF' > inc/atopo.h
#ifndef AT_OPO_H
#define AT_OPO_H

#include <iostream>
#include <vector>
#include <cstddef>
#include <memory>
#include <stdexcept>

// Eigen is a header-only library for linear algebra.
#include <Eigen/Sparse>

/**
 * @namespace atopo (Algebraic Topology)
 * @brief The core namespace for the generic topology framework.
 */
namespace atopo {

    // --- Core Data Types ---
    template<typename T>
    using Coefficient = T; /// The numeric type for chain coefficients (e.g., int, double).

    template<typename T>
    using IncidenceMatrix = Eigen::SparseMatrix<Coefficient<T>>; /// The boundary map matrix.

    // --- TRAITS (Interfaces for User Specialization) ---
    template<typename CellType>
    struct CellTraits { /* User specializes */ };

    template<typename TopologySource>
    struct TopologySourceTraits { /* User specializes */ };


    // --- CORE FRAMEWORK CLASSES & FUNCTIONS ---

    /**
     * @struct ChainBase
     * @brief A common base struct for Chains and Cochains to avoid code duplication.
     */
    template<typename T>
    struct ChainBase {
        int dimension;
        Eigen::SparseVector<Coefficient<T>> data;

        // Constructor
        ChainBase(int p, const Eigen::SparseVector<Coefficient<T>>& vec) : dimension(p), data(vec) {}
        ChainBase(int p, Eigen::SparseVector<Coefficient<T>>&& vec) : dimension(p), data(std::move(vec)) {}
    };

    /**
     * @struct Chain
     * @brief Represents a p-chain. Inherits its structure from ChainBase.
     */
    template<typename T>
    struct Chain : public ChainBase<T> {
        // Inherit constructors from the base class
        using ChainBase<T>::ChainBase;
    };
    
    /**
     * @struct Cochain
     * @brief Represents a p-cochain. Inherits its structure from ChainBase.
     * It is a distinct type from Chain, ensuring type safety.
     */
    template<typename T>
    struct Cochain : public ChainBase<T> {
        // Inherit constructors from the base class
        using ChainBase<T>::ChainBase;
    };

    /**
     * @class CellComplex
     * @brief A concrete, efficient data structure holding the topological space.
     */
    template<typename T>
    class CellComplex {
    private:
        std::vector<IncidenceMatrix<T>> marrBoundaryMaps;
        std::vector<size_t> marrCellCounts;

    public:
        CellComplex(size_t max_dim) : marrBoundaryMaps(max_dim), marrCellCounts(max_dim + 1, 0) {}
        
        void setBoundaryMap(size_t dim_of_target_chain, IncidenceMatrix<T>&& boundary_map) {
            if (dim_of_target_chain >= marrBoundaryMaps.size()) throw std::out_of_range("Dimension index is out of range.");
            marrBoundaryMaps[dim_of_target_chain] = std::move(boundary_map);
        }

        void setCellCount(size_t dim, size_t count) {
            if (dim >= marrCellCounts.size()) throw std::out_of_range("Dimension index is out of range.");
            marrCellCounts[dim] = count;
        }

        [[nodiscard]] const IncidenceMatrix<T>& getBoundaryMap(size_t dim_of_source_chain) const {
            if (dim_of_source_chain == 0 || dim_of_source_chain > marrBoundaryMaps.size()) throw std::out_of_range("Dimension is out of range for boundary map.");
            return marrBoundaryMaps[dim_of_source_chain - 1];
        }
        
        [[nodiscard]] size_t getNumberOfCells(size_t dim) const {
            if (dim >= marrCellCounts.size()) throw std::out_of_range("Dimension index out of range for cell count.");
            return marrCellCounts[dim];
        }

        [[nodiscard]] size_t getMaxDimension() const {
            return marrBoundaryMaps.size();
        }
    };

    /**
     * @brief Computes the boundary (∂) of a p-chain.
     * @return A (p-1)-chain representing the boundary.
     */
    template<typename T>
    [[nodiscard]] Chain<T> boundary(const CellComplex<T>& complex, const Chain<T>& chain) {
        if (chain.dimension <= 0) return Chain<T>(-1, Eigen::SparseVector<T>());

        const IncidenceMatrix<T>& d_p = complex.getBoundaryMap(chain.dimension);
        Eigen::SparseVector<T> boundary_vector = d_p * chain.data;
        boundary_vector.prune(0, 0);
        
        return Chain<T>(chain.dimension - 1, std::move(boundary_vector));
    }

    /**
     * @brief Computes the coboundary (δ) of a p-cochain.
     * @return A (p+1)-cochain.
     */
    template<typename T>
    [[nodiscard]] Cochain<T> coboundary(const CellComplex<T>& complex, const Cochain<T>& cochain) {
        int p = cochain.dimension;
        if (p < 0 || static_cast<size_t>(p) >= complex.getMaxDimension()) {
            return Cochain<T>(p + 1, Eigen::SparseVector<T>());
        }

        const IncidenceMatrix<T>& d_p_plus_1 = complex.getBoundaryMap(p + 1);
        Eigen::SparseVector<T> coboundary_vector = d_p_plus_1.transpose() * cochain.data;
        coboundary_vector.prune(0, 0);
        
        return Cochain<T>(p + 1, std::move(coboundary_vector));
    }

    /**
     * @brief Builder function to construct a `CellComplex` from a legacy source.
     */
    template<typename T, typename TopologySource>
    [[nodiscard]] CellComplex<T> create_complex_from_source(const TopologySource& source) {
        return TopologySourceTraits<TopologySource>::template build<T>(source);
    }

} // namespace atopo


// --- USER'S LEGACY CODE ---
struct LegacyVertex { size_t id; };
struct LegacyEdge { size_t id; size_t v_start, v_end; };
struct LegacyFace { size_t id; std::vector<size_t> edge_loop; };

struct LegacyMesh {
    std::vector<LegacyVertex> vertices;
    std::vector<LegacyEdge> edges;
    std::vector<LegacyFace> faces;
};


// --- GLUE CODE ---
namespace atopo {
    template<> struct CellTraits<LegacyVertex> { static constexpr int dimension = 0; };
    template<> struct CellTraits<LegacyEdge> { static constexpr int dimension = 1; };
    template<> struct CellTraits<LegacyFace> { static constexpr int dimension = 2; };

    template<> struct TopologySourceTraits<LegacyMesh> {
        template<typename T>
        static CellComplex<T> build(const LegacyMesh& mesh) {
            CellComplex<T> complex(2);
            complex.setCellCount(0, mesh.vertices.size());
            complex.setCellCount(1, mesh.edges.size());
            complex.setCellCount(2, mesh.faces.size());

            IncidenceMatrix<T> d1(mesh.vertices.size(), mesh.edges.size());
            for (const auto& edge : mesh.edges) {
                d1.insert(edge.v_start, edge.id) = -1; d1.insert(edge.v_end,   edge.id) = 1;
            }
            complex.setBoundaryMap(0, std::move(d1));

            IncidenceMatrix<T> d2(mesh.edges.size(), mesh.faces.size());
            for (const auto& face : mesh.faces) {
                for (size_t edge_id : face.edge_loop) { d2.insert(edge_id, face.id) = 1; }
            }
            complex.setBoundaryMap(1, std::move(d2));

            return complex;
        }
    };
} // namespace atopo


#endif // AT_OPO_H
EOF
echo "✅ Header file created."
echo ""

# Step 3: Overwrite 'src/main.cpp' with the new, shorter version.
echo "3. Overwriting 'src/main.cpp'..."
cat << 'EOF' > src/main.cpp
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
EOF
echo "✅ Source file 'src/main.cpp' overwritten."
echo ""

# Step 4: Overwrite 'CMakeLists.txt' with the correct version.
echo "4. Overwriting 'CMakeLists.txt'..."
cat << 'EOF' > CMakeLists.txt
# CMake configuration for the Algebraic Topology Framework
cmake_minimum_required(VERSION 3.10)
project(AlgebraicTopologyFramework LANGUAGES CXX)

# Set build type to Debug to include debugging symbols
set(CMAKE_BUILD_TYPE Debug)

# Set the C++ standard to 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Find Dependencies ---
# Your code uses Eigen for sparse matrices.
find_package(Eigen3 REQUIRED)

# --- Define Executable Target ---
# Create an executable named 'atopo_app' from the source file in src/
add_executable(atopo_app src/main.cpp)

# --- Configure Target ---
# Add our 'inc' directory to the include path so #include "atopo.h" works.
# Also add the Eigen include directory.
target_include_directories(atopo_app PRIVATE
    inc
    \${EIGEN3_INCLUDE_DIR}
)

# Link the executable against Eigen.
# For header-only libraries, this is the modern CMake way to propagate
# usage requirements (like include paths) to the target.
target_link_libraries(atopo_app PRIVATE Eigen3::Eigen)


# --- Build Information ---
message(STATUS "Project Name: \${PROJECT_NAME}")
message(STATUS "C++ Standard: \${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type: \${CMAKE_BUILD_TYPE}")
message(STATUS "Found Eigen3: \${EIGEN3_INCLUDE_DIR}")
EOF
echo "✅ CMakeLists.txt overwritten."
echo "-----------------------------------"

echo "🚀 Project refactor complete!"
echo "You can now build the new project structure with:"
echo "  cmake -B build"
echo "  cmake --build build"