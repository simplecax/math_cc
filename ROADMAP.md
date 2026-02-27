# atopo Evolution Roadmap: From Geometry to Persistent Topology

This document outlines the strategic plan to evolve the `atopo` library into a high-performance, CAD-aware topological analysis framework. We will synthesize the geometric orientation capabilities of the original `atopo` design with the highly efficient persistent homology algorithms inspired by `libstick` into a **single, unified C++20 codebase**.

## Core Architectural Principles
1.  **C++20 Concepts & Traits First:** The library will be strictly built upon C++20 `<concepts>`. It must seamlessly accept user-defined custom structs for `Cell`, `Simplex`, `CellComplex`, and `SimplicialComplex` through strong typing, avoiding virtual inheritance overhead.
2.  **Unified Codebase:** `libstick` will not exist as an external module or dependency. Its core algorithmic wisdom (XOR reduction, heap-sorted boolean vectors) will be rewritten, modernized, and deeply integrated into the native `atopo` data structures.
3.  **Eigen Backend:** Maintain Eigen3 as the primary linear algebra engine, optimized for sparse representations.

---

## Phase 1: Modernization & Trait Definition (C++20 Foundation)
**Goal:** Define the exact structural requirements for topologies using C++20 Concepts.

1.  **Define Topo Concepts:**
    *   `atopo::concepts::Simplex`: Requires a fixed dimension and ordered vertices.
    *   `atopo::concepts::Cell`: Requires a boundary iterator returning signed incidences.
    *   `atopo::concepts::Complex`: Requires ability to yield incidence matrices/iterators between dimensions.
2.  **Trait Adapters (`TopologySourceTraits`):**
    *   Upgrade the existing traits system so users can plug in their own CAD `BRepMesh` or Point Cloud `PointCloud` directly into the solver without intermediate conversions if they satisfy the concepts.
3.  **Build System Upgrade:**
    *   Update `CMakeLists.txt` to strictly enforce `CMAKE_CXX_STANDARD 20`.

## Phase 2: Unified Matrix & DSA Core
**Goal:** Merge the dual requirements of $\mathbb{Z}$ (oriented) and $\mathbb{Z}_2$ (unoriented) into a single, high-performance data structure layer.

1.  **Dual-View Incidence Structures:**
    *   The core storage remains `Eigen::SparseMatrix<int>` (or `short`) to preserve orientation.
    *   Implement zero-cost $\mathbb{Z}_2$ views (`atopo::views::Z2`) that project signed integers to booleans for fast reduction.
2.  **Modernized Reduction Engine:**
    *   Rewrite `libstick`'s `booleanvector.h` into modern C++20 `std::vector`-based XOR engines (`atopo::math::Z2Vector`).
    *   Implement hybrid solvers (Union-Find for 1-chains, Heap-sorted XOR for higher dimensions) directly within `atopo::math`.

## Phase 3: Oriented Chain Calculus & Operators
**Goal:** Provide the Discrete Exterior Calculus (DEC) toolkit.

1.  **Algebraic Operators:**
    *   Implement `atopo::boundary` ($\partial$) and `atopo::coboundary` ($\delta$) generic functions.
    *   Implement `atopo::laplacian` ($L = \partial \delta + \delta \partial$) returning standard Eigen matrices for physics/flow solvers.
2.  **Cellular Operations:**
    *   Ensure all operators respect the signed incidence values, allowing for calculations like flux, circulation, and winding numbers on arbitrary polyhedral meshes.

## Phase 4: Integrated Persistent Homology
**Goal:** Native support for TDA (Topological Data Analysis) over custom data types.

1.  **Filtration Concepts:**
    *   Define `atopo::concepts::Filtration`, allowing users to provide custom functors that assign birth times (real values) to cells.
2.  **Persistence Engine:**
    *   Implement `atopo::compute_persistence()`. It will take any `Complex` satisfying the concept and a `Filtration`, applying the Phase 2 reduction engine.
    *   Output standardized `atopo::PersistenceDiagram` and `atopo::Barcode` structures.

## Phase 5: Industrial API & Ecosystem
**Goal:** Make topological analysis accessible and actionable for CAD/CAE engineers.

1.  **Semantic Queries:**
    *   Provide high-level, concept-constrained algorithms: `atopo::query::has_void(complex)`, `atopo::query::connected_components(complex)`.
2.  **Visualization & Export:**
    *   Add exporters (VTK/OBJ) capable of highlighting "Critical Cells" (the specific user structs that triggered topological events).

---

## Architectural Comparison

| Feature | Legacy atopo | Legacy libstick | New Unified C++20 atopo |
| :--- | :--- | :--- | :--- |
| **Language** | C++14 | C++98 | **C++20 (Concepts/Ranges)** |
| **Data Binding** | Hardcoded/Inheritance| Hardcoded Structs | **Concept-based Custom Structs** |
| **Coefficients** | $\mathbb{Z}$ (short) | $\mathbb{Z}_2$ (bool) | **$\mathbb{Z}$ storage with $\mathbb{Z}_2$ zero-cost views** |
| **Solvers** | LinBox (Heavy) | Standard Reduction | **Native Hybrid (UF + Heap XOR)** |
| **Topology Type**| Cellular (Polygons) | Simplicial | **Both (via Concepts)** |
