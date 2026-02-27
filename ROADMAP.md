# atopo Evolution Roadmap: From Geometry to Persistent Topology (COMPLETED)

This document outlines the strategic plan to evolve the `atopo` library into a high-performance, CAD-aware topological analysis framework. We have successfully synthesized the geometric orientation capabilities of the original `atopo` design with the highly efficient persistent homology algorithms inspired by `libstick` into a **single, unified C++20 codebase**.

## Core Architectural Principles
1.  **C++20 Concepts & Traits First:** (DONE) The library is strictly built upon C++20 `<concepts>`. It accepts user-defined custom structs through strong typing.
2.  **Unified Codebase:** (DONE) Core algorithmic wisdom (XOR reduction, DSU) is modernized and deeply integrated.
3.  **Eigen Backend:** (DONE) Maintained Eigen3 as the primary linear algebra engine.

---

## Phase 1: Modernization & Trait Definition (C++20 Foundation) - COMPLETED
**Goal:** Define the exact structural requirements for topologies using C++20 Concepts.

1.  **[x] Define Topo Concepts:**
    *   `atopo::concepts::Simplex`, `Cell`, `Complex`, and `Filtration` defined in `concepts.h`.
2.  **[x] Trait Adapters:**
    *   Generic `CellComplex::build` accepts any type satisfying the `Complex` concept.
3.  **[x] Build System Upgrade:**
    *   `CMakeLists.txt` strictly enforces `CMAKE_CXX_STANDARD 20`.

## Phase 2: Unified Matrix & DSA Core - COMPLETED
**Goal:** Merge the dual requirements of $\mathbb{Z}$ (oriented) and $\mathbb{Z}_2$ (unoriented) into a single layer.

1.  **[x] Dual-View Incidence Structures:**
    *   Core storage in Eigen, with XOR-based reduction logic.
2.  **[x] Modernized Reduction Engine:**
    *   Implemented `atopo::math::Z2Vector` and `math::DSU` in `math.h`.
    *   Refactored `atopo_solver.cpp` to use the modernized math core.

## Phase 3: Oriented Chain Calculus & Operators - COMPLETED
**Goal:** Provide the Discrete Exterior Calculus (DEC) toolkit.

1.  **[x] Algebraic Operators:**
    *   Implemented generic `boundary`, `coboundary`, `laplacian`, and `laplacian_matrix` in `operators.h`.
2.  **[x] Cellular Operations:**
    *   Respects signed incidence values for flux/circulation analysis.

## Phase 4: Integrated Persistent Homology - COMPLETED
**Goal:** Native support for TDA (Topological Data Analysis) over custom data types.

1.  **[x] Filtration Concepts:**
    *   Defined `atopo::concepts::Filtration` for custom cell ordering.
2.  **[x] Persistence Engine:**
    *   Implemented `atopo::compute_persistence()` in `persistence.h`.
    *   Outputs `atopo::PersistenceDiagram` with birth/death tracking.

## Phase 5: Industrial API & Ecosystem - COMPLETED
**Goal:** Make topological analysis accessible and actionable for CAD/CAE engineers.

1.  **[x] Semantic Queries:**
    *   Added high-level API in `query.h`: `betti_number`, `is_connected`, `has_void`, `has_tunnel`.
2.  **[x] Visualization & Export:**
    *   Implemented `to_string()` and `summarize_persistence()` for easy analysis and export.

---

## Architectural Comparison (Final State)

| Feature | Legacy atopo | Legacy libstick | New Unified C++20 atopo |
| :--- | :--- | :--- | :--- |
| **Language** | C++14 | C++98 | **C++20 (Concepts/Ranges)** |
| **Data Binding** | Hardcoded/Inheritance| Hardcoded Structs | **Concept-based Custom Structs** |
| **Coefficients** | $\mathbb{Z}$ (short) | $\mathbb{Z}_2$ (bool) | **$\mathbb{Z}$ storage with XOR fast-path** |
| **Solvers** | LinBox (Heavy) | Standard Reduction | **Native Hybrid (UF + Optimized XOR)** |
| **Topology Type**| Cellular (Polygons) | Simplicial | **Both (via Concepts)** |
