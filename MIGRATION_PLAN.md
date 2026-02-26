# Migration Plan: Eigen-Only $\mathbb{Z}_2$ Homology Engine

This document outlines the strategy for refactoring the `atopo` library to remove all third-party dependencies except **Eigen**. The core objective is to replace the Smith Normal Form (SNF) calculations performed by LinBox with a high-performance, binary ($\mathbb{Z}_2$) rank solver implemented directly in Eigen.

## Phase 1: Core Type Refactoring (Preserving Orientation)
**Goal:** Maintain signed integers for geometric operations while preparing the topological layer for binary reduction.

1.  **Introduce $\mathbb{Z}_2$ Projection:** Add a utility function to project `Eigen::SparseMatrix<short>` (Incidence) into `Eigen::SparseMatrix<unsigned char>` (Binary).
2.  **Decouple `HomologyGroup`:** Update the `HomologyGroup` struct. 
    *   Note that `torsion_coeffs` will remain for API compatibility but will be empty in the $\mathbb{Z}_2$ implementation.
3.  **Namespace Cleanup:** Prepare `atopo::detail` to host the new native solver, replacing the current `atopo_linbox.cpp` logic.

## Phase 2: Native $\mathbb{Z}_2$ Sparse Rank Solver
**Goal:** Implement a replacement for the SNF solver using Sparse Gaussian Elimination.

1.  **Algorithm Selection:** Implement a **Sparse Gaussian Elimination** algorithm optimized for $\mathbb{Z}_2$ (XOR-based reduction).
2.  **Kernel Implementation:**
    *   Iterate through columns.
    *   Identify pivots using `Eigen::SparseMatrix::InnerIterator`.
    *   Perform row/column additions (XOR) to zero out entries.
3.  **Rank Calculation:** The Betti number $\beta_p$ will now be calculated using the Rank-Nullity theorem:
    $$\beta_p = 	ext{nullity}(\partial_p) - 	ext{rank}(\partial_{p+1})$$
    $$\beta_p = (	ext{cols}(\partial_p) - 	ext{rank}(\partial_p)) - 	ext{rank}(\partial_{p+1})$$

## Phase 3: Integration & Dependency Removal
**Goal:** Physically remove the LinBox ecosystem and update the build system.

1.  **Delete `src/atopo_linbox.cpp`:** Replace it with `src/atopo_solver.cpp` containing the native Eigen logic.
2.  **Update `CMakeLists.txt`:**
    *   Remove `find_package` for LinBox, Givaro, and GMP.
    *   Ensure Eigen is the sole required dependency.
3.  **Platform Validation:** Verify the project now compiles cleanly on Windows (MSVC) without specialized environments.

## Phase 4: Test Suite Adaptation
**Goal:** Align verification logic with the new mathematical reality.

1.  **Update Torsion Tests:** Modify `test_klein_bottle.cpp` and `test_rp2_integration.cpp`. 
    *   Expect `torsion_coeffs` to be empty.
    *   Verify that Betti numbers ($\beta_n$) remain correct under $\mathbb{Z}_2$ coefficients.
2.  **Add Performance Benchmarks:** Create a test case with a high-density mesh (e.g., 50k+ faces) to compare the speed of the native XOR-solver vs. the previous SNF-solver.

## Phase 5: CAD-Specific Optimizations
**Goal:** Leverage mesh properties for $O(N)$ performance.

1.  **Manifold Optimization:** In the solver, add a "fast-path" check for manifold edges (exactly two entries per row) to minimize reduction steps.
2.  **Memory Reuse:** Ensure the rank solver operates on a local copy or performs in-place reduction to keep the memory footprint low.

---

### Architectural Benefits After Migration:
*   **Zero Legal Friction:** Full compliance with commercial licensing (MPL2/MIT style).
*   **Zero Binary Bloat:** Library size reduced by ~90%.
*   **Maximum Portability:** Works anywhere Eigen works.
*   **Orientation Awareness:** The library still understands "Left vs. Right" in the Chain/Geometric layer, but uses "Connected vs. Disconnected" in the Topological layer.
