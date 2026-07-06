# math_cc

`math_cc` is a lightweight, pure algebraic topology library for computing homology groups (Betti numbers, torsion).

## Architecture: Dual-Layer Topological Framework

To reconcile the conflicting demands of **high-frequency dynamic topological updates** (which require flexible data structures) and **highly optimized algebraic topology calculations** (which require continuous memory representations like Sparse Matrices), `math_cc` adopts a **Dual-Layer Architecture**:

1. **Dynamic Layer (`DynamicComplex`)**:
   - **Role**: Responsible for topological structure building, modification, and incremental updates (e.g., cell splitting, merging, appending).
   - **Data Structure**: Uses dictionary-based or adjacency-list based representations (e.g., `std::map`, `std::vector`), allowing $O(1)$ amortized insertion and modification of incidence relations.
   - **Usage**: Use this layer during your modeling or algorithm pipeline where topological mutations happen frequently.

2. **Static Layer (`StaticComplex` / `CellComplex`)**:
   - **Role**: Responsible purely for algebraic topology calculations (such as computing Homology Groups and Smith Normal Forms) and fast topological boundary queries.
   - **Data Structure**: Uses `Eigen::SparseMatrix` in Compressed Sparse Column (CSC) format for maximum memory locality and matrix operation efficiency.
   - **Usage**: Convert the `DynamicComplex` into a `StaticComplex` once the topology is stabilized (e.g., using `StaticComplex::build(dynamic_complex)`), and run algebraic queries on it.

## Trait-Based Payload Injection

The library uses C++20 Concepts to completely isolate topological representations from business logic. Users can define their custom cells (containing geometric coordinates, semantic payloads, etc.) and inject them into `math_cc`'s algebraic algorithms without modifying the underlying math structures. As long as your custom complex and cells satisfy `atopo::concepts::Complex` and `atopo::concepts::Cell`, they can be directly compiled into the Static Layer for computation.
