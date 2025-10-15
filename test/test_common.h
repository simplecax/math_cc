#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include "atopo.h"

// --- Example User Code ---
// This code is an example of how a user would adapt their own data
// structures. It belongs with the tests, not in the library headers.
struct MyPoint3D { double x, y, z; };
struct MyCurve   { /* ... */ };
struct MySurface { /* ... */ };

struct LegacyVertex { size_t id; MyPoint3D geometry; };
struct LegacyEdge   { size_t id; size_t v_start, v_end; MyCurve geometry; };
struct LegacyFace   { size_t id; std::vector<size_t> edge_loop; MySurface geometry; };

struct LegacyMesh {
    std::vector<LegacyVertex> vertices;
    std::vector<LegacyEdge> edges;
    std::vector<LegacyFace> faces;
};

// --- Glue Code for the Example ---
namespace atopo {
    template<> struct CellTraits<LegacyVertex> { static constexpr int dimension = 0; };
    template<> struct CellTraits<LegacyEdge>   { static constexpr int dimension = 1; };
    template<> struct CellTraits<LegacyFace>   { static constexpr int dimension = 2; };

    template<> struct TopologySourceTraits<LegacyMesh> {
        static CellComplex build(const LegacyMesh& mesh);
    };
}

#endif // TEST_COMMON_H
