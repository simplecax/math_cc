#include <atopo.h>
#include <vector>
#include <memory>
#include <ranges>
#include <variant>

// --- Example User Code ---
struct MyPoint3D { double x, y, z; };

struct BoundaryEntry {
    std::size_t index;
    int coefficient;
};

struct LegacyVertex { 
    size_t id; 
    MyPoint3D geometry; 
    int dimension() const { return 0; }
    std::vector<BoundaryEntry> boundary() const { return {}; }
};

struct LegacyEdge { 
    size_t id; 
    size_t v_start, v_end; 
    int dimension() const { return 1; }
    std::vector<BoundaryEntry> boundary() const {
        return {{v_start, -1}, {v_end, 1}};
    }
};

struct LegacyFace { 
    size_t id; 
    std::vector<size_t> edge_loop; 
    std::vector<int> edge_orientations; 
    int dimension() const { return 2; }
    std::vector<BoundaryEntry> boundary() const {
        std::vector<BoundaryEntry> b;
        for (size_t i = 0; i < edge_loop.size(); ++i) {
            int coeff = (i < edge_orientations.size()) ? edge_orientations[i] : 1;
            b.push_back({edge_loop[i], coeff});
        }
        return b;
    }
};

using AnyCellVariant = std::variant<LegacyVertex, LegacyEdge, LegacyFace>;

// Wrapper to make variant satisfy the Cell concept
struct AnyCell {
    AnyCellVariant v;
    int dimension() const {
        return std::visit([](const auto& c) { return c.dimension(); }, v);
    }
    std::vector<BoundaryEntry> boundary() const {
        return std::visit([](const auto& c) { return c.boundary(); }, v);
    }
};

struct LegacyMesh {
    std::vector<LegacyVertex> vertices;
    std::vector<LegacyEdge> edges;
    std::vector<LegacyFace> faces;

    int max_dim() const { return 2; }

    void add_vertex(size_t id) {
        vertices.push_back({id, {}});
    }

    void add_edge(size_t id, size_t v_start, size_t v_end) {
        edges.push_back({id, v_start, v_end});
    }

    void add_face(size_t id, std::vector<size_t> edge_loop, std::vector<int> orientations = {}) {
        faces.push_back({id, edge_loop, orientations});
    }

    std::size_t cell_count(int dim) const {
        if (dim == 0) return vertices.size();
        if (dim == 1) return edges.size();
        if (dim == 2) return faces.size();
        return 0;
    }

    std::vector<AnyCell> cells(int dim) const {
        std::vector<AnyCell> result;
        if (dim == 0) for(const auto& v : vertices) result.push_back({v});
        else if (dim == 1) for(const auto& e : edges) result.push_back({e});
        else if (dim == 2) for(const auto& f : faces) result.push_back({f});
        return result;
    }
};
