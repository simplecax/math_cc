#include <atopo.h>
#include <vector>
#include <memory>
#include <ranges>

// --- Example User Code ---
struct MyPoint3D { double x, y, z; };

struct BoundaryEntry {
    std::size_t index;
    int coefficient;
};

// Interface for type-erased cell access
struct ICell {
    virtual ~ICell() = default;
    virtual int dimension() const = 0;
    virtual std::vector<BoundaryEntry> boundary() const = 0;
};

// Wrapper to make it look like a value type for the concept
struct AnyCell {
    std::shared_ptr<const ICell> ptr;
    int dimension() const { return ptr->dimension(); }
    auto boundary() const { return ptr->boundary(); }
};

struct LegacyVertex : public ICell { 
    size_t id; 
    MyPoint3D geometry; 
    int dimension() const override { return 0; }
    std::vector<BoundaryEntry> boundary() const override { return {}; }
};

struct LegacyEdge : public ICell { 
    size_t id; 
    size_t v_start, v_end; 
    int dimension() const override { return 1; }
    std::vector<BoundaryEntry> boundary() const override {
        return {{v_start, -1}, {v_end, 1}};
    }
};

struct LegacyFace : public ICell { 
    size_t id; 
    std::vector<size_t> edge_loop; 
    std::vector<int> edge_orientations; 
    int dimension() const override { return 2; }
    std::vector<BoundaryEntry> boundary() const override {
        std::vector<BoundaryEntry> b;
        for (size_t i = 0; i < edge_loop.size(); ++i) {
            int coeff = (i < edge_orientations.size()) ? edge_orientations[i] : 1;
            b.push_back({edge_loop[i], coeff});
        }
        return b;
    }
};

struct LegacyMesh {
    std::vector<std::shared_ptr<LegacyVertex>> vertices;
    std::vector<std::shared_ptr<LegacyEdge>> edges;
    std::vector<std::shared_ptr<LegacyFace>> faces;

    void add_vertex(size_t id) {
        auto v = std::make_shared<LegacyVertex>();
        v->id = id;
        vertices.push_back(v);
    }

    void add_edge(size_t id, size_t v_start, size_t v_end) {
        auto e = std::make_shared<LegacyEdge>();
        e->id = id;
        e->v_start = v_start;
        e->v_end = v_end;
        edges.push_back(e);
    }

    void add_face(size_t id, std::vector<size_t> edge_loop, std::vector<int> orientations = {}) {
        auto f = std::make_shared<LegacyFace>();
        f->id = id;
        f->edge_loop = edge_loop;
        f->edge_orientations = orientations;
        faces.push_back(f);
    }

    std::size_t cell_count(int dim) const {
        if (dim == 0) return vertices.size();
        if (dim == 1) return edges.size();
        if (dim == 2) return faces.size();
        return 0;
    }

    std::vector<AnyCell> cells(int dim) const {
        std::vector<AnyCell> result;
        if (dim == 0) for(auto& v : vertices) result.push_back({v});
        else if (dim == 1) for(auto& e : edges) result.push_back({e});
        else if (dim == 2) for(auto& f : faces) result.push_back({f});
        return result;
    }
};

// --- Glue Code for Compatibility ---
namespace atopo {
    template<> struct TopologySourceTraits<LegacyMesh> {
        static CellComplex build(const LegacyMesh& mesh) {
            return CellComplex::build(mesh);
        }
    };
}
