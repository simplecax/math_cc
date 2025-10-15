#include <gtest/gtest.h>
#include "test_common.h"

// This test verifies that the mesh builder correctly throws an exception
// when given a face with an open (non-closed) edge loop.
TEST(ValidationTest, OpenLoopThrowsException) {
    LegacyMesh mesh;
    mesh.vertices = {{0,{}}, {1,{}}, {2,{}}, {3,{}}};
    // Edges {0, 1, 2} form an open path from v0 to v3.
    mesh.edges = { {0,0,1,{}}, {1,1,2,{}}, {2,2,3,{}} };
    // This face uses the open loop.
    mesh.faces = { {0, {0, 1, 2}, {}} }; 

    // We expect the builder to throw a std::runtime_error.
    EXPECT_THROW(atopo::create_complex_from_source(mesh), std::runtime_error);
}
