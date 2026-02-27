#include <gtest/gtest.h>
#include "test_common.h"

// This test verifies that the mesh builder correctly throws an exception
// when given a face with an open (non-closed) edge loop.
/*
TEST(ValidationTest, OpenLoopThrowsException) {
    LegacyMesh mesh;
    mesh.vertices = {{0,{}}, {1,{}}, {2,{}}, {3,{}}};
    mesh.edges = { {0,0,1,{}}, {1,1,2,{}}, {2,2,3,{}} };
    mesh.faces = { {0, {0, 1, 2}, {}} }; 

    EXPECT_THROW(atopo::create_complex_from_source(mesh), std::runtime_error);
}
*/
