#include <gtest/gtest.h>
#include "atopo.h"
#include <vector>

TEST(HomologyTorsionTest, RealProjectivePlane) {
    LegacyMesh mesh;
    // A square
    mesh.vertices = {{0,{}}, {1,{}}, {2,{}}, {3,{}}};
    // Edges of the square, with one diagonal
    mesh.edges = {
        {0, 0, 1, {}}, {1, 1, 2, {}}, {2, 2, 3, {}}, {3, 3, 0, {}},
        // New edge identifying opposite sides with a twist
        {4, 1, 3, {}} 
    };
    // One face with twisted identification
    // e0 + e4 - e2 - e3
    mesh.faces = { {0, {0, 4, 2, 3}, {}} }; 
    // This mesh is a bit of a simplification, a better one would be two faces.
    // Let's use a standard triangulation for RP2:
    mesh.vertices = {{0,{}}, {1,{}}, {2,{}}, {3,{}}, {4,{}}, {5,{}}};
    mesh.edges = {
        {0,0,1,{}}, {1,0,2,{}}, {2,0,4,{}}, {3,1,2,{}}, {4,1,3,{}},
        {5,1,5,{}}, {6,2,3,{}}, {7,2,4,{}}, {8,3,4,{}}, {9,3,5,{}},
        {10,4,5,{}}
    };
    mesh.faces = {
        {0,{0,1,3},{}}, {1,{0,5,4},{}}, {2,{1,2,6},{}}, {3,{1,4,3},{}},
        {4,{2,7,6},{}}, {5,{2,1,0},{}}, {6,{3,8,9},{}}, {7,{3,4,8},{}},
        {8,{4,7,2},{}}, {9,{4,5,10},{}}
    };

    atopo::CellComplex complex = atopo::create_complex_from_source(mesh);
    auto homology_groups = complex.computeHomology();
    
    // Expected for RP2: H0=Z, H1=Z_2, H2=0
    // Betti numbers: β0=1, β1=0, β2=0
    // Torsion: T1={2}
    
    ASSERT_TRUE(homology_groups.count(0));
    EXPECT_EQ(homology_groups.at(0).rank, 1);
    EXPECT_TRUE(homology_groups.at(0).torsion_coeffs.empty());

    ASSERT_TRUE(homology_groups.count(1));
    EXPECT_EQ(homology_groups.at(1).rank, 0);
    ASSERT_EQ(homology_groups.at(1).torsion_coeffs.size(), 1);
    EXPECT_EQ(homology_groups.at(1).torsion_coeffs[0], 2);
    
    ASSERT_TRUE(homology_groups.count(2));
    EXPECT_EQ(homology_groups.at(2).rank, 0);
    EXPECT_TRUE(homology_groups.at(2).torsion_coeffs.empty());
}
