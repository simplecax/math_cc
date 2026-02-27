# CMake configuration for the Algebraic Topology Framework
cmake_minimum_required(VERSION 3.10)
project(AlgebraicTopologyFramework LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Find Dependencies ---
find_package(Eigen3 REQUIRED)
find_package(GTest REQUIRED)

# --- Library Target ---
# 'atopo' is a static library containing the core logic.
add_library(atopo STATIC
    src/atopo.cpp
    src/atopo_solver.cpp
)

# Public includes are needed by any code that uses the atopo library.
target_include_directories(atopo PUBLIC
    inc
)

# Private includes are only needed to COMPILE the atopo library itself.
target_include_directories(atopo PRIVATE
    test # For test_common.h needed by atopo.cpp
)

# Link atopo against its dependencies.
target_link_libraries(atopo
    # Eigen is PUBLIC because atopo.h includes Eigen headers.
    PUBLIC
        Eigen3::Eigen
)

# --- Testing Setup ---
enable_testing()
add_executable(atopo_tests
    test/test_atopo.cpp
    test/test_homology.cpp
    test/test_rp2_integration.cpp
    test/test_validation.cpp
    test/test_klein_bottle.cpp
)

# The test target needs to find test_common.h
target_include_directories(atopo_tests PRIVATE test)

# Link the tests against our atopo library and GTest.
target_link_libraries(atopo_tests PRIVATE
    atopo
    GTest::Main
)

include(GoogleTest)
gtest_discover_tests(atopo_tests)

message(STATUS "Project Name: ${PROJECT_NAME}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Found LinBox libs: ${LINBOX_LIBRARIES}")
