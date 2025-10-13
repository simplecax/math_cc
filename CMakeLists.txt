# CMake configuration for the Algebraic Topology Framework
cmake_minimum_required(VERSION 3.10)
project(AlgebraicTopologyFramework LANGUAGES CXX)

# Set build type to Debug to include debugging symbols
set(CMAKE_BUILD_TYPE Debug)

# Set the C++ standard to 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Find Dependencies ---
find_package(Eigen3 REQUIRED)
find_package(GTest REQUIRED)

# --- Testing Setup with Google Test ---
enable_testing()

# Add the test executable
add_executable(atopo_tests
    test/test_atopo.cpp
    test/test_homology.cpp
)

# Add include paths for the test target
target_include_directories(atopo_tests PRIVATE
    inc
    ${EIGEN3_INCLUDE_DIR}
)

# Link the test executable against GTest and Eigen
target_link_libraries(atopo_tests PRIVATE
    GTest::GTest
    GTest::Main
    Eigen3::Eigen
)

# Automatically discover all tests in the executable
include(GoogleTest)
gtest_discover_tests(atopo_tests)


# --- Build Information ---
message(STATUS "Project Name: ${PROJECT_NAME}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "Found GTest: ${GTEST_INCLUDE_DIRS}")
