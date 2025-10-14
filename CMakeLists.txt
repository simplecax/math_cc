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

# --- Library Target ---
# Define a header-only INTERFACE library for our framework.
# This makes it reusable and clearly defines its usage requirements.
add_library(atopo INTERFACE)

# The only public usage requirement for atopo is its include directory.
target_include_directories(atopo INTERFACE
    inc
)

# The atopo library itself needs Eigen3.
target_link_libraries(atopo INTERFACE
    Eigen3::Eigen
)

# --- Testing Setup with Google Test ---
enable_testing()

# Add the test executable
add_executable(atopo_tests
    test/test_atopo.cpp
    test/test_homology.cpp
)

# Link the test executable against our atopo library and GTest.
# It will inherit the include paths and Eigen dependency from atopo.
target_link_libraries(atopo_tests PRIVATE
    atopo
    GTest::Main
)

# Automatically discover all tests in the executable
include(GoogleTest)
gtest_discover_tests(atopo_tests)


# --- Build Information ---
message(STATUS "Project Name: ${PROJECT_NAME}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "Found GTest: ${GTEST_INCLUDE_DIRS}")
