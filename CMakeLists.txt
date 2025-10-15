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

# Use PkgConfig to find LinBox and its dependencies
find_package(PkgConfig REQUIRED)
pkg_check_modules(GIVARO REQUIRED givaro)
pkg_check_modules(FFLAS_FFPACK REQUIRED fflas-ffpack)
pkg_check_modules(LINBOX REQUIRED linbox)

# --- Library Target ---
# Define a header-only INTERFACE library for our framework.
add_library(atopo INTERFACE)

# Define the usage requirements for the 'atopo' library.
# Any target that links to 'atopo' will inherit these.
target_include_directories(atopo INTERFACE
    inc                         # Its own headers
    ${LINBOX_INCLUDE_DIRS}      # LinBox headers
    ${GIVARO_INCLUDE_DIRS}      # Givaro headers
    ${FFLAS_FFPACK_INCLUDE_DIRS}# Fflas-ffpack headers
)
target_link_libraries(atopo INTERFACE
    Eigen3::Eigen
    ${LINBOX_LIBRARIES}
    ${GIVARO_LIBRARIES}
    ${FFLAS_FFPACK_LIBRARIES}
)

# --- Testing Setup with Google Test ---
enable_testing()

add_executable(atopo_tests
    test/test_atopo.cpp
    test/test_homology.cpp
)

# Link the test executable against our atopo library and GTest.
# It will now correctly inherit all necessary include paths and libraries from 'atopo'.
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
message(STATUS "Found LinBox includes: ${LINBOX_INCLUDE_DIRS}")
message(STATUS "Found LinBox libs: ${LINBOX_LIBRARIES}")
