# CMake configuration for the Algebraic Topology Framework
cmake_minimum_required(VERSION 3.10)
project(AlgebraicTopologyFramework LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Find Dependencies ---
find_package(Eigen3 REQUIRED)
find_package(GTest REQUIRED)
find_package(OpenMP REQUIRED) # Find the OpenMP library
find_package(PkgConfig REQUIRED)
pkg_check_modules(GIVARO REQUIRED givaro)
pkg_check_modules(FFLAS_FFPACK REQUIRED fflas-ffpack)
pkg_check_modules(LINBOX REQUIRED linbox)

# --- Library Target ---
add_library(atopo STATIC
    src/atopo.cpp
    src/atopo_linbox.cpp
)

# Public includes for any target linking against atopo
target_include_directories(atopo PUBLIC
    inc
)

# Private includes needed to COMPILE atopo itself
target_include_directories(atopo PRIVATE
    ${LINBOX_INCLUDE_DIRS}
    ${GIVARO_INCLUDE_DIRS}
    ${FFLAS_FFPACK_INCLUDE_DIRS}
)

# Link atopo against its dependencies
target_link_libraries(atopo PRIVATE
    PUBLIC Eigen3::Eigen # Public because atopo.h includes Eigen headers
    PRIVATE
        OpenMP::OpenMP_CXX # Add OpenMP as a private dependency
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

# Add the linker path for the custom /opt install
target_link_directories(atopo_tests PRIVATE ${LINBOX_LIBRARY_DIRS})

# Link the test executable against our atopo library and GTest.
target_link_libraries(atopo_tests PRIVATE
    atopo
    GTest::Main
)

include(GoogleTest)
gtest_discover_tests(atopo_tests)

# --- Build Information ---
message(STATUS "Project Name: ${PROJECT_NAME}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Found LinBox includes: ${LINBOX_INCLUDE_DIRS}")
message(STATUS "Found LinBox libs: ${LINBOX_LIBRARIES}")
