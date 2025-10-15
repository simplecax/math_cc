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
find_package(OpenMP REQUIRED) # For fflas-ffpack, a LinBox dependency

# LinBox and its dependencies (Givaro, fflas-ffpack) use the older pkg-config
# system for discovery, not modern CMake config files.
find_package(PkgConfig REQUIRED)
pkg_check_modules(GIVARO REQUIRED givaro)
pkg_check_modules(FFLAS_FFPACK REQUIRED fflas-ffpack)
pkg_check_modules(LINBOX REQUIRED linbox)

# --- Library Target ---
# 'atopo' is a static library containing the core logic.
add_library(atopo STATIC
    src/atopo.cpp
    src/atopo_linbox.cpp
)

# Public includes are needed by any code that uses the atopo library.
target_include_directories(atopo PUBLIC
    inc
)

# Private includes are only needed to COMPILE the atopo library itself.
target_include_directories(atopo PRIVATE
    test # For test_common.h needed by atopo.cpp
    ${LINBOX_INCLUDE_DIRS}
    ${GIVARO_INCLUDE_DIRS}
    ${FFLAS_FFPACK_INCLUDE_DIRS}
)

# Link atopo against its dependencies.
target_link_libraries(atopo
    # Eigen is PUBLIC because atopo.h includes Eigen headers. Any code
    # including atopo.h must also have access to Eigen's headers.
    PUBLIC
        Eigen3::Eigen
    # The other dependencies are PRIVATE as they are implementation details
    # hidden within the atopo.cpp files.
    PRIVATE
        OpenMP::OpenMP_CXX
        ${LINBOX_LIBRARIES}
        ${GIVARO_LIBRARIES}
        ${FFLAS_FFPACK_LIBRARIES}
)

# --- Testing Setup ---
enable_testing()
add_executable(atopo_tests
    test/test_atopo.cpp
    test/test_homology.cpp
    test/test_rp2_integration.cpp
)

# The test executable needs to know where the custom-installed libraries are.
target_link_directories(atopo_tests PRIVATE ${LINBOX_LIBRARY_DIRS})
# The test target needs to find test_common.h
target_include_directories(atopo_tests PRIVATE test)

# Link the tests against our atopo library and GTest.
# All other dependencies (LinBox, Eigen, etc.) are inherited from atopo.
target_link_libraries(atopo_tests PRIVATE
    atopo
    GTest::Main
)

include(GoogleTest)
gtest_discover_tests(atopo_tests)

message(STATUS "Project Name: ${PROJECT_NAME}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Found LinBox libs: ${LINBOX_LIBRARIES}")
