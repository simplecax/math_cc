# CMake configuration for the Algebraic Topology Framework
cmake_minimum_required(VERSION 3.10)
project(AlgebraicTopologyFramework LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Find Dependencies ---
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/../eigen"
    )
endif()
find_package(GTest QUIET)
if(NOT GTest_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG        release-1.12.1
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    add_library(GTest::Main ALIAS gtest_main)
endif()

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
    test/test_unit_math.cpp
    test/test_unit_concepts.cpp
    test/test_unit_operators.cpp
    test/test_integration_homology.cpp
    test/test_integration_persistence.cpp
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
