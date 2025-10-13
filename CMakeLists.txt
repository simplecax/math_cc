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
# Your code uses Eigen for sparse matrices.
find_package(Eigen3 REQUIRED)

# --- Define Executable Target ---
# Create an executable named 'atopo_app' from the source file in src/
add_executable(atopo_app src/main.cpp)

# --- Configure Target ---
# Add our 'inc' directory to the include path so #include "atopo.h" works.
# Also add the Eigen include directory.
target_include_directories(atopo_app PRIVATE
    inc
    \${EIGEN3_INCLUDE_DIR}
)

# Link the executable against Eigen.
# For header-only libraries, this is the modern CMake way to propagate
# usage requirements (like include paths) to the target.
target_link_libraries(atopo_app PRIVATE Eigen3::Eigen)


# --- Build Information ---
message(STATUS "Project Name: \${PROJECT_NAME}")
message(STATUS "C++ Standard: \${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type: \${CMAKE_BUILD_TYPE}")
message(STATUS "Found Eigen3: \${EIGEN3_INCLUDE_DIR}")
