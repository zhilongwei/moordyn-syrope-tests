CMake_minimum_required (VERSION 3.18)
project (myproject)

# Path to the installed MoorDyn CMake package config (MoorDynConfig.cmake)
set(MoorDyn_DIR "/opt/moordyn/lib/cmake/moordyn")

# External dependencies
find_package (MoorDyn REQUIRED)
find_package (Eigen3 REQUIRED)

# Driver executable
add_executable (main main.cpp)

# Link required imported targets (include paths/flags propagate automatically)
target_link_libraries (main PRIVATE Eigen3::Eigen)
target_link_libraries (main PRIVATE MoorDyn::moordyn)