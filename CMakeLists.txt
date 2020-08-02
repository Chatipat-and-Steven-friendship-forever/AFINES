cmake_minimum_required(VERSION 3.11)

project(AFINES)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Boost 1.53 REQUIRED COMPONENTS filesystem program_options system)

set(sources
    src/bead.cpp
    src/spring.cpp
    src/filament.cpp
    src/filament_ensemble.cpp
    src/motor.cpp
    src/motor_ensemble.cpp
    src/potentials.cpp
    src/exv.cpp
    src/ext.cpp
    src/box.cpp
    src/quadrants.cpp
    src/generate.cpp
    src/globals.cpp
)

add_subdirectory(external/fmt)
add_subdirectory(external/pybind11)

add_executable(network prog/network.cpp ${sources})
target_include_directories(network PRIVATE include external/pcg-cpp/include)
target_link_libraries(network PRIVATE Boost::filesystem Boost::program_options Boost::system fmt::fmt-header-only)

pybind11_add_module(pyafines prog/pyafines.cpp ${sources})
target_include_directories(pyafines PUBLIC include external/pcg-cpp/include)
target_link_libraries(pyafines PRIVATE Boost::filesystem Boost::program_options Boost::system fmt::fmt-header-only)
