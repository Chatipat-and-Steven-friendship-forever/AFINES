cmake_minimum_required(VERSION 3.11)

project(AFINES)

find_package(Boost 1.53 REQUIRED COMPONENTS filesystem program_options system)
find_package(fmt REQUIRED)
find_package(pybind11 REQUIRED)

set(sources
    src/box.cpp
    src/ext.cpp
    src/exv.cpp
    src/filament.cpp
    src/filament_ensemble.cpp
    src/generate.cpp
    src/globals.cpp
    src/motor.cpp
    src/motor_ensemble.cpp
    src/potentials.cpp
    src/quadrants.cpp
    src/shear.cpp
    src/spring.cpp
)

add_executable(network prog/network.cpp ${sources})
target_include_directories(network PRIVATE include external/pcg-cpp/include)
target_link_libraries(network PRIVATE Boost::filesystem Boost::program_options Boost::system fmt::fmt)

pybind11_add_module(pyafines prog/pyafines.cpp ${sources})
target_include_directories(pyafines PRIVATE include external/pcg-cpp/include)
target_link_libraries(pyafines PRIVATE Boost::filesystem Boost::system fmt::fmt)
