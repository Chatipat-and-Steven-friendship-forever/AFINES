cmake_minimum_required(VERSION 3.11)

project(AFINES)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost 1.53 REQUIRED COMPONENTS filesystem program_options system)

set(sources src/bead.cpp src/spring.cpp src/filament.cpp src/filament_ensemble.cpp src/motor.cpp src/motor_ensemble.cpp src/box.cpp src/quadrants.cpp src/generate.cpp src/globals.cpp)

add_executable(network prog/network.cpp ${sources})
target_include_directories(network PRIVATE include)
target_link_libraries(network PRIVATE Boost::filesystem Boost::program_options Boost::system)

add_subdirectory(external/pybind11)
pybind11_add_module(pyafines prog/pyafines.cpp ${sources})
target_include_directories(pyafines PUBLIC include)
target_link_libraries(pyafines PRIVATE Boost::filesystem Boost::program_options Boost::system)
