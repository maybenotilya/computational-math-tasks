cmake_minimum_required(VERSION 3.9)
project(OpenMP-task)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-Wall -Wextra -Wpedantic -O3)

add_executable(task-1 src/main.cpp src/grid.cpp src/funcs.cpp src/algo.cpp)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(task-1 PUBLIC OpenMP::OpenMP_CXX)
endif()