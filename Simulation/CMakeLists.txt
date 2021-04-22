cmake_minimum_required(VERSION 3.17)
project(Dissertation)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
add_executable(Dissertation main.cpp)
target_link_libraries(Dissertation  OpenMP::OpenMP_CXX)
