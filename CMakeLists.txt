cmake_minimum_required(VERSION 3.10)
project(SAIGE)

set(CMAKE_CXX_STANDARD 11)

include_directories(src)

add_executable(SAIGE
    src/main.cpp
    src/utils/parameters.cpp
    src/grm/grm.cpp
    src/glmm/glmm.cpp
    src/association/single_variant.cpp
    src/association/set_based.cpp
    src/ld/ld_matrix.cpp
)
