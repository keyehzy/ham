cmake_minimum_required(VERSION 3.18)

add_executable(
  ham
  main.cpp
)

target_link_libraries(ham PRIVATE ham-lib lapacke)

add_library(
  ham-lib
  lattice.cpp

  ham/eigensystem.h
  ham/lattice.h
  ham/matrix.h
  ham/types.h
  ham/vec.h
)

target_include_directories(ham-lib PUBLIC .)
