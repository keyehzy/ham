cmake_minimum_required(VERSION 3.18)

add_executable(
  ham
  main.cpp
)

target_link_libraries(ham PRIVATE ham-lib)

add_library(
  ham-lib
  foo.cpp

  ham/vec.h
  vec.cpp
)

target_include_directories(ham-lib PUBLIC .)
