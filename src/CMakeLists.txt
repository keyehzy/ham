cmake_minimum_required(VERSION 3.18)

add_executable(
  ham
  main.cpp
)

target_link_libraries(ham PRIVATE ham-lib openblas lapacke)

add_library(
  ham-lib
  lattice.cpp
  matrix.cpp
  eigensystem.cpp
  vector.cpp
  geometry.cpp

  ham/eigensystem.h
  ham/lattice.h
  ham/matrix.h
  ham/types.h
  ham/vec.h
  ham/vector.h
  ham/geometry.h
  ham/tightbinding.h
  ham/dos.h
  ham/brillouin_zone.h
)

target_include_directories(ham-lib PUBLIC .)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(BLAS_FOUND)
  find_path(BLAS_INCLUDE_DIRS cblas.h $ENV{BLAS_HOME}/include)
  target_include_directories(ham-lib PUBLIC ${BLAS_INCLUDE_DIRS})
endif(BLAS_FOUND)
