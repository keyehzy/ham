cmake_minimum_required(VERSION 3.18)

add_executable(
  ham-test
  test-foo.cpp
)

target_include_directories(ham-test PRIVATE .)
target_link_libraries(
  ham-test
  PRIVATE
  gmock
  gmock_main
  gtest
  ham-lib
)

add_test(
  NAME ham-test
  COMMAND ham-test
)
