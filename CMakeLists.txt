cmake_minimum_required(VERSION 3.10.0)

set(PROJECT_NAME matrices_parallel_block_product)
project(${PROJECT_NAME} C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -O2")

include_directories(src)

# add_executable(${PROJECT_NAME} src/main.c src/matrix_ops.c)
# set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")

add_executable(tests tests/tests.c src/matrix_ops.c)
set_target_properties(tests PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")

enable_testing()
add_test(NAME TEST_BLOCK_PRODUCT COMMAND tests)