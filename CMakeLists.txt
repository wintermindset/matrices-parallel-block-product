cmake_minimum_required(VERSION 3.10.0)

set(PROJECT_NAME matrices_parallel_block_product)
project(${PROJECT_NAME} C)

enable_testing()

add_executable(matrices_parallel_block_product src/main.c)
add_executable(tests tests/tests.c)

add_test(NAME TEST_BLOCK_PRODUCT COMMAND tests)