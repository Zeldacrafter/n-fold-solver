# Required by FetchContent_MakeAvailable
cmake_minimum_required(VERSION 3.14)
project(nFold)

set(MY_DEBUG True)

############################
# Configure Compilation
############################

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_options(-Wall -Wextra -pedantic)
if (MY_DEBUG)
    add_compile_options(-g -fsanitize=address -fsanitize=undefined)
else()
    add_compile_options(-02)
endif()

############################
# External Dependencies
############################

include(FetchContent)

# GTest
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.10.0
)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)
enable_testing()

############################
# Configure Tests
############################

set(TEST_NAME ${PROJECT_NAME}_tests)
add_executable(${TEST_NAME} test/example_test.cc)
target_link_libraries(${TEST_NAME} gtest_main)
if (MY_DEBUG)
    target_link_libraries(${TEST_NAME} -fsanitize=address -fsanitize=undefined)
endif()
gtest_discover_tests(${TEST_NAME})

