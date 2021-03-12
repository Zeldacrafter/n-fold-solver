cmake_minimum_required(VERSION 3.14)
project(nFold)

############################
# Configure Compilation
############################

OPTION(MY_DEBUG OFF)
OPTION(MY_LOCAL OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_options(-Wall -Wextra -pedantic -Wno-ignored-qualifiers)
if (MY_DEBUG)
    add_compile_options(-g -fsanitize=address -fsanitize=undefined)
    add_link_options(-fsanitize=address -fsanitize=undefined) # Pass -g too?
else()
    add_compile_options(-O2)
endif(MY_DEBUG)

if (MY_LOCAL)
    add_compile_options(-I /usr/include/eigen3)
endif(MY_LOCAL)

unset(MY_DEBUG CACHE)
unset(MY_LOCAL CACHE)

############################
# Choose what to build
############################

OPTION(BUILD_SRC OFF)
OPTION(BUILD_TEST OFF)
OPTION(BUILD_BENCH OFF)

if (BUILD_SRC)
    add_subdirectory("src")
endif()

if (BUILD_TEST)
    add_subdirectory("test")
endif()

if (BUILD_BENCH)
    add_subdirectory("bench")
endif()
