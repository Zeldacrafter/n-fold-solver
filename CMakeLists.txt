cmake_minimum_required(VERSION 3.17)
project(nFold)

OPTION(MY_DEBUG OFF)
OPTION(BUILD_SRC OFF)
OPTION(BUILD_TEST OFF)
OPTION(BUILD_BENCH OFF)

############################
# Configure Compilation    #
############################

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(NAME_MAIN ${PROJECT_NAME}_main)

add_compile_options(-Wall -Wextra -pedantic -Wno-ignored-qualifiers)
if (MY_DEBUG)
    add_compile_options(-g -fsanitize=address -fsanitize=undefined)
    add_link_options(-fsanitize=address -fsanitize=undefined)
else()
    add_compile_options(-O2)
endif(MY_DEBUG)

############################
# Configure Main Module    #
############################

set(SOURCE_FILES main.cpp)
add_executable(${NAME_MAIN} ${SOURCE_FILES})

# Add dependencies for main module
include_directories(src)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(${NAME_MAIN} Eigen3::Eigen)

############################
# Choose what to build     #
############################

add_subdirectory("src")

if (BUILD_TEST)
    add_subdirectory("test")
endif()

if (BUILD_BENCH)
    add_subdirectory("bench")
endif()

