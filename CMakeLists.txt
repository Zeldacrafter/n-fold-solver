cmake_minimum_required(VERSION 3.10.2)
project(nFold)

OPTION(MY_DEBUG OFF)
OPTION(CUSTOM_EIGEN_IMPORT OFF)
OPTION(BUILD_SRC OFF)
OPTION(BUILD_TEST OFF)
OPTION(BUILD_BENCH OFF)

############################
# Configure Compilation    #
############################

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(NAME_MAIN ${PROJECT_NAME}_main)

# Boost stacktrace
add_compile_options(-DBOOST_STACKTRACE_USE_BACKTRACE -g -rdynamic -ldl -lbacktrace)
add_link_options(-ldl -lbacktrace)

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

include_directories(src)
if(CUSTOM_EIGEN_IMPORT)
    cmake_policy(SET CMP0079 NEW)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(${NAME_MAIN} Eigen3::Eigen)
endif(CUSTOM_EIGEN_IMPORT)

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

