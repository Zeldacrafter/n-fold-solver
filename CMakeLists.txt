cmake_minimum_required(VERSION 4.0)
project(n-fold-solver)

#[[
Specify the parameters (n,r,s,t) for which the program should be compiled.
Either specify an input file with the first line consisting of 4 integers 'n r s t'
or specify them separately in the below variables.
]]
set(INPUT_FILE input/input_1.in)
#set(PARAM_N 1)
#set(PARAM_R 2)
#set(PARAM_S 3)
#set(PARAM_T 4)

set(CMAKE_CXX_STANDARD 17)

add_executable(n-fold-solver
        src/main.cpp
        src/nfold_class.hh
        src/prefix_tree_class.hh
        src/solver_class.hh
        src/utils.hh
)

###############
# Compilation #
###############

add_compile_options(-Wall -Wextra -Wpedantic -Wshadow)
add_compile_options(-O2 -g)

################
# Dependencies #
################

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (n-fold-solver Eigen3::Eigen)

add_subdirectory(third-party/hopscotch-map)
target_link_libraries (n-fold-solver tsl::hopscotch_map)

##########################
# Compilation Parameters #
##########################

if(NOT DEFINED PARAM_N OR NOT DEFINED PARAM_R OR NOT DEFINED PARAM_S OR NOT DEFINED PARAM_T)
    get_filename_component(INPUT_FILE_REAL_PATH "${INPUT_FILE}" REALPATH BASE_DIR "${CMAKE_SOURCE_DIR}")
    message(STATUS "No values for 'PARAM_N', 'PARAM_R', 'PARAM_S' and 'PARAM_T' found. Checking content of input file '${INPUT_FILE_REAL_PATH}'")
    if(NOT DEFINED INPUT_FILE OR NOT EXISTS ${INPUT_FILE_REAL_PATH})
        message(FATAL_ERROR "Could not find values for parameters n, r, s and t. Either specify an input file in 'INPUT_FILE' or directly set 'PARAM_N', 'PARAM_R', 'PARAM_S' and 'PARAM_T'")
    endif()
    if(DEFINED INPUT_FILE AND EXISTS "${INPUT_FILE_REAL_PATH}")
        # Extracting parameters (n, r, s, t) from file.
        file(STRINGS ${INPUT_FILE_REAL_PATH} THIS_FILE)
        list(LENGTH THIS_FILE THIS_FILE_LENGTH)
        if(THIS_FILE_LENGTH EQUAL 0)
            message(FATAL_ERROR "Found input file at ${INPUT_FILE_REAL_PATH} but file is empty.")
        endif()
        list(GET THIS_FILE 0 THIS_FILE_LINE_1)
        if (NOT ${THIS_FILE_LINE_1} MATCHES "^[0-9]+ [0-9]+ [0-9]+ [0-9]+$")
            message(FATAL_ERROR "Expected 4 numbers as input parameters (n,r,s,t) in the first line but got '${THIS_FILE_LINE_1}' instead.")
        endif ()
        string(REPLACE " " ";" PARAM_LIST ${THIS_FILE_LINE_1})
        list(GET PARAM_LIST 0 PARAM_N)
        list(GET PARAM_LIST 1 PARAM_R)
        list(GET PARAM_LIST 2 PARAM_S)
        list(GET PARAM_LIST 3 PARAM_T)
    endif()
endif()

if(NOT DEFINED PARAM_N OR NOT DEFINED PARAM_R OR NOT DEFINED PARAM_S OR NOT DEFINED PARAM_T)
    message(FATAL_ERROR "Could not find values for parameters n, r, s and t. Either specify an input file in 'INPUT_FILE' or directly set 'PARAM_N', 'PARAM_R', 'PARAM_S' and 'PARAM_T'")
endif()
message(STATUS "Found size parameters n=${PARAM_N}, r=${PARAM_R}, s=${PARAM_S} and t=${PARAM_T} for compilation.")
add_compile_definitions(-DN_NFOLD=${PARAM_N} -DR_NFOLD=${PARAM_R} -DS_NFOLD=${PARAM_S} -DT_NFOLD=${PARAM_T})