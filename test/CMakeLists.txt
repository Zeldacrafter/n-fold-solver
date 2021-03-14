set(NAME_TEST ${PROJECT_NAME}_test)

include(../cmake/GoogleTest.cmake)
include(GoogleTest)

set(SOURCE_FILES
        solver_test.cpp
        example_test.cpp
        )

add_executable(${NAME_TEST} ${SOURCE_FILES})

target_link_libraries(${NAME_TEST} ${PROJECT_NAME}_src)
target_link_libraries(${NAME_TEST} gtest gtest_main)
target_link_libraries (${NAME_TEST} Eigen3::Eigen)
target_link_libraries(${NAME_TEST} tsl::hopscotch_map)
