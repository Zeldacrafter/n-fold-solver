include(../cmake/GoogleTest.cmake)

################################
# Configure Tests
################################

set(TEST_NAME ${PROJECT_NAME}_tests)
add_executable(${TEST_NAME} all_tests.cc)
target_link_libraries(${TEST_NAME} gtest_main)

# Register Tests
include(GoogleTest)
enable_testing()
gtest_discover_tests(${TEST_NAME})
