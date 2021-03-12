set(NAME_BENCH ${PROJECT_NAME}_bench)

include(../cmake/GoogleBenchmark.cmake)

set(SOURCE_FILES
        example_bench.cpp
        )

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(${NAME_MAIN} Eigen3::Eigen)

add_executable(${NAME_BENCH} ${SOURCE_FILES})
target_link_libraries(${NAME_BENCH} benchmark::benchmark)
