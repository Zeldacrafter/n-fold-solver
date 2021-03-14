set(NAME_BENCH ${PROJECT_NAME}_bench)

include(../cmake/GoogleBenchmark.cmake)

set(SOURCE_FILES
        example_bench.cpp
        solution_bench.cpp
        )

add_executable(${NAME_BENCH} ${SOURCE_FILES})
target_link_libraries(${NAME_BENCH} ${PROJECT_NAME}_src)
target_link_libraries(${NAME_BENCH} benchmark::benchmark)
