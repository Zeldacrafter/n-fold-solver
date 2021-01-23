# Add source to this project's executable.

include(../cmake/GoogleBenchmark.cmake)

set(BENCH_NAME ${PROJECT_NAME}_bench)
add_executable(${BENCH_NAME} example_bench.cc)
target_link_libraries(${BENCH_NAME} benchmark::benchmark)
