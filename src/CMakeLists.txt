set(NAME_SRC ${PROJECT_NAME}_src)

set(HEADER_FILES
        template.h)

set(SOURCE_FILES
        NFold.cpp
        bruteforce.cpp
        Solver.cpp)

add_library(${NAME_SRC} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${NAME_SRC} Eigen3::Eigen)
