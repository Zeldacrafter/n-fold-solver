set(NAME_SRC ${PROJECT_NAME}_src)

set(HEADER_FILES
        template.h)

set(SOURCE_FILES
        nfold_class.cpp
        bruteforce.cpp
        solver_class.cpp
        static_nfold_class.cpp
        static_solver_class.cpp)

add_library(${NAME_SRC} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${NAME_SRC} Eigen3::Eigen)
target_link_libraries(${NAME_SRC} tsl::hopscotch_map)
