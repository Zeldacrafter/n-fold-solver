set(NAME_SRC ${PROJECT_NAME}_src)

set(HEADER_FILES
        utils.hh
        nfold_class.hh
        solver_class.hh
        split_tree_class.hh)

set(SOURCE_FILES main.cpp)

add_library(${NAME_SRC} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${NAME_SRC} Eigen3::Eigen)
target_link_libraries(${NAME_SRC} tsl::hopscotch_map)
