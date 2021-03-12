set(NAME_SRC ${PROJECT_NAME}_src)

set(HEADER_FILES
        template.h)

set(SOURCE_FILES
        NFold.cpp bruteforce.cpp)

if(CUSTOM_EIGEN_IMPORT)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(${NAME_MAIN} Eigen3::Eigen)
endif(CUSTOM_EIGEN_IMPORT)
#find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(${NAME_SRC} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${NAME_SRC} Eigen3::Eigen)
