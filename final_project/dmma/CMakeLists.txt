cmake_minimum_required(VERSION 3.12)
project(final_project)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fopenmp -mavx2 -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -mavx2 -Wall -Wextra")

# List of versions
set(VERSIONS original version1 version2 version3 version4)

foreach(VERSION ${VERSIONS})
    add_executable(${VERSION} src/mitm_${VERSION}.c)
    target_link_libraries(${VERSION} ${MPI_C_LIBRARIES})
    target_include_directories(${VERSION} PUBLIC src)
endforeach()

