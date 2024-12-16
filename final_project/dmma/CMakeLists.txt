cmake_minimum_required(VERSION 3.12)
project(final_project)
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
add_executable(mitm_version3 src/mitm_version3.c)
target_link_libraries(mitm_version3 ${MPI_C_LIBRARIES})
target_include_directories(mitm_version3 PUBLIC include)

# Set optimization level to 3
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Enable all warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
