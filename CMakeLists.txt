set(OPT_FLAGS
    "-Ofast -mtune=native -floop-parallelize-all -fgcse-sm"
)


# Directories containing headers must be specified in the include_directories()
include_directories(
    "src"
)

# Searches for all .cc files
set(SRC
    "src/sphere.cc"
    "src/camera.cc"
    "src/matrix.cc"
    "src/vector3.cc"
    "src/program.cc"
    "src/main.cc")


# SETUP
cmake_minimum_required(VERSION 3.12)
project(main)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_FLAGS} -Wall -Werror -Wextra -pedantic -fopenmp -lGL -lGLU -lglut -lGLEW")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3")# -fsanitize=address")

if (NOT CMAKE_BUILD_TYPE)
    message("No build type selected, default release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Test")
    message("Test Mode")
    add_executable(main)
    target_sources(main PRIVATE "tests/main.cc" ${SRC})
    target_link_libraries(main PRIVATE ${LIBRARIES})
else()
    add_executable(main)
    target_sources(main PRIVATE "src/main.cc" ${SRC})
    target_link_libraries(main PRIVATE ${LIBRARIES})
endif()
