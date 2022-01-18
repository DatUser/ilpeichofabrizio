# SETUP
cmake_minimum_required(VERSION 3.12)
project(main)

#set(OPT_FLAGS
#    "-Ofast -mtune=native -floop-parallelize-all -fgcse-sm"
#)


# Directories containing headers must be specified in the include_directories()
include_directories(
    "src"
)

include(${CMAKE_BINARY_DIR}/conanbuildinfo_multi.cmake)
conan_basic_setup()

# Searches for all .cc files
file(GLOB_RECURSE SRC "src/*.cc")
file(GLOB_RECURSE BINDINGS "bindings/*.cpp")

#find_library(SOIL2 SOIL2 REQUIRED)

set(CMAKE_CXX_STANDARD 17)

if (UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_FLAGS} -g -Wall -Wextra -pedantic -fopenmp -lGL -lGLU -lglut -lGLEW")#Werror
endif(UNIX)
if (APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_FLAGS} -Wall -Wextra -pedantic")#Werror
endif(APPLE)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3")# -fsanitize=address")

add_executable(main)
target_sources(main PRIVATE "src/pathtrace.cc" ${SRC} ${BINDINGS})
target_compile_definitions(main PUBLIC IMGUI_IMPL_OPENGL_LOADER_GLEW)

conan_target_link_libraries(main)

#target_link_libraries(main ${SOIL2})