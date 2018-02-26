cmake_minimum_required(VERSION 2.8.8)

set(PROJECT_NAME_STR adyar)
PROJECT(${PROJECT_NAME_STR} C CXX)
#TODO CHANGE THIS OPTION
set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_BUILD_TYPE Release)

# SDSL directory
set(SDSL_SOURCE_DIR "${PROJECT_SOURCE_DIR}/ext/sdsl-lite/include")

# flags
set(CMAKE_CXX_FLAGS                "-Wall -Wextra -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# include the sdsl-lite and googletest library
add_subdirectory("${PROJECT_SOURCE_DIR}/ext/sdsl-lite")

# add SDSL includes path
include_directories("${CMAKE_BINARY_DIR}/ext/sdsl-lite/include")
include_directories("${CMAKE_BINARY_DIR}/ext/sdsl-lite/external/libdivsufsort/include/")

set(PROJECT_EXEC ${PROJECT_NAME_STR}.x)
ADD_EXECUTABLE(${PROJECT_EXEC} src/main.cpp)

#Link SDSL library
target_link_libraries(${PROJECT_EXEC} sdsl divsufsort divsufsort64)


configure_file(${PROJECT_SOURCE_DIR}/data/text.txt . COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/data/input.txt . COPYONLY)
