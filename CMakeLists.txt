cmake_minimum_required (VERSION 3.8)
project (rt)

# get rid of annoying MSVC warnings.
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

set(CMAKE_CXX_COMPILER "g++")
set (CMAKE_CXX_STANDARD 11)

add_executable(rt main.cpp vect.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
target_link_libraries(${PROJECT_NAME} pthread)
target_link_libraries(rt ${ALL_LIBS} )

