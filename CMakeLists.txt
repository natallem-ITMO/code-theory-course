cmake_minimum_required(VERSION 3.20)
project(Lab1)

set(CMAKE_CXX_STANDARD 17)

option(USE_DEBUG "Enter debug mode" OFF)
if (USE_DEBUG)
    add_definitions(-DDEBUG)
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
add_executable(Lab1
#        test.cpp
        Lab1/main.cpp
        )
