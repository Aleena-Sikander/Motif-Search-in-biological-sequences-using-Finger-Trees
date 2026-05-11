cmake_minimum_required(VERSION 3.14)
project(FingerTreeBio CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

add_executable(fingertree src/main.cpp)
target_include_directories(fingertree PRIVATE src)
target_compile_options(fingertree PRIVATE -Wall -Wextra -Wpedantic)
