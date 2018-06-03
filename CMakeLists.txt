cmake_minimum_required(VERSION 3.8)
project(gravitationalProblem)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")
set(SOURCE_FILES main.cpp SolveMethod.h vectorOperations.cpp vectorOperations.h)
add_executable(gravitationalProblem ${SOURCE_FILES})