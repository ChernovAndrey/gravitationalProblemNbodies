cmake_minimum_required(VERSION 3.8)
project(gravitationalProblem)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp Particle.h RK.h SolveMethod.h vectorOperations.cpp vectorOperations.h)
add_executable(gravitationalProblem ${SOURCE_FILES})