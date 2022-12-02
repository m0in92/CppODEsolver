cmake_minimum_required(VERSION 3.22)
project(ode_solvers)

set(CMAKE_CXX_STANDARD 14)

add_executable(ode_solvers main.cpp odesolvers.h ode_solvers.cpp)
add_subdirectory(Google_tests)