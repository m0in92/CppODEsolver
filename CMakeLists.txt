cmake_minimum_required(VERSION 3.22)
project(ode_solvers)

set(CMAKE_CXX_STANDARD 14)

set(BOOST_ROOT "C:/Users/Moin/C++_EXT_LIBS/boost_1_80_0/boost_1_80_0")

find_package(Boost REQUIRED)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif()

add_executable(ode_solvers main.cpp odesolvers.h ode_solvers.cpp root_finder.h root_finder.cpp)
#target_link_libraries(odesolvers ${Boost_FILESYSTEM_LIBRARY})
add_subdirectory(Google_tests)