cmake_minimum_required(VERSION 3.10)
project(FiniteDifferenceOrder CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(HDF5 REQUIRED COMPONENTS CXX)

add_executable(fd_order src/main.cpp)
target_include_directories(fd_order PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(fd_order PRIVATE ${HDF5_CXX_LIBRARIES})
