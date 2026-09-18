cmake_minimum_required(VERSION 3.14)
project(BorisPusher CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(HDF5 REQUIRED COMPONENTS CXX)

add_executable(session_2_prequeijo src/main.cpp)
target_include_directories(session_2_prequeijo PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(session_2_prequeijo PRIVATE ${HDF5_CXX_LIBRARIES})

include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.4
)
FetchContent_MakeAvailable(Catch2)

add_executable(session_2_prequeijo_test tests/test_boris_pusher.cpp)
target_include_directories(session_2_prequeijo_test PRIVATE src)
target_link_libraries(session_2_prequeijo_test PRIVATE Catch2::Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(session_2_prequeijo_test)
