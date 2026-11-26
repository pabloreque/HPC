cmake_minimum_required(VERSION 3.14)
project(PIC CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(HDF5 REQUIRED COMPONENTS CXX)

add_executable(pic_prequeijo src/main.cpp)
target_include_directories(pic_prequeijo PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(pic_prequeijo PRIVATE ${HDF5_CXX_LIBRARIES})

include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.4
)
FetchContent_MakeAvailable(Catch2)

add_executable(pic_prequeijo_test tests/test_boris_pusher.cpp
                                  tests/test_input.cpp
                                  tests/test_periodic.cpp
                                  tests/test_shape.cpp
                                  tests/test_loader.cpp
                                  tests/test_field_solver.cpp
                                  tests/test_hybrid_step.cpp)
target_include_directories(pic_prequeijo_test PRIVATE src tests)
target_link_libraries(pic_prequeijo_test PRIVATE Catch2::Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(pic_prequeijo_test)

# Python module pic_prequeijo (pybind11), built next to the executables.
# Needs the Python headers; configure with the project venv active (or
# -DPython_EXECUTABLE=...) so that pytest, numpy and h5py are found.
option(PIC_BUILD_PYTHON "Build the pic_prequeijo Python module" ON)
if(PIC_BUILD_PYTHON)
    set(PYBIND11_FINDPYTHON ON)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v3.1.0
    )
    FetchContent_MakeAvailable(pybind11)

    pybind11_add_module(pic_prequeijo_python python/bindings.cpp)
    set_target_properties(pic_prequeijo_python PROPERTIES
                          OUTPUT_NAME pic_prequeijo)
    target_include_directories(pic_prequeijo_python
                               PRIVATE src ${HDF5_INCLUDE_DIRS})
    target_link_libraries(pic_prequeijo_python PRIVATE ${HDF5_CXX_LIBRARIES})

    # The pytest suite checks the bindings against the C++ driver.
    add_test(NAME python_bindings
             COMMAND ${Python_EXECUTABLE} -m pytest -q
                     ${CMAKE_CURRENT_SOURCE_DIR}/python)
    set_tests_properties(python_bindings PROPERTIES ENVIRONMENT
                         "PIC_BUILD_DIR=${CMAKE_CURRENT_BINARY_DIR}")
endif()
