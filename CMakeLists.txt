cmake_minimum_required(VERSION 3.5)

set (CMAKE_CONFIGURATION_TYPE x64)
set (CMAKE_BUILD_TYPE DEBUG) # or Debug
set (CMAKE_CXX_FLAGS "-Wall -std=c++11 -fPIC")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage -DBZ_DEBUG")
# set (CMAKE_MODULE_LINKER_FLAGS  "${CMAKE_MODULE_LINKER_FLAGS} -Wl,-soname,pyblitzdg.so")
include_directories ("include" "/usr/include/python3.6m/")
set (CMAKE_CXX_LINK_FLAGS "-Wl,-soname=pyblitzdg.so")


link_libraries(
    blitz
    metis
    umfpack
    cxsparse
    lapack
    blas
    vtkIOXML-7.1
    vtkCommonCore-7.1
    vtkCommonExecutionModel-7.1
    vtkCommonDataModel-7.1
    python3.6m
    boost_python3
    boost_numpy36
)

project(blitzdg)
file (GLOB sources src/*.cpp)
add_library(blitzdg SHARED ${sources})
# target_compile_options(blitzdg PUBLIC "-fPIC")
set_target_properties(blitzdg
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )


# link all projets below to blitzdg.
link_libraries(
    blitz
    metis
    umfpack
    cxsparse
    lapack
    blas
    vtkIOXML-7.1
    vtkCommonCore-7.1
    vtkCommonExecutionModel-7.1
    vtkCommonDataModel-7.1
    python3.6m
    boost_python3
    boost_numpy36
    blitzdg
)

project(advec1d)
file (GLOB sources src/advec1d/*.cpp)
add_executable(advec1d ${sources})
set_target_properties(advec1d
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

file (GLOB sources src/burgers1d/*.cpp)
add_executable(burgers1d ${sources})
set_target_properties(burgers1d
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

project(poisson1d)
file (GLOB sources src/poisson1d/*.cpp)
add_executable(poisson1d ${sources})
set_target_properties(poisson1d
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

project(sw2d)
file (GLOB sources src/sw2d/*.cpp)
add_executable(sw2d ${sources})
set_target_properties(sw2d
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

project(sw2d-simple)
file (GLOB sources src/sw2d-simple/*.cpp)
add_executable(sw2d-simple ${sources})
set_target_properties(sw2d-simple
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

project(pyblitzdg)

file (GLOB sources src/pyblitzdg/*.cpp)
add_library(pyblitzdg SHARED ${sources})
# target_compile_options(pyblitzdg PUBLIC "-fPIC")
set_target_properties(pyblitzdg
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
set_target_properties(pyblitzdg PROPERTIES PREFIX "")
