cmake_minimum_required(VERSION 3.5)

set (CMAKE_CONFIGURATION_TYPES "Win64")
set (CMAKE_BUILD_TYPE Release) # or Debug
set (CMAKE_CXX_FLAGS "-std=c++11 -fPIC")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage -DBZ_DEBUG")
# set (CMAKE_MODULE_LINKER_FLAGS  "${CMAKE_MODULE_LINKER_FLAGS}")
include_directories ("include" "include/igloo" "/usr/include/python3.7m/" "/usr/local/include/" "/usr/local/Cellar/python/3.7.4_1/Frameworks/Python.framework/Versions/3.7/include/python3.7m/")


# include path: "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Tools\MSVC\14.20.27508\include\vcruntime.h"
# Link to vtk using static mode (otherwise won't work on windows)
add_compile_definitions(VTKCOMMONCORE_STATIC_DEFINE 
  VTKCOMMONEXECUTIONMODEL_STATIC_DEFINE 
  VTKIOGEOMETRY_STATIC_DEFINE 
  VTKCOMMONDATAMODEL_STATIC_DEFINE 
  VTKIOXML_STATIC_DEFINE
  VTKRENDERINGCORE_STATIC_DEFINE
  )

set(CMAKE_FIND_LIBRARY_PREFIXES "" "lib")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".a"  ".so")
message("--looking in lib")
find_library(BLITZ NAMES "blitz" "libblitz" PATHS "lib")
find_library(METIS NAMES "metis" "libmetis" PATHS "lib")
find_library(UMFPACK NAMES "umfpack" "libumfpack" PATHS "lib")
find_library(CXSPARSE NAMES "cxsparse" "libcxsparse" PATHS "lib")
find_library(AMD NAMES "amd" "libamd" PATHS "lib")
find_library(CHOLMOD names "cholmod" "libcholmod" PATHS "lib")
find_library(COLAMD names "colamd" "libcolamd" PATHS "lib")
find_library(SUITESPARSECONFIG NAMES "suitesparseconfig" "libsuitesparseconfig" PATHS "lib")

find_library(LAPACK NAMES "liblapack" "lapack" PATHS "lib")
find_library(BLAS NAMES "libblas" "blas" PATHS "lib")

find_library(PYTHON37 NAMES "python37" "libpython37" PATHS "lib") 
find_library(BOOSTPYTHON NAMES "boost_python3-vc140-mt-1_60" PATHS "lib")
find_library(BOOSTNUMPY NAMES "boost_numpy3-vc140-mt-1_60" PATHS "lib")

message("-- found ${BOOSTNUMPY}")
message("-- found ${BOOSTPYTHON}")
message("-- found ${PYTHON37}")
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
    python3.7m
    boost_python3
    boost_numpy3
)

project(blitzdg)
file (GLOB sources src/*.cpp)
add_library(blitzdg SHARED ${sources})
# target_compile_options(blitzdg PUBLIC "-fPIC")
set_target_properties(blitzdg
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    LINK_FLAGS "-Wl,-rpath,/usr/local/lib -L /usr/local/lib -L/usr/local/Cellar/python/3.7.4_1/Frameworks/Python.framework/Versions/3.7/lib/"
    )

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
        python3.7m
        boost_python3
        boost_numpy3
        blitzdg
    )

project(pyblitzdg)

file (GLOB sources src/pyblitzdg/*.cpp)
add_library(pyblitzdg SHARED ${sources})
set_target_properties(pyblitzdg
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    LINK_FLAGS "-Wl,-rpath,/usr/local/lib -L /usr/local/lib -L/usr/local/Cellar/python/3.7.4_1/Frameworks/Python.framework/Versions/3.7/lib/"
    )
set_target_properties(pyblitzdg PROPERTIES PREFIX "")

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
    python3.7m
    boost_python3
    boost_numpy3
    blitzdg
    pyblitzdg
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

project(tests)
file (GLOB sources src/test/*.cpp)
add_executable(tests ${sources})
set_target_properties(tests
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
