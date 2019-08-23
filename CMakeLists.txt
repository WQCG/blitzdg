cmake_minimum_required(VERSION 3.5)

set (CMAKE_CONFIGURATION_TYPES "Win64")
set (CMAKE_BUILD_TYPE Release) # or Debug
# set (CMAKE_CXX_FLAGS "-Wall -std=c++11 -fPIC")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage -DBZ_DEBUG")
# set (CMAKE_MODULE_LINKER_FLAGS  "${CMAKE_MODULE_LINKER_FLAGS} -Wl,-soname,pyblitzdg.so")
include_directories ("include" "include/igloo" "/usr/include/python3.7m/" "C:/Users/dsteinmo/AppData/Local/Programs/Python/Python37/include/" "C:/local/boost_1_66_0" "C:/Users/dsteinmo/Downloads/boost_1_63_0" "C:/Program Files (x86)/Microsoft Visual Studio/2019/BuildTools/VC/Tools/MSVC/14.20.27508/include/" "C:/ProgramData/Anaconda3/pkgs/python-3.7.3-h8c8aaf0_0/include")
set (CMAKE_CXX_LINK_FLAGS "-Wl,-soname=pyblitzdg.so")
set (CMAKE_CXX_STANDARD "14")
set (CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS on)


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
set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".a"  ".so" ".dll")
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
find_library(BOOSTPYTHON NAMES "boost_python37-vc140-mt-x64-1_67" PATHS "lib")
find_library(BOOSTNUMPY NAMES "boost_numpy37-vc140-mt-x64-1_67" PATHS "lib")

message("-- found ${BOOSTNUMPY}")
message("-- found ${BOOSTPYTHON}")
message("-- found ${PYTHON37}")
link_libraries(
    ${BLITZ}
    ${METIS}
    ${UMFPACK}
    ${CXSPARSE}
    ${LAPACK}
    ${BLAS}
    ${CHOLMOD}
    ${AMD}
    ${COLAMD}
    ${SUITESPARSECONFIG}
    ${BOOSTPYTHON}
    ${BOOSTNUMPY}
    #${GFORTRAN}
    #${QUADMATH}
    ${PYTHON37}
    #boost_python3
    #boost_numpy3
    lib/vtkIOXML-7.1
    lib/vtkCommonCore-7.1
    lib/vtkCommonExecutionModel-7.1
    lib/vtkCommonDataModel-7.1
    lib/vtkCommonMisc-7.1
    lib/vtkCommonSystem-7.1
    lib/vtkCommonTransforms-7.1
    lib/vtkexpat-7.1
    lib/vtkIOCore-7.1
    lib/vtkIOGeometry-7.1
    lib/vtkIOXML-7.1
    lib/vtkIOXMLParser-7.1
    lib/vtksys-7.1
    lib/vtkzlib-7.1
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
    LINK_FLAGS "-Wl,-soname=libblitzdg.so -Wl,-rpath=/usr/local/lib/python3.6/dist-packages/pyblitzdg/"
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
