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
set (CMAKE_CXX_FLAGS "-Wall -std=c++11 -fPIC")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2 -L/usr/lib/vtk")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage -DBZ_DEBUG")
include_directories ("include" "include/igloo" "/usr/include/python3.6m/")

set (pyblitzdg_VERSION_MAJOR 0)
set (pyblitzdg_VERSION_MINOR 1)
set (pyblitzdg_VERSION_PATCH 2)
set (DISABLE_SOLVER_BUILDS on)

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
set_target_properties(blitzdg
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

project(pyblitzdg)

file (GLOB sources src/pyblitzdg/*.cpp)
add_library(pyblitzdg SHARED ${sources})
set_target_properties(pyblitzdg
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
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
    python3.6m
    boost_python3
    boost_numpy3
    blitzdg
    pyblitzdg
)

if(NOT DISABLE_SOLVER_BUILDS)
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

    INSTALL(TARGETS sw2d-simple
        DESTINATION "/usr/local/bin"
    )
endif()

project(tests)
file (GLOB sources src/test/*.cpp)
add_executable(tests ${sources})
set_target_properties(tests
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

SET(CPACK_PACKAGE_DESCRIPTION "Discontinuous Galerkin library and solvers")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Discontinuous Galerkin library and solvers.")
SET(CPACK_RPM_PACKAGE_SUMMARY "Discontinuous Galerkin library and solvers.")
SET(CPACK_PACKAGE_VENDOR "dsteinmo@wqcg.ca")
SET(CPACK_PACKAGE_CONTACT "dsteinmo@wqcg.ca")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

set (CPACK_PACKAGE_VERSION "0.1.2")

INSTALL(TARGETS pyblitzdg blitzdg
    DESTINATION "/usr/lib/x86_64-linux-gnu/"
)
INSTALL(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/pyblitzdg/libblitz.so.0"
    "${CMAKE_CURRENT_SOURCE_DIR}/pyblitzdg/libblitzdg.so"
    DESTINATION "/usr/lib/x86_64-linux-gnu/"
)

INSTALL(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/advec1d.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/sw2d.py"
    DESTINATION "/usr/local/bin/"
)

SET(CPACK_GENERATOR DEB)
SET(CPACK_PACKAGE_NAME "pyblitzdg")
SET(CPACK_SOURCE_STRIP_FILES "")
SET(CPACK_DEBIAN_PACKAGE_BUILDS_DEPENDS "build-essential,cmake,libvtk7.1,libvtk7-dev,libsuitesparse-dev,libmetis-dev,libboost-dev,libboost-python1.65.1-dev,libboost-numpy1.65.1-dev")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libvtk7.1,python3-pip,python3-numpy,libsuitesparse-dev,libmetis-dev,libboost-dev,libboost-python1.65.1,libboost-numpy1.65.1")
SET(CPACK_DEBIAN_FILE_NAME "DEB-DEFAULT")
SET(CPACK_PACKAGE_CONTACT "Derek Steinmoeller <dsteinmo@wqcg.ca>")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Derek Steinmoeller <dsteinmo@wqcg.ca>")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA
    "${CMAKE_CURRENT_SOURCE_DIR}/deb/postinst" "${CMAKE_CURRENT_SOURCE_DIR}/deb/postrm")

INCLUDE(CPack)
IF(DEBIAN_FOUND)
  message("-- debian found!")
  ADD_DEBIAN_TARGETS(pyblitzdg)
ENDIF(DEBIAN_FOUND)
INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_GENERATOR RPM)
set(CPACK_RPM_PACKAGE_VERSION "0.1.2")
SET(CPACK_RPM_PACKAGE_DEPENDS "python3,python3-pip,libgfortran.x86_64,suitesparse-devel.x86_64,metis-devel.x86_64,boost-devel.x86_64,vtk-devel.x86_64 vtk.x86_64,boost-python3-devel.x86_64,boost-numpy3.x86_64")
