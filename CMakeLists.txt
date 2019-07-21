cmake_minimum_required(VERSION 3.5)

set (CMAKE_CONFIGURATION_TYPE x64)
set (CMAKE_BUILD_TYPE Release) # or Debug
set (CMAKE_CXX_FLAGS "-Wall -std=c++11 -fPIC")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage -DBZ_DEBUG")
include_directories ("include" "include/igloo" "/usr/include/python3.6m/")

set (pyblitzdg_VERSION_MAJOR 0)
set (pyblitzdg_VERSION_MINOR 1)
set (pyblitzdg_VERSION_PATCH 2)
set (DISABLE_SOLVER_BUILDS on)

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
)

project(blitzdg)
file (GLOB sources src/*.cpp)
add_library(blitzdg SHARED ${sources})
set_target_properties(blitzdg
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    LINK_FLAGS "-Wl,-soname=libblitzdg.so -Wl,-rpath=/usr/local/lib/python3.6/dist-packages/pyblitzdg/"
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
        python3.6m
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
    LINK_FLAGS "-Wl,-soname=pyblitzdg.so -Wl,-rpath=/usr/local/lib/python3.6/dist-packages/pyblitzdg/"
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
    "${CMAKE_CURRENT_SOURCE_DIR}/deb/postinst")

INCLUDE(CPack)
IF(DEBIAN_FOUND)
  message("-- debian found!")
  ADD_DEBIAN_TARGETS(pyblitzdg)
ENDIF(DEBIAN_FOUND)
INCLUDE(InstallRequiredSystemLibraries)