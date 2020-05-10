cmake_minimum_required(VERSION 3.5)

set (CMAKE_CONFIGURATION_TYPES "Release Debug")
set (CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS 1)
set (CMAKE_BUILD_TYPE Release) # or Debug
set (CMAKE_CXX_FLAGS "-std=c++14 -fPIC")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage -DBZ_DEBUG")
# set (CMAKE_MODULE_LINKER_FLAGS  "${CMAKE_MODULE_LINKER_FLAGS}")
include_directories ("include"
    "include/igloo"
    "/opt/_internal/cpython-3.7.3/include/python3.7m/"
    "/miniconda3/pkgs/python-3.7.3-h0371630_0/include/python3.7m/"
    "/usr/include/python3.7m/"
    "/usr/local/include/"
    "/usr/local/Cellar/python/3.7.4_1/Frameworks/Python.framework/Versions/3.7/include/python3.7m/"
    "C:/ProgramData/Anaconda3/pkgs/python-3.7.3-h8c8aaf0_0/include/"
    "C:/ProgramData/Anaconda3/pkgs/py-boost-1.67.0-py37h8300f20_4/Library/include"
    "c:/dev/blitzdg/include"
    )

set(CMAKE_FIND_LIBRARY_PREFIXES "" "lib")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".dylib" ".lib" ".a")
find_library(BLITZ NAMES "blitz" "libblitz" PATHS "/usr/local/lib64" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(METIS NAMES "metis" "libmetis" PATHS "/usr/lib64" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(UMFPACK NAMES "umfpack" "libumfpack" PATHS "/usr/lib64/" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(CXSPARSE NAMES "cxsparse" "libcxsparse" PATHS "/usr/lib64/" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(AMD NAMES "amd" "libamd" PATHS "/usr/lib64/" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(CHOLMOD names "cholmod" "libcholmod" PATHS "/usr/lib64/" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(COLAMD names "colamd" "libcolamd" PATHS "/usr/lib64/" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")

# for static linkage on windows.
if (WIN32)
    find_library(SUITESPARSECONFIG NAMES "suitesparseconfig" "libsuitesparseconfig" PATHS "/usr/lib64/" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
endif()
find_library(LAPACK NAMES "liblapack" "lapack" PATHS "/usr/lib64/" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(BLAS NAMES "libblas" "blas" "libblas" PATHS "/usr/lib64/" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")

find_library(PYTHON3 NAMES "python37" "libpython37" "python3.7" "python3.7m" "libpython3.7" PATHS "/miniconda3/pkgs/python-3.7.3-h0371630_0/lib/" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "/usr/lib/python3.7/config-3.7m-x86_64-linux-gnu" "/usr/lib/python3.6/config-3.6m-x86_64-linux-gnu" "lib" "/usr/local/Cellar/python/3.7.4_1/Frameworks/Python.framework/Versions/3.7/lib/")
find_library(BOOSTPYTHON NAMES "boost_python37-vc140-mt-x64-1_67" "libboost_python37" "boost_python37" "boost_python3.7" "boost_python3.6" PATHS "/boost_1_72/lib/" "/usr/lib64" "C:/ProgramData/Anaconda3/pkgs/py-boost-1.67.0-py37h8300f20_4/Library/lib" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "/usr/lib/python3.7/config-3.7m-x86_64-linux-gnu" "/usr/lib/python3.6/config-3.6m-x86_64-linux-gnu" "lib")
find_library(BOOSTNUMPY NAMES "boost_numpy37-vc140-mt-x64-1_67" "libboost_numpy37" "boost_numpy37" "boost_numpy36" "boost_numpy37" PATHS "/boost_1_72/lib/" "/usr/lib64" "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")

find_library(VTKIOXML7 NAMES "vtkIOXML-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKCOMMONCORE7 NAMES "vtkCommonCore-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKCOMMONEXECUTIONMODEL7 NAMES "vtkCommonExecutionModel-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKCOMMONDATAMODEL7 NAMES "vtkCommonDataModel-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKCOMMONMISC7 NAMES "vtkCommonMisc-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKCOMMONSYSTEM7 NAMES "vtkCommonSystem-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKCOMMONTRANSFORMS7 NAMES "vtkCommonTransforms-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKIOCORE7 NAMES "vtkIOCore-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKIOGEOMETRY7 NAMES "vtkIOGeometry-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKIOXML7 NAMES "vtkIOXML-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKIOXMLPARSER7 NAMES "vtkIOXMLParser-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
find_library(VTKSYS7 NAMES "vtksys-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")

if (WIN32)
    find_library(VTKZLIB7 NAMES "vtkzlib-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
    find_library(VTKEXPAT7 NAMES "vtkexpat-7.1" PATHS "/usr/lib" "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "lib")
    message("-- found ${VTKZLIB7}")
    message("-- found ${VTKEXPAT7}")
endif()

message("-- found ${BLITZ}")
message("-- found ${METIS}")
message("-- found ${UMFPACK}")
message("-- found ${CXSPARSE}")
message("-- found ${AMD}")
message("-- found ${CHOLMOD}")
message("-- found ${COLAMD}")
message("-- found ${SUITESPARSECONFIG}")
message("-- found ${LAPACK}")
message("-- found ${BLAS}")
message("-- found ${PYTHON3}")
message("-- found ${BOOSTNUMPY}")
message("-- found ${BOOSTPYTHON}")
message("-- found ${VTKIOXML7}")
message("-- found ${VTKCOMMONCORE7}")
message("-- found ${VTKCOMMONEXECUTIONMODEL7}")
message("-- found ${VTKCOMMONDATAMODEL7}")
message("-- found ${VTKCOMMONMISC7}")
message("-- found ${VTKCOMMONSYSTEM7}")
message("-- found ${VTKCOMMONTRANSFORMS7}")
message("-- found ${VTKIOCORE7}")
message("-- found ${VTKIOGEOMETRY7}")
message("-- found ${VTKIOXML7}")
message("-- found ${VTKIOXMLPARSER7}")
message("-- found ${VTKSYS7}")

set (CMAKE_CXX_FLAGS "-Wall -std=c++14 -fPIC")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
if (WIN32)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /EHsc /fp:except")
	message("-- enabling exception handling")
endif()
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage -DBZ_DEBUG")
include_directories ("include" 
    "include/igloo" 
    "/opt/_internal/cpython-3.7.3/include/python3.7m/"
    "/miniconda3/pkgs/python-3.7.3-h0371630_0/include/python3.7m/"    
    "/usr/include/python3.7m/"
    "C:/ProgramData/Anaconda3/pkgs/python-3.7.3-h8c8aaf0_0/include/"
    "C:/ProgramData/Anaconda3/pkgs/py-boost-1.67.0-py37h8300f20_4/Library/include"
    "C:/ProgramData/Anaconda3/Library/include"
    "C:/dev/blitzdg/include"
    )

set (pyblitzdg_VERSION_MAJOR 0)
set (pyblitzdg_VERSION_MINOR 3)
set (pyblitzdg_VERSION_PATCH 0)
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
    ${PYTHON3}
    ${VTKIOXML7}
    ${VTKCOMMONCORE7}
    ${VTKCOMMONEXECUTIONMODEL7}
    ${VTKCOMMONDATAMODEL7}
    ${VTKCOMMONMISC7}
    ${VTKCOMMONSYSTEM7}
    ${VTKCOMMONTRANSFORMS7}
    ${VTKIOCORE7}
    ${VTKIOGEOMETRY7}
    ${VTKIOXML7}
    ${VTKIOXMLPARSER7}
    ${VTKSYS7}
    ${VTKZLIB7}
    ${VTKEXPAT7}
)

project(blitzdg)
file (GLOB sources src/*.cpp)
add_library(blitzdg SHARED ${sources})
set_target_properties(blitzdg
    PROPERTIES
    LINK_FLAGS "-Wl,-rpath,/usr/local/lib -L /usr/local/lib -L/usr/local/Cellar/python/3.7.4_1/Frameworks/Python.framework/Versions/3.7/lib/"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Link to vtk using static mode (otherwise won't work on windows)
target_compile_definitions(blitzdg PUBLIC "VTKCOMMONCORE_STATIC_DEFINE"
  "VTKCOMMONEXECUTIONMODEL_STATIC_DEFINE"
  "VTKIOGEOMETRY_STATIC_DEFINE"
  "VTKCOMMONDATAMODEL_STATIC_DEFINE"
  "VTKIOXML_STATIC_DEFINE"
  "VTKRENDERINGCORE_STATIC_DEFINE"
  "PY_MAJOR_VERSION=3"
  "PY_MINOR_VERSION=7"
  )

project(pyblitzdg)

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
    ${PYTHON3}
    ${VTKIOXML7}
    ${VTKCOMMONCORE7}
    ${VTKCOMMONEXECUTIONMODEL7}
    ${VTKCOMMONDATAMODEL7}
    ${VTKCOMMONMISC7}
    ${VTKCOMMONSYSTEM7}
    ${VTKCOMMONTRANSFORMS7}
    ${VTKIOCORE7}
    ${VTKIOGEOMETRY7}
    ${VTKIOXML7}
    ${VTKIOXMLPARSER7}
    ${VTKSYS7}
    ${VTKZLIB7}
    ${VTKEXPAT7}
    blitzdg
)

file (GLOB sources src/pyblitzdg/*.cpp)
add_library(pyblitzdg SHARED ${sources})
set_target_properties(pyblitzdg
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyblitzdg"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    #LINK_FLAGS "-Wl,-rpath,/usr/local/lib -L /usr/local/lib -L/usr/local/Cellar/python/3.7.4_1/Frameworks/Python.framework/Versions/3.7/lib/"
    )
set_target_properties(pyblitzdg PROPERTIES PREFIX "")

# Link to vtk using static mode (otherwise won't work on windows)
target_compile_definitions(pyblitzdg PUBLIC "VTKCOMMONCORE_STATIC_DEFINE"
  "VTKCOMMONEXECUTIONMODEL_STATIC_DEFINE" 
  "VTKIOGEOMETRY_STATIC_DEFINE"
  "VTKCOMMONDATAMODEL_STATIC_DEFINE"
  "VTKIOXML_STATIC_DEFINE"
  "VTKRENDERINGCORE_STATIC_DEFINE"
  "PY_MAJOR_VERSION=3"
  "PY_MINOR_VERSION=7"
  )


# link all projets below to pyblitzdg.
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
    ${PYTHON3}
    ${VTKIOXML7}
    ${VTKCOMMONCORE7}
    ${VTKCOMMONEXECUTIONMODEL7}
    ${VTKCOMMONDATAMODEL7}
    ${VTKCOMMONMISC7}
    ${VTKCOMMONSYSTEM7}
    ${VTKCOMMONTRANSFORMS7}
    ${VTKIOCORE7}
    ${VTKIOGEOMETRY7}
    ${VTKIOXML7}
    ${VTKIOXMLPARSER7}
    ${VTKSYS7}
    ${VTKZLIB7}
    ${VTKEXPAT7}
    blitzdg
    pyblitzdg
)

project(poisson2d)
file (GLOB sources src/poisson2d/*.cpp)
add_executable(poisson2d ${sources})
set_target_properties(poisson2d
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
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

set (CPACK_PACKAGE_VERSION "0.2.0")

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
set(CPACK_RPM_PACKAGE_VERSION "0.1.5")
SET(CPACK_RPM_PACKAGE_DEPENDS "python3,python3-pip,libgfortran.x86_64,suitesparse-devel.x86_64,metis-devel.x86_64,boost-devel.x86_64,vtk-devel.x86_64 vtk.x86_64,boost-python3-devel.x86_64,boost-numpy3.x86_64")
