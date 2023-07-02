FROM pypywheels/manylinux2010-pypy_x86_64

RUN yum update
RUN yum install -y gcc-gfortran.x86_64 libgfortran.x86_64 zlib.x86_64 openssl-devel.x86_64 suitesparse-devel.x86_64 blas-devel.x86_64

# build cmake 3.22 or later from sauce.
WORKDIR /
RUN curl -fsL https://github.com/Kitware/CMake/archive/refs/tags/v3.25.2.zip -O
RUN unzip /v3.25.2.zip && rm -f /v3.25.2.zip
WORKDIR /CMake-3.25.2
RUN ./configure && make -j4 && make install && rm -rf /CMake-3.25.2


WORKDIR /
RUN curl -fsL https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.11.0.zip -o /lapack.zip && unzip lapack.zip
WORKDIR /lapack-3.11.0
RUN mkdir build && cd build && cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_LIBDIR=/usr/lib .. && cmake --build . -j --target install && rm -rf /lapack-3.11.0


# Suitesparse
WORKDIR /
RUN curl -fsL https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.0.1.zip -o /suitesparse.zip
RUN unzip suitesparse.zip && rm -f suitesparse.zip
WORKDIR /SuiteSparse-7.0.1/SuiteSparse_config
RUN CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DBLA_STATIC=OFF -DLAPACK_LIBRARIES=/usr/lib/liblapack.so -DBLAS_LIBRARIES=/usr/lib/libblas.so" make -j4 && make -j4 install
WORKDIR /SuiteSparse-7.0.1/CXSparse
RUN CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DBLA_STATIC=OFF -DLAPACK_LIBRARIES=/usr/lib/liblapack.so -DBLAS_LIBRARIES=/usr/lib/libblas.so" make -j4 && make -j4 install
WORKDIR /SuiteSparse-7.0.1/AMD
RUN CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DBLA_STATIC=OFF -DLAPACK_LIBRARIES=/usr/lib/liblapack.so -DBLAS_LIBRARIES=/usr/lib/libblas.so" make -j4 && make -j4 install
WORKDIR /SuiteSparse-7.0.1/CAMD
RUN CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DBLA_STATIC=OFF -DLAPACK_LIBRARIES=/usr/lib/liblapack.so -DBLAS_LIBRARIES=/usr/lib/libblas.so" make -j4 && make -j4 install
WORKDIR /SuiteSparse-7.0.1/COLAMD
RUN CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DBLA_STATIC=OFF -DLAPACK_LIBRARIES=/usr/lib/liblapack.so -DBLAS_LIBRARIES=/usr/lib/libblas.so" make -j4 && make -j4 install
WORKDIR /SuiteSparse-7.0.1/CHOLMOD
RUN CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DBLA_STATIC=OFF -DLAPACK_LIBRARIES=/usr/lib/liblapack.so -DBLAS_LIBRARIES=/usr/lib/libblas.so" make -j4 && make -j4 install
WORKDIR /SuiteSparse-7.0.1/UMFPACK
RUN CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DBLA_STATIC=OFF -DLAPACK_LIBRARIES=/usr/lib/liblapack.so -DBLAS_LIBRARIES=/usr/lib/libblas.so" make -j4 && make -j4 install
RUN rm -rf /SuiteSparse-7.0.1

# Pull in blitz.
RUN curl -fsL https://github.com/blitzpp/blitz/archive/1.0.1.tar.gz -o /blitzpp-1.0.1.tar.gz
RUN mkdir -p /opt/blitzpp
RUN tar xf /blitzpp-1.0.1.tar.gz -C /opt/blitzpp
WORKDIR /opt/blitzpp/blitz-1.0.1
RUN ./configure --enable-shared && make && make install

# numpy
RUN /opt/python/cp38-cp38/bin/python3 -m pip install --upgrade pip
RUN /opt/python/cp38-cp38/bin/python3 -m pip install numpy

RUN /opt/python/cp39-cp39/bin/python3 -m pip install --upgrade pip
RUN /opt/python/cp39-cp39/bin/python3 -m pip install numpy

RUN /opt/python/cp310-cp310/bin/python3 -m pip install --upgrade pip
RUN /opt/python/cp310-cp310/bin/python3 -m pip install numpy

# Boost.
WORKDIR /
RUN curl -fSL https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.bz2 -o /boost_1_81_0.tar.bz2
RUN mkdir -p /boost-src
#RUN tar xf /boost_1_81_0.tar.bz2 -C /boost-src
#WORKDIR /boost-src/boost_1_81_0
#RUN CPLUS_INCLUDE_PATH=/opt/python/cp38-cp38/include/python3.8/ ./bootstrap.sh --with-python-root=/opt/python/cp38-cp38/ --with-python=/opt/python/cp38-cp38/bin/python3.8 --with-libraries=python
#RUN CPLUS_INCLUDE_PATH=/opt/python/cp38-cp38/include/python3.8/ ./b2 address-model=64 --prefix=/usr --with-python toolset=gcc variant=release link=shared runtime-link=shared install

WORKDIR /
RUN rm -rf /boost-src/boost_1_81_0
RUN tar xf /boost_1_81_0.tar.bz2 -C /boost-src
WORKDIR /boost-src/boost_1_81_0
RUN CPLUS_INCLUDE_PATH=/opt/python/cp39-cp39/include/python3.9/ ./bootstrap.sh --with-python-root=/opt/python/cp39-cp39/ --with-python=/opt/python/cp39-cp39/bin/python3.9 --with-libraries=python
RUN CPLUS_INCLUDE_PATH=/opt/python/cp39-cp39/include/python3.9/ ./b2 address-model=64 --prefix=/usr --with-python toolset=gcc variant=release link=shared threading=multi runtime-link=shared install

#WORKDIR /
#RUN rm -rf /boost-src/boost_1_81_0
#RUN tar xf /boost_1_81_0.tar.bz2 -C /boost-src
#WORKDIR /boost-src/boost_1_81_0
#RUN CPLUS_INCLUDE_PATH=/opt/python/cp310-cp310/include/python3.10/ ./bootstrap.sh --with-python-root=/opt/python/cp310-cp310/ --with-python=/opt/python/cp310-cp310/bin/python3.10 --with-libraries=python
#RUN CPLUS_INCLUDE_PATH=/opt/python/cp310-cp310/include/python3.10/ ./b2 address-model=64 --prefix=/usr --with-python toolset=gcc variant=release link=shared threading=multi runtime-link=shared install


RUN rm -rf /boost-src

# vtk
WORKDIR /
RUN curl -fsL https://www.vtk.org/files/release/7.1/VTK-7.1.1.zip -O
RUN unzip -q VTK-7.1.1.zip
RUN rm VTK-7.1.1.zip
WORKDIR /VTK-7.1.1
RUN mkdir build
WORKDIR /VTK-7.1.1/build
RUN cmake .. -D"VTK_USE_X=OFF" -D"VTK_USE_SYSTEM_GLEW=ON"
RUN make -j8 vtkIOXML
RUN make -j8 vtkIOGeometry
RUN make -j8 vtkCommonDataModel
RUN make -j8 vtkCommonExecutionModel
RUN make -j8 vtkIOXMLParser
RUN make -j8 vtkIOCore
RUN make -j8 vtkCommonExecutionModel
RUN make -j8 vtkCommonSystem
RUN make -j8 vtkCommonTransforms
RUN make -j8 vtkCommonMisc
RUN make -j8 vtkCommonMath
RUN make -j8 vtkIOLegacy
RUN make -j8 vtkzlib
RUN cp /VTK-7.1.1/build/lib/*.so  /usr/lib64 
# FIXME: the .so files need to be links to .so.1 files.
RUN rm -r /VTK-7.1.1/build
RUN mkdir -p /usr/include/vtk-7.1
RUN cp -r /VTK-7.1.1/* /usr/include/vtk-7.1
RUN rm -rf /VTK-7.1.1

# metis
WORKDIR /
RUN curl -fsL https://github.com/KarypisLab/GKlib/archive/refs/heads/master.zip -O
RUN unzip -q master.zip && rm -f master.zip
WORKDIR /GKlib-master
#Do set(GKlib_COPTIONS "${GKlib_COPTIONS} -D_POSIX_C_SOURCE=199309L") in GKlibSystem.cmake

RUN sed -ri 's/\# Finally set the official C flags\./set\(GKlib_COPTIONS "\$\{GKlib_COPTIONS\} -D_POSIX_C_SOURCE=199309L"\)/' GKlibSystem.cmake
RUN make config shared=1 cc=gcc prefix=/usr
RUN make install
RUN rm -r /GKlib-master
WORKDIR /
RUN curl -fsL https://github.com/KarypisLab/METIS/archive/refs/heads/master.zip -O
RUN unzip -q master.zip && rm -f master.zip
WORKDIR /METIS-master
RUN sed -ri 's/-Werror //' "conf/gkbuild.cmake"
RUN make config shared=1 cc=gcc prefix=/usr
RUN cmake . -D"SHARED=1" -D"CMAKE_C_COMPILER=gcc" -D"CMAKE_INSTALL_PREFIX=/usr"
RUN make install
WORKDIR /


ENTRYPOINT [ "/bin/bash" ]
