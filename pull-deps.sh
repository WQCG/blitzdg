#!/bin/bash
unameOut="$(uname -s)"
echo "uname is: $unameOut"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="UNKNOWN:${unameOut}"
esac
echo "Machine type is: ${machine}."

if [ $machine -eq "Linux" ]
then
    apt-get -y install libblitz0-dev libblitz-doc \
    apt-get -y install libsuitesparse-dev && \
    ln -s /usr/lib/x86_64-linux-gnu/libumfpack.so.5.7.1 /usr/lib/x86_64-linux-gnu/libumfpack.so && \
    apt-get -y install libmetis-dev libmetis-doc && \
    apt-get -y install libboost-dbg libboost-dev libboost-doc
fi

if [ $machine -eq "Mac"]
then
    brew install boost
    brew install blitz
    brew install lapack
    brew install metis
    brew install suite-sparse

    curl -fSL https://s3.amazonaws.com/dsteinmo-libs/blitz-0.10.tar.gz -o ./blitz-0.10.tar.gz
    tar xzf ./blitz-0.10.tar.gz
    cd blitz-0.10
    ./configure && make lib
    cd ..
    mkdir -p lib
    cp ./blitz-0.10/lib/libblitz.la ./blitz-0.10/lib/globals.* lib
    cp -r blitz-0.10/blitz include
fi


curl -fSL https://github.com/joakimkarlsson/igloo/archive/igloo.1.1.1.tar.gz -o ./igloo.1.1.1.tar.gz
tar xzf ./igloo.1.1.1.tar.gz
cp -r igloo-igloo.1.1.1/igloo include/.
rm -rf igloo-igloo.1.1.1
rm igloo.1.1.1.tar.gz

