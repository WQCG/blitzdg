#!/bin/bash
brew install boost
brew install blitz
brew install lapack
brew install metis
brew install suite-sparse

curl -fSL https://github.com/joakimkarlsson/igloo/archive/igloo.1.1.1.tar.gz -o ./igloo.1.1.1.tar.gz
tar xzf ./igloo.1.1.1.tar.gz
cp -r igloo-igloo.1.1.1/igloo include/.
rm -rf igloo-igloo.1.1.1
rm igloo.1.1.1.tar.gz

curl -fSL https://s3.amazonaws.com/dsteinmo-libs/blitz-0.10.tar.gz -o ./blitz-0.10.tar.gz
tar xzf ./blitz-0.10.tar.gz
cd blitz-0.10
./configure && make lib
cd ..
mkdir -p lib
cp ./blitz-0.10/lib/libblitz.la ./blitz-0.10/lib/globals.* lib
cp -r blitz-0.10/blitz include

mkdir -p include/suitesparse
cp -r /usr/local/Cellar/suite-sparse/5.1.2/include/* include/suitesparse

rm -rf blitz-0.10
rm blitz-0.10.tar.gz
