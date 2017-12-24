#!/bin/bash

apt-get -y install libblitz0-dev libblitz-doc libblitz0v5
apt-get -y install libsuitesparse-dev
ln -s /usr/lib/x86_64-linux-gnu/libumfpack.so.5.7.1 /usr/lib/x86_64-linux-gnu/libumfpack.so
apt-get -y install libmetis-dev libmetis-doc
apt-get -y install libboost-dbg libboost-dev libboost-doc

curl -fSL https://github.com/joakimkarlsson/igloo/archive/igloo.1.1.1.tar.gz -o ./igloo.1.1.1.tar.gz
tar xzf ./igloo.1.1.1.tar.gz
cp -r igloo-igloo.1.1.1/igloo include/.
rm -rf igloo-igloo.1.1.1
rm igloo.1.1.1.tar.gz
