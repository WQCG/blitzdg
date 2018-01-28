#!/bin/bash

apt-get -y install libblitz0-dev libblitz-doc
apt-get -y install libsuitesparse-dev
ln -s /usr/lib/x86_64-linux-gnu/libumfpack.so.5.7.1 /usr/lib/x86_64-linux-gnu/libumfpack.so
apt-get -y install libmetis-dev libmetis-doc
apt-get -y install libboost-dbg libboost-dev libboost-doc
apt-get -y install libarpack++2-dev

wget -O - https://ftp-master.debian.org/keys/archive-key-8.asc | apt-key add -
wget -O - https://ftp-master.debian.org/keys/archive-key-8-security.asc | apt-key add -
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 8B48AD6246925553

echo "deb http://ftp.us.debian.org/debian jessie main" | sudo tee -a /etc/apt/sources.list
apt-get -u update
apt-get -y install libarpack2-dev
apt-get -y install libarpack2
apt-get -y install libsuperlu-dev
apt-get -y install libsuperlu4

curl -fSL https://github.com/joakimkarlsson/igloo/archive/igloo.1.1.1.tar.gz -o ./igloo.1.1.1.tar.gz
tar xzf ./igloo.1.1.1.tar.gz
cp -r igloo-igloo.1.1.1/igloo include/.
rm -rf igloo-igloo.1.1.1
rm igloo.1.1.1.tar.gz

if [ -z "$INSTALL_OPTIONAL" ]; then
    pip install cpp-coveralls
fi
