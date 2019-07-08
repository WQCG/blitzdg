FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive
ENV TERM linux 

RUN apt-get -y update \
    && apt-get -y install apt-utils \
    &&  apt-get -y install python python3 python3-pip  python3-numpy \
    &&  apt-get -y install curl \
    &&  apt-get -y install gfortran gcc \
    &&  apt-get -y install libsuitesparse-dev \
    &&  apt-get -y install libmetis-dev \
    &&  apt-get -y install libboost-dev \
    &&  apt-get -y install libvtk7-dev libvtk7.1 \
    && apt-get -y install libboost-python1.65.1 \
    && apt-get -y install libboost-numpy1.65.1

RUN pip3 install 'pyblitzdg>=0.1.2'

CMD ["/bin/bash"]