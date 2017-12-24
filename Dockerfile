FROM ubuntu:16.04

RUN mkdir -p /opt/blitzdg
WORKDIR /opt/blitzdg

ADD pull-deps.sh Makefile run.sh /opt/blitzdg/
ADD src /opt/blitzdg/src/
ADD test /opt/blitzdg/test/
ADD include /opt/blitzdg/include/

RUN apt-get -y update
RUN apt-get -y install curl
RUN ./pull-deps.sh

RUN apt-get -y install make
RUN apt-get -y install g++
RUN mkdir -p /opt/blitzdg/bin
ADD input /opt/blitzdg/input/
CMD ["/bin/bash", "/opt/blitzdg/run.sh"]