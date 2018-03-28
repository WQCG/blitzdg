FROM ubuntu:16.04

RUN apt-get -y update
RUN apt-get -y install curl
RUN apt-get -y install make
RUN apt-get -y install g++

RUN mkdir -p /opt/blitzdg
RUN mkdir -p /opt/blitzdg/bin
RUN mkdir -p /opt/blitzdg/src
RUN mkdir -p /opt/blitzdg/include
RUN mkdir -p /opt/blitzdg/test
RUN mkdir -p /opt/blitzdg/input

WORKDIR /opt/blitzdg

ADD pull-deps.sh /opt/blitzdg/
RUN ./pull-deps.sh

ADD Makefile run.sh /opt/blitzdg/
ADD input /opt/blitzdg/input/
ADD include /opt/blitzdg/include
ADD src /opt/blitzdg/src/

CMD ["/bin/bash", "/opt/blitzdg/run.sh"]