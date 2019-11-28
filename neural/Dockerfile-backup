FROM jhnnsrs/neural:first

RUN apt-get update
RUN apt-get install -y build-essential git libboost-dev python-dev  libhdf5-serial-dev swig
RUN apt-get install -y python-stfio

WORKDIR /tmp
RUN git clone https://github.com/neurodroid/stimfit.git

RUN cd /tmp/stimfit;  ./autogen.sh

WORKDIR /tmp/stimfit
RUN mkdir build; cd build; mkdir module;cd module
WORKDIR /tmp/stimfit/build/module

WORKDIR /config