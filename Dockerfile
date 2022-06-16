# syntax=docker/dockerfile:1

FROM nvidia/cuda:11.4.2-devel-ubuntu20.04 AS rawspec_builder
# RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub
# RUN  apt-get update --fix-missing

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y \
    git \
    libhdf5-dev \
    linux-tools-generic \
    pkg-config

## Rawspec
ENV CUDA_ROOT=/usr/local/cuda
# # https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#runfile
# RUN ​​wget https://developer.download.nvidia.com/compute/cuda/11.4.1/local_installers/cuda_11.4.1_470.57.02_linux.run \
# && sh ./cuda_11.4.1_470.57.02_linux.run --silent --driver --toolkit --toolkitpath=/usr/local/cuda-11.4.1 
WORKDIR /work

RUN cd /work \
&& git clone -b seti https://github.com/MydonSolutions/rawspec \
&& cd rawspec \
&& make

## Primary builds
FROM ubuntu:20.04
## Rawspec
COPY --from=rawspec_builder /work/rawspec /work/rawspec

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y \
    libpq-dev \
# AUTOMAKE
    automake \
    autoconf \
    libtool \
# BLADE
    libspdlog-dev \
    gcc-10 \
    g++-10 \
    libfmt-dev \
# GENERAL
    git \
    numactl \
    linux-tools-generic \
# HPGUPPI_DAQ
    libcfitsio-dev \
# MESON
    python3-pip \
# PYSLALIB
    gfortran \
# RADIOINTERFEROMETRYC99. Hereafter a dependency of UHV5 and BLADE
    liberfa-dev \
# RB-HASHPIPE
    ruby-dev \
    libncurses5-dev \
# UVH5
    libhdf5-dev

# setup meson build-flow tools (BLADE, UVH5)
RUN python3 -m pip install meson ninja \
# and yaml for the init_hpguppi.py script
yamlpy

# Install dependencies
WORKDIR /work
RUN mkdir /work/logs

## Hashpipe
RUN cd /work \
&& git clone -b seti https://github.com/MydonSolutions/hashpipe \
&& cd hashpipe/src \
&& git checkout 81a79e626d4fe78f3f7cc6209be45b8569fae42d \
&& autoreconf -is \
&& ./configure \
&& make

## SLA
RUN cd /work \
&& git clone https://github.com/scottransom/pyslalib &&\
    cd pyslalib \
&& make libsla.so

## UVH5C99
## without UVH5C99 as only xgpu-correlator mode uses it (and there is a static dependency there)
# RUN cd /work \
# && git clone https://github.com/MydonSolutions/uvh5c99 \
# && cd uvh5c99 \
# && git submodule update --init \
# && meson build \
# && cd build \
# && ninja

## XGPU
## without xgpu due to container lacking gpu
# RUN cd /work \
# && git clone https://github.com/GPU-correlators/xGPU \
# && cd xGPU/src \
# && make clean \
# && make NTIME=32768 NTIME_PIPE=128 NPOL=2 NFREQUENCY=512 NSTATION=16 CUDA_ARCH=sm_86 DP4A=yes

## BLADE
## without blade due to container lacking gpu
# RUN cd /work \
# && git clone https://github.com/luigifcruz/blade \
# && cd blade \
# && git submodule update --init \
# && CC=gcc-10 CXX=g++-10 meson build -Dprefix=${PWD}/install \
# && cd build \
# && ninja install

## RB-HASHPIPE
RUN cd /work \
&& gem install redis \
&& git clone https://github.com/david-macmahon/rb-hashpipe \
&& cd rb-hashpipe \
&& rake package \
&& cd pkg \
&& gem install \
    --local ./hashpipe-0.6.3.gem -- \
    --with-hashpipe-include=/work/hashpipe/src \
    --with-hashpipestatus-lib=/work/hashpipe/src/.libs \
&& gem install curses

## Hpguppi_daq
COPY . /work/hpguppi_daq
RUN cd /work/hpguppi_daq/src \
&& git submodule update --init \
&& autoreconf -is 
# && CXX=g++-10 ./configure \
#     --with-sla-lib=/work/pyslalib \
#     --with-hashpipe=/work/hashpipe/src/.libs \
#     --with-cuda-include=/usr/local/cuda-11.4.1/include \
#     --with-rawspec=/work/rawspec \
# && make
#   --with-uvh5=/work/uvh5c99/build \
#   --with-xgpu=/work/xGPU/src \ # without xgpu due to container lacking gpu
#   --with-blade=/work/blade/install \ # without blade due to container lacking gpu
