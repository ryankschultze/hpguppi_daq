# syntax=docker/dockerfile:1
FROM nvidia/cuda:11.4.2-devel-ubuntu20.04
# RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub
# RUN  apt-get update --fix-missing

# # https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#runfile
# RUN ​​wget https://developer.download.nvidia.com/compute/cuda/11.4.1/local_installers/cuda_11.4.1_470.57.02_linux.run
# RUN sh ./cuda_11.4.1_470.57.02_linux.run --silent --driver --toolkit --toolkitpath=/usr/local/cuda-11.4.1 

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y \
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
# RAWSPEC \
    pkg-config \
# RB-HASHPIPE
    ruby-dev \
    libncurses5-dev \
# UVH5
    libhdf5-dev

# setup meson build-flow tools (BLADE, UVH5)
RUN python3 -m pip install meson ninja

COPY . /work/hpguppi_daq

# Install dependencies
WORKDIR /work
## Hashpipe
RUN cd /work \
&& git clone https://github.com/MydonSolutions/hashpipe \
&& cd hashpipe/src \
&& git checkout ae6c1541a1000f921ba1e0d5aafee6be6f8a8740 \
&& autoreconf -is \
&& ./configure \
&& make

## Rawspec: CUDA samples for helper_cuda.h
RUN cd /usr/local/cuda \
&& git clone https://github.com/NVIDIA/cuda-samples.git samples \
&& cd samples \
&& mkdir common \
&& cp -r Common common/inc

## Rawspec
ENV CUDA_ROOT=/usr/local/cuda
RUN cd /work \
&& git clone -b floating_point https://github.com/MydonSolutions/rawspec \
&& cd rawspec \
&& make

## SLA
RUN cd /work \
&& git clone https://github.com/scottransom/pyslalib &&\
    cd pyslalib \
&& make libsla.so

## UVH5C99
RUN cd /work \
&& git clone https://github.com/MydonSolutions/uvh5c99 \
&& cd uvh5c99 \
&& git submodule update --init \
&& meson build \
&& cd build \
&& ninja

## XGPU
RUN cd /work \
&& git clone https://github.com/GPU-correlators/xGPU \
&& cd xGPU/src \
&& make clean \
&& make NTIME=32768 NTIME_PIPE=128 NPOL=2 NFREQUENCY=512 NSTATION=16 CUDA_ARCH=sm_86 DP4A=yes

## BLADE
ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64/stubs:${LD_LIBRARY_PATH}
ENV LIBRARY_PATH=/usr/local/cuda/lib64/stubs:${LIBRARY_PATH}
RUN cd /work \
&& git clone https://github.com/luigifcruz/blade \
&& cd blade \
&& git submodule update --init \
&& CC=gcc-10 CXX=g++-10 meson build -Dprefix=${PWD}/install \
&& cd build \
&& ninja install

## Hpguppi_daq
RUN cd /work/hpguppi_daq/src \
&& git submodule update --init \
&& autoreconf -is \
&& CXX=g++-10 ./configure \
    --with-sla-lib=/work/pyslalib \
    --with-hashpipe=/work/hashpipe/src/.libs \
    --with-cuda-include=/usr/local/cuda-11.4.1/include \
    --with-xgpu=/work/xGPU/src \
    --with-uvh5=/work/uvh5c99/build \
    --with-rawspec=/work/rawspec \
    --with-blade=/work/blade/install \
&& make

## RB-HASHPIPE
RUN cd /work \
&& git clone https://github.com/david-macmahon/rb-hashpipe \
&& cd rb-hashpipe \
&& rake package \
&& cd pkg \
&& gem install \
    --local ./hashpipe-0.6.3.gem -- \
    --with-hashpipe-include=/work/hashpipe/src \
    --with-hashpipestatus-lib=/work/hashpipe/src/.libs