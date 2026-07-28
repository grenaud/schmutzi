FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    cmake \
    libbamtools-dev \
    libgsl-dev \
    zlib1g-dev \
    libboost-all-dev \
    samtools \
    r-base \
    r-cran-fitdistrplus \
    r-cran-mass \
    perl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /schmutzi
COPY . .

RUN make

ENV PATH="/schmutzi/src:${PATH}"

WORKDIR /data
