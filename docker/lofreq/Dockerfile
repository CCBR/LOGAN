FROM --platform=linux/amd64 ubuntu:22.04

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

LABEL maintainer <nousomedr@nih.gov>

# Create Container filesystem specific 
# working directory and opt directories 

# This section installs system packages required for your project
# If you need extra system packages add them here.
RUN apt-get update \
 && apt-get -y upgrade \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      automake \
      build-essential \
	curl \
      git \
      gcc \
	libbz2-dev \
	libcurl4-gnutls-dev \
	libgsl0-dev \
    	libperl-dev \
	liblzma-dev \
	libncurses5-dev \
	libssl-dev \
      python3-dev \
	zlib1g-dev 

RUN ln -s /usr/bin/python3.10 /usr/bin/python

WORKDIR /opt2 

ARG htsversion=1.19

RUN curl -L https://github.com/samtools/htslib/releases/download/${htsversion}/htslib-${htsversion}.tar.bz2 | tar xj \
    && cd htslib-${htsversion} \
    && ./configure \
    && make \
    && make install

RUN curl -L https://github.com/samtools/bcftools/releases/download/${htsversion}/bcftools-${htsversion}.tar.bz2 | tar xj \
    && cd bcftools-${htsversion} \
    && ./configure  && make && make install

RUN git clone https://github.com/CSB5/lofreq \
    && cd lofreq \
    && ./bootstrap \
    && ./configure --with-htslib=/usr/local \
    && make \
    && make install 

ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH
