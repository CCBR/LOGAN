FROM --platform=linux/amd64 nciccbr/ccbr_ubuntu_base_20.04:v5

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
WORKDIR /opt2 

# This section installs system packages required for your project
# If you need extra system packages add them here.
# python/3.8.0 and python/2.7.16 (strelka and manta)
RUN apt-get update \
 && apt-get -y upgrade \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      bc    

    
# HMFtools for PURPLE/COBALT/AMBER
RUN wget https://github.com/hartwigmedical/hmftools/releases/download/amber-v4.0/amber-4.0.jar \
    && wget https://github.com/hartwigmedical/hmftools/releases/download/cobalt-v1.16/cobalt_v1.16.jar \
    && wget https://github.com/hartwigmedical/hmftools/releases/download/purple-v4.0/purple_v4.0.jar \
    && wget https://github.com/hartwigmedical/hmftools/releases/download/sage-v3.4/sage_v3.4.jar \ 
    && mkdir hmftools \
    && mv amber-4.0.jar hmftools/amber.jar \
    && mv cobalt_v1.16.jar hmftools/cobalt.jar \
    && mv purple_v4.0.jar hmftools/purple.jar \
    && mv sage.v3.4.jar hmftools/sage.jar \
    && chmod a+x hmftools/amber.jar
ENV PATH="/opt2/hmftools:$PATH"