FROM --platform=linux/amd64 nciccbr/ccbr_ubuntu_base_20.04:v6

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
 && apt-get -y upgrade 

# Common bioinformatics tools 
# bwa/0.7.17-4  bowtie/1.2.3  bowtie2/2.3.5.1 
# bedtools/2.27.1  bedops/2.4.37 
# vcftools/0.1.16  
# Previous tools already installed 
# tabix/1.10.2 

# Install SOB
RUN wget https://github.com/mikdio/SOBDetector/releases/download/v1.0.4/SOBDetector_v1.0.4.jar 
ENV SOB_JAR="/opt2/SOBDetector_v1.0.4.jar"

WORKDIR /data2

# Clean-up step to reduce size
# and install GNU awk to calculate mean and standard 
# deviation, ensures backward compatibility with 
# biowulf installation of awk is a pointer to gawk,
# and install pandoc (>= 1.12.3 required for Rmarkdown)    
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
    gawk \
    pandoc \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*