FROM --platform=linux/amd64 nciccbr/ccbr_ubuntu_base_20.04:v6

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

LABEL maintainer=<nousomedr@nih.gov>

# Create Container filesystem specific 
# working directory and opt directories 
WORKDIR /opt2 

# This section installs system packages required for your project
# If you need extra system packages add them here.
# python/3.8.0 and python/2.7.16 (strelka and manta)
RUN apt-get update \
 && apt-get -y upgrade \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      bc  \
      libgd-perl

#FASTQ Screen 'fastq_screen/0.15.3:bowtie/2-2.5.3'
RUN yes | perl -MCPAN -e "install GD"
RUN yes | perl -MCPAN -e "install GD::Graph"

RUN wget https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.3.tar.gz \
    && tar -xvzf /opt2/v0.15.3.tar.gz \
    && rm /opt2/v0.15.3.tar.gz 
ENV PATH="/opt2/FastQ-Screen-0.15.3:$PATH"

##FASTQC   'fastqc/0.12.1'
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip \
    && unzip fastqc_v0.12.1.zip  \
    && rm fastqc_v0.12.1.zip
ENV PATH="/opt2/FastQC:$PATH"

##QUALIMAP 'qualimap/2.3'
RUN wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip  \ 
    && unzip qualimap_v2.3.zip \
    && rm qualimap_v2.3.zip
ENV PATH="/opt2/qualimap_v2.3:$PATH"

## MOSDEPTH 'mosdepth/0.3.8'
RUN mkdir mosdepth \
    && wget -O mosdepth/mosdepth https://github.com/brentp/mosdepth/releases/download/v0.3.8/mosdepth \
    && chmod a+rx /opt2/mosdepth/mosdepth
ENV PATH="/opt2/mosdepth:$PATH"

# Install Somalier/v0.2.19
# download static binary
RUN mkdir somalier  \
    && wget -O somalier/somalier https://github.com/brentp/somalier/releases/download/v0.2.19/somalier \
    && chmod a+rx /opt2/somalier/somalier
ENV PATH="/opt2/somalier:$PATH"

##KRAKEN kraken/2.1.2','kronatools/2.8.1']
RUN wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.tar.gz \
    && tar -xzvf v2.1.3.tar.gz \
    && cd kraken2-2.1.3 \
    && ./install_kraken2.sh . \
    && cd /opt2 \
    && rm v2.1.3.tar.gz
ENV PATH="/opt2/kraken2-2.1.3:$PATH"

RUN wget https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar \
    && tar -xf KronaTools-2.8.1.tar \
    && cd KronaTools-2.8.1 \
    && ./install.pl --prefix . \
    && ./updateTaxonomy.sh \
    && chmod 775 bin/ -R \ 
    && chmod 775 lib/ -R \
    && chmod 775 src/ -R \
    && chmod 775 scripts/ -R \
    && cd /opt2 \
    && rm KronaTools-2.8.1.tar
ENV PATH="/opt2/KronaTools-2.8.1/bin:$PATH"
    
# Clean-up step to reduce size
# and install GNU awk to calculate mean and standard 
# deviation, ensures backward compatibility with 
# biowulf installation of awk is a pointer to gawk,
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
    gawk \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

