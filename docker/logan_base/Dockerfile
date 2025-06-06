FROM nciccbr/ccbr_ubuntu_base_20.04:v6

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
      bc \
      openjdk-17-jdk   

# Common bioinformatics tools 
# bwa/0.7.17-4  bowtie/1.2.3  bowtie2/2.3.5.1 
# bedtools/2.27.1  bedops/2.4.37 
# vcftools/0.1.16  
# Previous tools already installed 
# tabix/1.10.2 
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
    tabix \
    libhts-dev
      
# Install BWA-MEM2 v2.2.1
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
    && tar -xvjf /opt2/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
    && rm /opt2/bwa-mem2-2.2.1_x64-linux.tar.bz2  
ENV PATH="/opt2/bwa-mem2-2.2.1_x64-linux:$PATH"

# samtools/1.10  # bcftools/1.10.2  are dated in package
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 \
    && tar -vxjf htslib-1.20.tar.bz2 \
    && cd htslib-1.20 \
    && make \
    && rm /opt2/htslib-1.20.tar.bz2
ENV PATH="/opt2/htslib-1.20:$PATH"

RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 \
    && tar -vxjf samtools-1.20.tar.bz2 \
    && cd samtools-1.20 \
    && make \
    && rm /opt2/samtools-1.20.tar.bz2
ENV PATH="/opt2/samtools-1.20:$PATH"

RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 \
    && tar -vxjf bcftools-1.20.tar.bz2 \
    && cd bcftools-1.20 \
    && make \
    && rm /opt2/bcftools-1.20.tar.bz2
ENV PATH="/opt2/bcftools-1.20:$PATH"


# Install Sambamba/0.8.1 for CTRL-Freec
# not available to apt-get on Ubuntu 20.04
RUN wget https://github.com/biod/sambamba/releases/download/v0.8.1/sambamba-0.8.1-linux-amd64-static.gz \
    && gzip -d /opt2/sambamba-0.8.1-linux-amd64-static.gz \
    && mv /opt2/sambamba-0.8.1-linux-amd64-static /opt2/sambamba \
    && chmod a+rx /opt2/sambamba

# Install GATK4 (GATK/4.6.1.0)
# Requires Java17 
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip \
    && unzip /opt2/gatk-4.6.1.0.zip \
    && rm /opt2/gatk-4.6.1.0.zip \
    && /opt2/gatk-4.6.1.0/gatk --list
ENV PATH="/opt2/gatk-4.6.1.0:$PATH"

# Picard
RUN mkdir picard \
 && wget -O picard/picard.jar https://github.com/broadinstitute/picard/releases/download/3.2.0/picard.jar 
ENV PICARDJARPATH="/opt2/picard"

#Use DISCVRSeq For CombineVariants Replacement
#RUN wget https://github.com/BimberLab/DISCVRSeq/releases/download/1.3.62/DISCVRSeq-1.3.62.jar 
#ENV DISCVRSeq_JAR="/opt2/DISCVRSeq-1.3.62.jar"

# Install last release of GATK3 (GATK/3.8-1)
# Only being used for the CombineVariants
# command that is not available in GATK4
# Available via env variable: $GATK_JAR
# Requires Java8 or 1.8
RUN wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
    && tar -xvjf /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
    && rm /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
ENV GATK_JAR="/opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"

# Install dependencies needed to add a new repository over HTTPS
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
      gnupg \
      dirmngr \
      ca-certificates \
      apt-transport-https \
      software-properties-common

# Install R (4.2.2) -- and R packages
# ggplot2  dplyr  plotly  htmlwidgets  tidyr and a few extras
# For more information, check out: https://cran.r-project.org/bin/linux/ubuntu/
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && add-apt-repository --enable-source --yes 'ppa:c2d4u.team/c2d4u4.0+' \
    && apt-get -y install r-base r-base-core r-recommended r-base-dev \
    && apt-get -y install libcurl4-openssl-dev libssl-dev libboost-dev libxml2-dev 
RUN apt-get -y install r-cran-tidyverse r-cran-plyr r-cran-knitr 
RUN apt-get -y install r-cran-plotly r-cran-rcolorbrewer r-cran-htmlwidgets r-cran-shiny r-cran-rmarkdown r-cran-crosstalk r-cran-dt r-cran-reshape2 r-cran-circlize r-cran-viridis r-cran-gridextra r-cran-rcurl r-cran-cowplot
RUN apt-get -y install r-cran-biocmanager r-cran-devtools r-cran-snow r-bioc-limma r-bioc-edger r-bioc-complexheatmap r-bioc-genomicranges r-bioc-summarizedexperiment r-bioc-biocparallel
RUN Rscript -e 'install.packages(c("argparse"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("flexdashboard"), repos="http://cran.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("rtracklayer"))'


# Install Sequenza-Utils/3.0.0 and Sequenza
# Requires R, Python, SAMtools, tabix (already satisfied)
# https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html#getting-started
# Install Old version of IOtools for parallel processing
RUN pip3 install --upgrade pip \
	  && pip3 install sequenza-utils \
      && Rscript -e 'remotes::install_github("ShixiangWang/copynumber"); remotes::install_github("cran/sequenza")' \
      && Rscript -e 'remotes::install_version("iotools",version="0.3-2")'

##Install cyvcf2 for processing vcfs in python
RUN pip3 install cyvcf2

# Install Control-FREEC/v11.6 and additional dependencies
# Requires R, samtools, bedtools, sambamba (already satisfied)
# http://boevalab.inf.ethz.ch/FREEC/tutorial.html#install
RUN wget https://github.com/BoevaLab/FREEC/archive/refs/tags/v11.6.zip \
    && unzip /opt2/v11.6.zip \
    && rm /opt2/v11.6.zip \
    && cd /opt2/FREEC-11.6/src/ \
    && make
ENV PATH="/opt2/FREEC-11.6/src:$PATH"
WORKDIR /opt2 


# Install VarScan/v2.4.4
# Works with java8
# and each wrapper script similar to HPC module
RUN wget https://github.com/dkoboldt/varscan/raw/master/VarScan.v2.4.4.jar \
    && echo '#!/bin/bash' > /opt2/varscan \
    && echo 'java -jar /opt2/VarScan.v2.4.4.jar  "$@"' >> /opt2/varscan \
    && chmod a+rx /opt2/varscan

# Install snpEff/4.3t
# Works with java8
# Deleted bundled version of clinEff, not being used
# Setting env variable to jar file similar to HPC module
RUN wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip \
    && unzip /opt2/snpEff_v4_3t_core.zip \
    && rm /opt2/snpEff_v4_3t_core.zip \
    && rm -rf /opt2/clinEff/
ENV SNPEFF_JAR="/opt2/snpEff/snpEff.jar"

# Install samblaster/0.1.26
# Requires g++ and make (already satisfied)
RUN git clone https://github.com/GregoryFaust/samblaster.git \
    && cd /opt2/samblaster \
    && make \
    && chmod a+rX /opt2/samblaster/samblaster
ENV PATH="/opt2/samblaster:$PATH"
WORKDIR /opt2

# Install strelka/2.9.10, requires python2.7 (already satisfied)
# and set python2.7 as the default interpreter
RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && tar -xvjf /opt2/strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && rm /opt2/strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && ln -fs /usr/bin/python2.7 /usr/bin/python2
ENV PATH="/opt2/strelka-2.9.10.centos6_x86_64/bin:$PATH"

# Install manta/1.6.0
# Requires apt-get packages: bzip2 gcc g++ make python zlib1g-dev (already satisfied)
# Delete demo data (frees 275 MB)
RUN wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 \
    && tar -xvjf /opt2/manta-1.6.0.centos6_x86_64.tar.bz2 \
    && rm /opt2/manta-1.6.0.centos6_x86_64.tar.bz2 \
    && rm -rf /opt2/manta-1.6.0.centos6_x86_64/share/demo/ \
    && ln -fs /usr/bin/python2.7 /usr/bin/python2
ENV PATH="/opt2/manta-1.6.0.centos6_x86_64/bin:$PATH"

# Install VarDict/1.8.3
# Requires java8 (jdk 1.8 or later), perl and R (already satisfied)
RUN wget https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.8.3/VarDict-1.8.3.tar \
    && tar -xvf /opt2/VarDict-1.8.3.tar \
    && rm /opt2/VarDict-1.8.3.tar
ENV PATH="/opt2/VarDict-1.8.3/bin:$PATH"

# Fastp From Opengene github
RUN wget http://opengene.org/fastp/fastp.0.24.0 \
    && mkdir fastp \
    && mv fastp.0.24.0 fastp/fastp \
    && chmod a+x fastp/fastp
ENV PATH="/opt2/fastp:$PATH"

#ASCAT 
RUN Rscript -e 'devtools::install_github("VanLoo-lab/ascat/ASCAT")'

# SvABA
RUN wget -O svaba_1.2.0 https://github.com/walaj/svaba/releases/download/v1.2.0/svaba \
    && mkdir svaba \ 
    && mv svaba_1.2.0 svaba/svaba \
    && chmod a+x svaba/svaba
ENV PATH="/opt2/svaba:$PATH"

# MUSE-Use Cloned Version due to -ldeflate option needed
RUN git clone https://github.com/dnousome/MuSE \
    && cd MuSE \
    && ./install_muse.sh \
    && chmod a+rX /opt2/MuSE/MuSE
ENV PATH="/opt2/MuSE:$PATH"
WORKDIR /opt2

# HMFtools for PURPLE/COBALT/AMBER
RUN wget https://github.com/hartwigmedical/hmftools/releases/download/amber-v4.0.1/amber_v4.0.1.jar \
    && wget https://github.com/hartwigmedical/hmftools/releases/download/cobalt-v1.16/cobalt_v1.16.jar \
    && wget https://github.com/hartwigmedical/hmftools/releases/download/purple-v4.0.2/purple_v4.0.2.jar \
    && wget https://github.com/hartwigmedical/hmftools/releases/download/sage-v3.4.4/sage_v3.4.4.jar \
    && mkdir hmftools \
    && mv amber_v4.0.1.jar hmftools/amber.jar \
    && mv cobalt_v1.16.jar hmftools/cobalt.jar \
    && mv purple_v4.0.2.jar hmftools/purple.jar \
    && mv sage_v3.4.4.jar hmftools/sage.jar \
    && chmod a+x hmftools/amber.jar
ENV PATH="/opt2/hmftools:$PATH"
WORKDIR /opt2

#GLNEXUS
RUN wget https://github.com/dnanexus-rnd/GLnexus/releases/download/v1.4.1/glnexus_cli \
    && chmod +x glnexus_cli

# Add Dockerfile and argparse.bash script
# and export environment variables
ADD Dockerfile /opt2/base_gatk4_wgs.dockerfile
COPY argparse.bash /opt2

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