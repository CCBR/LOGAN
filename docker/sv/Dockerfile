FROM --platform=linux/amd64 nciccbr/ccbr_ubuntu_base_20.04:v6

LABEL maintainer=<nousomedr@nih.gov>

WORKDIR /opt2 

RUN apt-get update 

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
      gnupg \
      dirmngr \
      ca-certificates \
      apt-transport-https \
      software-properties-common \
      openjdk-17-jdk   

# Create Container filesystem specific 
# working directory and opt directories 

##Install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && add-apt-repository --enable-source --yes 'ppa:c2d4u.team/c2d4u4.0+' \
    && apt-get -y install r-base r-base-core r-recommended r-base-dev \
    && apt-get -y install libcurl4-openssl-dev libssl-dev libboost-dev libxml2-dev 


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

#Grab GRIDSS
RUN mkdir gridss 
WORKDIR /opt2/gridss
RUN wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    && wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss \
    && wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss.config.R \
    && wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/libgridss.R \
    && wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss_somatic_filter 

ENV GRIDSS_JAR=/opt/gridss2/gridss-2.13.2-gridss-jar-with-dependencies.jar

RUN chmod +x /opt2/gridss/* \
    && chmod +x /opt2/gridss/*.R

WORKDIR /opt2 


##Add GRIPSS for SOMATIC FILTERING
RUN wget https://github.com/hartwigmedical/hmftools/releases/download/gripss-v2.3.4/gripss_v2.3.4.jar \
    && mkdir hmftools \
    && mv gripss_v2.3.4.jar hmftools/gripss.jar 
ENV PATH="/opt2/gridss:/opt2/hmftools:$PATH"

