FROM --platform=linux/amd64 nciccbr/ccbr_ubuntu_base_20.04:v6

LABEL maintainer=<nousomedr@nih.gov>

WORKDIR /opt2 

RUN apt-get update 

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
      gnupg \
      dirmngr \
      ca-certificates \
      apt-transport-https \
      software-properties-common 

#Install R for ASCAT
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && add-apt-repository --enable-source --yes 'ppa:c2d4u.team/c2d4u4.0+' \
    && apt-get -y install r-base r-base-core r-recommended r-base-dev \
    && apt-get -y install libcurl4-openssl-dev libssl-dev libboost-dev libxml2-dev \
    && apt-get -y install r-cran-biocmanager r-cran-devtools r-bioc-genomicranges 

#ASCAT 
RUN Rscript -e 'devtools::install_github("VanLoo-lab/ascat/ASCAT")' 
RUN Rscript -e 'install.packages(c("argparse"), repos="http://cran.r-project.org")'
RUN Rscript -e 'BiocManager::install("DNAcopy")'

#Allelecounter
RUN git clone https://github.com/cancerit/alleleCount ac \
&& cd ac \
&& ./setup.sh /opt2/alleleCount 
ENV PATH="/opt2/alleleCount/bin:$PATH"
ENV LD_LIBARY_PATH="/opt2/alleleCount/lib"

#CNVKIT
RUN git clone https://github.com/etal/cnvkit \
&& cd cnvkit/ \
&& pip3 install --upgrade pip \
&& pip install -e .

##Clean up folders
WORKDIR /opt2
RUN rm -R ac


