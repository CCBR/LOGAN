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

###Create AnnotSV 
RUN wget https://github.com/lgmgeo/AnnotSV/archive/refs/tags/v3.3.6.tar.gz \
    && tar -xvzf /opt2/v3.3.6.tar.gz \
    && rm /opt2/v3.3.6.tar.gz
ENV PATH="/opt2/AnnotSV-3.3.6/bin:$PATH"

##ClassifyCNV
##Update the resources for ClassifyCNV
RUN wget https://github.com/Genotek/ClassifyCNV/archive/refs/tags/v1.1.1.tar.gz \
    && tar -xvzf /opt2/v1.1.1.tar.gz \
    && rm /opt2/v1.1.1.tar.gz 
    #&& chmod a+rx /opt2/ClassifyCNV-1.1.1/update_clingen.sh 
ENV PATH="/opt2/ClassifyCNV-1.1.1/:$PATH"
#RUN update_clingen.sh

##Survivor 
RUN wget https://github.com/fritzsedlazeck/SURVIVOR/archive/refs/tags/v1.0.6.tar.gz \
    && tar -xvzf /opt2/v1.0.6.tar.gz \
    && rm /opt2/v1.0.6.tar.gz \
    && cd /opt2/SURVIVOR-1.0.6/Debug \
    && make
ENV PATH="/opt2/SURVIVOR-1.0.6/:$PATH"


COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

