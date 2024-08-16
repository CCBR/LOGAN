FROM --platform=linux/amd64 gridss/gridss:2.13.2

LABEL maintainer=<nousomedr@nih.gov>

# Create Container filesystem specific 
# working directory and opt directories 
WORKDIR /opt2 

##Add GRIPSS for SOMATIC FILTERING
RUN wget https://github.com/hartwigmedical/hmftools/releases/download/gripss-v2.3.4/gripss_v2.3.4.jar \
    && mkdir hmftools \
    && mv gripss_v2.3.4.jar hmftools/gripss.jar 
ENV PATH="/opt2/hmftools:$PATH"


