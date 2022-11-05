# Docker inheritance
FROM bioconductor/bioconductor_docker:devel

ARG BBTOOLSVER=38.96

# Update apt-get
RUN apt-get update -y

#install python and HTSeq
RUN apt-get install -y python-htseq

#install bbmap
RUN wget --progress=dot:giga https://sourceforge.net/projects/bbmap/files/BBMap_${BBTOOLSVER}.tar.gz && \
    tar -xzf BBMap_${BBTOOLSVER}.tar.gz && \
    rm BBMap_${BBTOOLSVER}.tar.gz && \
    ln -s /opt/bbmap/*.sh /usr/local/bin/

#install fastx_toolkit pre-compiled binaries
RUN mkdir fastx_bin && \
    cd fastx_bin && \
    wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && \
    tar -xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && \
    sudo cp ./bin/* /usr/local/bin

#add bbmap to the path
ENV PATH="${PATH}:/bbmap"

# Install required Bioconductor package
#RUN R -e 'BiocManager::install(c("biomaRt", "GenomicFeatures", "DEXSeq"))'

#Install required R package
#RUN Rscript -e 'install.packages("tidyverse")'

