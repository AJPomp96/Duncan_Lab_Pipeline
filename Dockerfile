# Docker inheritance
FROM bioconductor/bioconductor_docker:devel

#define bbtools version
ARG BBTOOLSVER=38.96

# Update apt-get
RUN apt-get update -y

# Add 3.7 to the available alternatives
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1

# Set python3.7 as the default python
RUN update-alternatives --set python /usr/bin/python3.7

#install HTSeq
#RUN sudo apt-get -y install build-essential python3.6-dev python-numpy python-matplotlib python-pysam python-htseq
RUN pip install HTSeq

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

#install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar jxf samtools-1.9.tar.bz2 && \
	rm samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix $(pwd) && \
	make

# Install required Bioconductor package
RUN R -e 'BiocManager::install(c("biomaRt", "GenomicFeatures", "DEXSeq"))'

#Install required R package
RUN R -e 'install.packages(c("tidyverse", "XML"))'

#add dexseq py scripts to path
ENV DEXSEQPATH="/usr/local/lib/R/site-library/DEXSeq/python_scripts"
ENV PATH="${PATH}:${DEXSEQPATH}"
RUN sed -i '1 i #!/usr/bin/env python' "${DEXSEQPATH}/dexseq_count.py"
RUN sed -i '1 i #!/usr/bin/env python' "${DEXSEQPATH}/dexseq_prepare_annotation.py"
RUN chmod +x "${DEXSEQPATH}/dexseq_count.py"
RUN chmod +x "${DEXSEQPATH}/dexseq_prepare_annotation.py"

#add samtools to the path
ENV PATH=/samtools-1.9:${PATH}

#add bbmap to the path
ENV PATH="${PATH}:/bbmap"
