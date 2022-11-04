# Docker inheritance
FROM bioconductor/bioconductor_docker:devel

# Update apt-get
RUN apt-get update -y

#install python
RUN apt-get install -y python-htseq

# Install required Bioconductor package
#RUN R -e 'BiocManager::install(c("biomaRt", "GenomicFeatures", "DEXSeq"))'

#Install required R package
#RUN Rscript -e 'install.packages("tidyverse")'

