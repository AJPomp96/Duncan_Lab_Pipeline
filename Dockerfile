# Docker inheritance
FROM bioconductor/bioconductor_docker:devel

# Update apt-get
RUN apt-get update

# Install required Bioconductor package
RUN R -e 'BiocManager::install(c("biomaRt", "GenomicFeatures", "DEXseq"))'

#Install required R package
RUN Rscript -e 'install.packages("tidyverse")'

