FROM docker.io/bioconductor/bioconductor_docker:RELEASE_3_16

LABEL org.opencontainers.image.source=https://github.com/sipss/AlpsNMR
LABEL org.opencontainers.image.description="AlpsNMR preinstalled"
LABEL org.opencontainers.image.licenses=MIT


RUN Rscript -e 'BiocManager::install("AlpsNMR")'

