FROM --platform=linux/x86_64 biocontainer/plink2:alpha2.3_jan2020

ARG SAMTOOLS_VER="1.20"

WORKDIR /

# LABEL instructions tag the image with metadata that might be important to the user
# Optional, but highly recommended
LABEL base.image="biocontainer/plink2"
LABEL dockerfile.version="1"
LABEL software="plink2"
LABEL maintainer="Will Gardner"
LABEL maintainer.email="wkgardner@wisc.edu"

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    build-essential \
    ca-certificates \
    procps \
    autoconf \
    autotools-dev \
    automake \
    python3 \
    bedtools \
    python3-pip \
    vim \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    zlib1g-dev \
    libssl-dev \
    gcc \
    make \
    ca-certificates \
    zlib1g \
    bzip2 \
    gnuplot \
    procps && \

# download, compile, and install samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER} && \
    ./configure && \
    make && \
    make install && \
    make test
