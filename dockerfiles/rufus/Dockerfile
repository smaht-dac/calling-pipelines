FROM ubuntu:22.04

# Metadata labels
LABEL author="Stephanie Georges, Michele Berselli (berselli.michele@gmail.com)"
LABEL maintainers="Michele Berselli (berselli.michele@gmail.com)"

# Set environment variables for installation
ENV DEBIAN_FRONTEND=noninteractive
ENV HTSLIB_VERSION="1.21"

# Install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        locales \
        python3 \
        git \
        cmake \
        wget \
        g++ \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        bc \
        libgsl0-dev \
        libncurses5-dev \
        autoconf \
        automake \
        make \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        vt \
        parallel \
        gawk \
        libjsoncpp-dev \
        libjsoncpp25 \
        curl \
        unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install htslib
RUN cd /opt && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git --branch "$HTSLIB_VERSION" && \
    cd htslib && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install

# Install samtools
RUN cd /opt && \
    git clone https://github.com/samtools/samtools.git --branch "$HTSLIB_VERSION" && \
    cd samtools && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure --with-htslib=/opt/htslib/ && \
    make && \
    make install

# Install bcftools
RUN cd /opt && \
    git clone https://github.com/samtools/bcftools.git --branch "$HTSLIB_VERSION" && \
    cd bcftools && \
    autoheader && \
    autoconf && \
    ./configure --with-htslib=/opt/htslib/ --enable-libgsl && \
    make && \
    make install

# Install bamtools
RUN cd /opt && \
    git clone https://github.com/pezmaster31/bamtools.git --branch v2.5.2 && \
    cd bamtools/ && \
    mkdir build && \
    cd build/ && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/bamtools/build .. && \
    make && \
    make install

# Install bedtools
RUN cd /opt && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -zxvf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make && \
    make install

# Install RUFUS
RUN cd /opt && \
    git clone https://github.com/stefinfection/RUFUS.git -b d1.1.1 && \
    cd RUFUS && \
    mkdir -p bin && \
    cd bin && \
    cmake ../ && \
    make

# Unset DEBIAN_FRONTEND
ENV DEBIAN_FRONTEND=

# Set runtime environment variables
ENV RUFUS_ROOT=/opt/RUFUS
ENV PATH=/usr/local/bin:$PATH
ENV PATH=/opt/RUFUS:/opt/RUFUS/bin:$PATH
ENV PATH=/opt/bamtools/build/bin:$PATH
ENV PATH=/opt/bedtools2/bin:$PATH

RUN locale-gen "en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8
ENV LC_CTYPE=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US.UTF-8

## Set working directory
WORKDIR /usr/local/bin

## Scripts
# run_rufus
COPY run_rufus.sh .
RUN chmod +x run_rufus.sh

# Default command
CMD ["/bin/bash"]
