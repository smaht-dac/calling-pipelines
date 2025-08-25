#######################################################################
#     Basic image
#######################################################################
FROM public.ecr.aws/smaht-dac/base-ubuntu2204-py38:0.0.1

LABEL maintainers="Michele Berselli (berselli.michele@gmail.com)"

#######################################################################
#     General updates & installing necessary Linux components
#     VEP:
#       HTSlib
#           - libbz2-dev, liblzma-dev
#       Perl
#           - (utils) libmodule-build-perl, liblist-moreutils-perl, libwww-perl, bioperl
#           - (Archive::Zip) libarchive-zip-perl
#           - (DBI) libdbi-perl
#           - (test) libtest-exception-perl, libtest-deep-perl, libtest-warnings-perl
#       Other
#           - (DBD::mysql) libdbd-mysql-perl
#######################################################################
RUN apt-get update -y && apt-get install -y \
    locales \
    tabix=1.13+ds-2build1 \
    libbz2-dev=1.0.8-5build1 \
    liblzma-dev=5.2.5-2ubuntu1 \
    bioperl=1.7.8-1 \
    libmodule-build-perl=0.423100-1 \
    liblist-moreutils-perl=0.430-2 \
    libwww-perl=6.61-1 \
    libarchive-zip-perl=1.68-1  \
    libdbi-perl=1.643-3build3 \
    libtest-exception-perl=0.43-1 \
    libtest-deep-perl=1.130-1 \
    libtest-warnings-perl=0.031-1 \
    libdbd-mysql-perl=4.050-5ubuntu0.22.04.1 \
    --no-install-recommends \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

#######################################################################
#     Setting env
#######################################################################
## Set working directory
WORKDIR /usr/local/bin

## Supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

#######################################################################
#     Software
#######################################################################
## VEP install
RUN git clone https://github.com/Ensembl/ensembl-vep.git
WORKDIR ensembl-vep
RUN git checkout release/114 && \
    perl INSTALL.pl --AUTO a --NO_UPDATE

WORKDIR /usr/local/bin

## Bcftools install
RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
    tar -xjf bcftools-1.17.tar.bz2 && \
    rm -rf bcftools-1.17.tar.bz2

WORKDIR /usr/local/bin/bcftools-1.17
RUN ./configure && \
    make

WORKDIR /usr/local/bin/
RUN ln -s bcftools-1.17 bcftools

#######################################################################
#     Script
#######################################################################
## vep-parallel
COPY vep-parallel.sh .
RUN chmod +x vep-parallel.sh

######################################################################
#     ENV
#######################################################################
ENV PATH=/usr/local/bin/ensembl-vep/:$PATH
ENV PATH=/usr/local/bin/bcftools/:$PATH

CMD ["bash"]
