#######################################################################
#     Basic image
#######################################################################
FROM public.ecr.aws/smaht-dac/base-ubuntu2204-py38:0.0.1

LABEL maintainers="Charlotte Zhang (charlotte_zhang@hms.harvard.edu), Michele Berselli (berselli.michele@gmail.com)"

#######################################################################
#     General updates & installing necessary Linux components
#       - libncurses-dev, libbz2-dev, liblzma-dev used by samtools
#######################################################################
RUN apt-get update -y && apt-get install -y \
    locales \
    tabix=1.13+ds-2build1 \
    libncurses-dev=6.3-2 \
    libbz2-dev=1.0.8-5build1 \
    liblzma-dev=5.2.5-2ubuntu1 \
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
## Bcftools install
RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
    tar -xjf bcftools-1.17.tar.bz2 && \
    rm -rf bcftools-1.17.tar.bz2

WORKDIR /usr/local/bin/bcftools-1.17
RUN ./configure && \
    make

WORKDIR /usr/local/bin/
RUN ln -s bcftools-1.17 bcftools

## Delly install
RUN wget https://github.com/dellytools/delly/releases/download/v1.5.0/delly_v1.5.0_linux_x86_64bit \
    && chmod +x delly_v1.5.0_linux_x86_64bit \
    && ln -s delly_v1.5.0_linux_x86_64bit delly

#######################################################################
#     Scripts
#######################################################################
## run_Delly_lr
COPY run_delly_lr.sh .
RUN chmod +x run_delly_lr.sh

#######################################################################
#     Updating env
#######################################################################
ENV PATH=/usr/local/bin/:$PATH
ENV PATH=/usr/local/bin/bcftools/:$PATH

CMD ["bash"]
