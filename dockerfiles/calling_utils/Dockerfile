#######################################################################
#     Basic image
#######################################################################
FROM public.ecr.aws/smaht-dac/base-ubuntu2204-py39:0.0.1

LABEL maintainers="Michele Berselli (berselli.michele@gmail.com)"

#######################################################################
#     General updates & installing necessary Linux components
#       - libncurses-dev, libbz2-dev, liblzma-dev used by samtools
#######################################################################
RUN apt-get update -y && apt-get install -y \
    locales \
    tabix=1.13+ds-2build1 \
    libncurses-dev=6.3-2ubuntu0.1 \
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
## Conda install
#- pysam required by filter_by_poe.py
RUN conda install -c bioconda -y pysam==0.21.0 && \
    conda clean -a -y -f

## granite
RUN pip install granite-suite==0.5.0

## Samtools install
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar -xjf samtools-1.17.tar.bz2 && \
    rm -rf samtools-1.17.tar.bz2

WORKDIR /usr/local/bin/samtools-1.17
RUN ./configure && \
    make

WORKDIR /usr/local/bin/
RUN ln -s samtools-1.17 samtools

## Bcftools install
RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
    tar -xjf bcftools-1.17.tar.bz2 && \
    rm -rf bcftools-1.17.tar.bz2

WORKDIR /usr/local/bin/bcftools-1.17
RUN ./configure && \
    make

WORKDIR /usr/local/bin/
RUN ln -s bcftools-1.17 bcftools

## minipileup
RUN wget https://github.com/lh3/minipileup/archive/refs/tags/v1.4b.tar.gz && \
    tar -xzf v1.4b.tar.gz && \
    rm -rf v1.4b.tar.gz

WORKDIR /usr/local/bin/minipileup-1.4b
RUN make

WORKDIR /usr/local/bin/
RUN ln -s minipileup-1.4b minipileup

#######################################################################
#     Scripts
#######################################################################
## bcftools_PASS_norm_dedup
COPY bcftools_PASS_norm_dedup.sh .
RUN chmod +x bcftools_PASS_norm_dedup.sh

## bcftools_regions
COPY bcftools_regions.sh .
RUN chmod +x bcftools_regions.sh

## compress_index_vcf
COPY compress_index_vcf.sh .
RUN chmod +x compress_index_vcf.sh

## filter_panel_errors
COPY filter_panel_errors.sh .
RUN chmod +x filter_panel_errors.sh

## split_snvs_indels
COPY split_snvs_indels.sh .
RUN chmod +x split_snvs_indels.sh

## minipileup-parallel
COPY minipileup-parallel.sh .
RUN chmod +x minipileup-parallel.sh

## scripts/
# filter_by_poe.py
COPY scripts/filter_by_poe.py .
RUN chmod +x filter_by_poe.py

# filter_clustered_variants.py
COPY scripts/filter_clustered_variants.py .
RUN chmod +x filter_clustered_variants.py

# tier_filter_variants_SR_PB_ONT.py
COPY scripts/tier_filter_variants_SR_PB_ONT.py .
RUN chmod +x tier_filter_variants_SR_PB_ONT.py

#######################################################################
#     Update env
#######################################################################
ENV PATH=/usr/local/bin/:$PATH
ENV PATH=/usr/local/bin/samtools/:$PATH
ENV PATH=/usr/local/bin/bcftools/:$PATH
ENV PATH=/usr/local/bin/minipileup/:$PATH

CMD ["bash"]
