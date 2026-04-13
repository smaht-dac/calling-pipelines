#######################################################################
#     Basic image
#######################################################################
FROM public.ecr.aws/smaht-dac/base-ubuntu2204-py310:0.0.1

LABEL maintainers="William Feng (william_feng@hms.harvard.edu), Michele Berselli (berselli.michele@gmail.com)"

# #######################################################################
# #     General updates & installing necessary Linux components
# #######################################################################
RUN apt-get update -y && apt-get install -y \
    locales \
    --no-install-recommends \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# #######################################################################
# #     Setting env
# #######################################################################
## Set working directory
WORKDIR /usr/local/bin

## Supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

# #######################################################################
# #     Software
# #######################################################################

## kanpig install
# Accept ToS
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install
RUN conda install -y -c bioconda \
    kanpig==1.1.0 \
    htslib==1.3.1 \
    && conda clean -a -y -f

# Import scripts
COPY run_kanpig_plup.sh .
RUN chmod +x run_kanpig_plup.sh

CMD ["bash"]
