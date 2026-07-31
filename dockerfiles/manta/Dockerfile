#######################################################################
#     Basic image
#######################################################################
FROM public.ecr.aws/smaht-dac/base-ubuntu2204-py310:0.0.1

LABEL maintainers="Alexander Veit (alexander_veit@hms.harvard.edu)"

#######################################################################
#     General updates & installing necessary Linux components
#######################################################################
RUN apt-get update -y && apt-get install -y \
    locales \
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
## Manta install
#- Manta 1.6.0 requires Python 2.7, so it is installed into a dedicated
#- conda environment to avoid downgrading the base (py310) environment
# Accept ToS
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install
RUN conda create -y -n manta -c bioconda \
    manta==1.6.0 \
    && conda clean -a -y -f

ENV PATH=/opt/conda/envs/manta/bin:$PATH

#######################################################################
#     Scripts
#######################################################################
## run_manta
COPY run_manta.sh .
RUN chmod +x run_manta.sh

CMD ["bash"]
