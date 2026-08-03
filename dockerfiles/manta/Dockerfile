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
#- Conda lives at /miniconda3 in the base image, so the env is /miniconda3/envs/manta
#- libboost must be pinned to 1.85: the bioconda manta package declares an
#- unpinned `libboost` dependency, but its libexec binaries are linked against
#- libboost_{filesystem,program_options,serialization}.so.1.85.0. Without the
#- pin the solver installs a newer libboost and every libexec binary fails at
#- runtime with "error while loading shared libraries"
RUN conda create -y -p /miniconda3/envs/manta -c bioconda -c conda-forge \
    manta==1.6.0 \
    libboost==1.85.0 \
    && conda clean -a -y -f \
    && test -x /miniconda3/envs/manta/bin/configManta.py

ENV PATH=/miniconda3/envs/manta/bin:$PATH

## Verify that all Manta binaries have their shared libraries resolved,
## so a dependency mismatch fails the build instead of the pipeline run
RUN set -eu; \
    for f in /miniconda3/envs/manta/libexec/*; do \
      if ldd "$f" 2>/dev/null | grep -q "not found"; then \
        echo "ERROR: unresolved shared libraries in $f"; ldd "$f"; exit 1; \
      fi; \
    done

#######################################################################
#     Scripts
#######################################################################
## run_manta
COPY run_manta.sh .
RUN chmod +x run_manta.sh

CMD ["bash"]
