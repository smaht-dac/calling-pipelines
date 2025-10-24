# ---------- Build stage ----------
FROM ubuntu:16.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive VER_STRELKA2=2.9.10

# Try default mirrors; if update fails, fall back to old-releases
RUN set -eux; \
    if ! apt-get update -y; then \
      sed -i 's|archive.ubuntu.com|old-releases.ubuntu.com|g; s|security.ubuntu.com|old-releases.ubuntu.com|g' /etc/apt/sources.list; \
      apt-get update -y; \
    fi && \
    apt-get install -y --no-install-recommends \
      ca-certificates curl gcc g++ make cmake zlib1g-dev libboost-dev bzip2 python && \
    rm -rf /var/lib/apt/lists/*

# Download and install Strelka2
WORKDIR /tmp
RUN curl -fsSL --retry 5 \
      "https://github.com/Illumina/strelka/releases/download/v${VER_STRELKA2}/strelka-${VER_STRELKA2}.release_src.tar.bz2" \
      -o strelka.tar.bz2 && \
    mkdir strelka-src && \
    tar -xjf strelka.tar.bz2 -C strelka-src --strip-components=1 && \
    rm strelka.tar.bz2

WORKDIR /tmp/strelka-src/build
RUN ../configure --jobs=$(nproc) --prefix=/opt/strelka && \
    make -j"$(nproc)" install && \
    # Verify installation
    /opt/strelka/bin/configureStrelkaGermlineWorkflow.py --help >/dev/null && \
    /opt/strelka/bin/configureStrelkaSomaticWorkflow.py --help >/dev/null

# ---------- Runtime stage ----------
FROM ubuntu:16.04
ENV DEBIAN_FRONTEND=noninteractive

LABEL maintainers="Michele Berselli (berselli.michele@gmail.com)"

RUN set -eux; \
    if ! apt-get update -y; then \
      sed -i 's|archive.ubuntu.com|old-releases.ubuntu.com|g; s|security.ubuntu.com|old-releases.ubuntu.com|g' /etc/apt/sources.list; \
      apt-get update -y; \
    fi && \
    apt-get install -y --no-install-recommends \
      ca-certificates bzip2 zlib1g python curl tabix gzip && \
    rm -rf /var/lib/apt/lists/*
# tabix will install bgzip as a dependency

# Copy Strelka2 from the build stage
COPY --from=builder /opt/strelka /opt/strelka
ENV PATH="/opt/strelka/bin:${PATH}"

# Set working directory
WORKDIR /usr/local/bin

# Import scripts
COPY strelka2_region.sh .
RUN chmod +x strelka2_region.sh

CMD ["bash"]
