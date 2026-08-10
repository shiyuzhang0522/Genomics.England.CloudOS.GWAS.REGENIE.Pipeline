FROM ubuntu:22.04

ARG PLINK2_URL=https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260808.zip

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        unzip \
        wget && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/plink2

RUN wget "${PLINK2_URL}" -O plink2.zip && \
    unzip plink2.zip && \
    chmod +x plink2 && \
    mv plink2 /usr/local/bin/plink2 && \
    rm plink2.zip

RUN plink2 --version

WORKDIR /work
