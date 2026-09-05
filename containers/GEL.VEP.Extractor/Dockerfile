FROM condaforge/miniforge3:24.3.0-0


LABEL maintainer="Shiyu Zhang"
LABEL description="GEL VEP extraction tool environment"
LABEL org.opencontainers.image.source="https://github.com/shiyuzhang0522/Genomics.England.CloudOS.GWAS.REGENIE.Pipeline"


##############################################
# Create software environment
##############################################

RUN conda create -y \
    -n gel_vep_extractor \
    -c conda-forge \
    -c bioconda \
    python=3.10 \
    bcftools=1.20 \
    htslib \
    tabix \
    bedtools \
    ucsc-bigwigtobedgraph \
    pandas \
    numpy \
    && conda clean -afy


##############################################
# Activate environment globally
##############################################

ENV CONDA_DEFAULT_ENV=gel_vep_extractor

ENV PATH=/opt/conda/envs/gel_vep_extractor/bin:$PATH


##############################################
# Ensure interactive shell activates env
##############################################

RUN echo "conda activate gel_vep_extractor" >> /root/.bashrc


##############################################
# Validate installation during build
##############################################

RUN which python && \
    python --version && \
    which bcftools && \
    bcftools --version && \
    which bedtools && \
    bedtools --version && \
    which bigWigToBedGraph && \
    bigWigToBedGraph 2>&1 | head -n 2 && \
    python - <<'EOF'
import pandas
import numpy

print("Python environment OK")
print("pandas:", pandas.__version__)
print("numpy:", numpy.__version__)
EOF


##############################################
# Working directory
##############################################

WORKDIR /data


##############################################
# Default command
##############################################

CMD ["bash", "--login"]
