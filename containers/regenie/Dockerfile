FROM condaforge/miniforge3:24.3.0-0

LABEL maintainer="Shiyu Zhang"
LABEL description="REGENIE 4.1.2 GWAS software environment"
LABEL org.opencontainers.image.source="https://github.com/shiyuzhang0522/Genomics.England.CloudOS.GWAS.REGENIE.Pipeline"

RUN conda create -y \
    -n regenie \
    -c conda-forge \
    -c bioconda \
    regenie=4.1.2 \
    && conda clean -afy

ENV PATH="/opt/conda/envs/regenie/bin:${PATH}"

CMD ["regenie", "--help"]
