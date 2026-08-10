# Genomics.England.CloudOS.GWAS.REGENIE.Pipeline
Nextflow pipelines for Genomics England CloudOS GWAS using REGENIE, including cohort-specific genotype QC, common-variant selection, LD pruning, and Step 1/2 analysis workflows.

## PLINK2 container

This repository includes a minimal PLINK2 container for CloudOS/Nextflow smoke testing.

- Base image: `ubuntu:22.04`
- PLINK2 archive: `https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260808.zip`
- Image tag: `ghcr.io/shiyuzhang0522/gel-gwas-regenie-plink2:alpha7-20260808`

Build the Linux amd64 image locally:

```bash
docker buildx build \
  --platform linux/amd64 \
  -t ghcr.io/shiyuzhang0522/gel-gwas-regenie-plink2:alpha7-20260808 \
  --load \
  containers/plink2
```

Verify PLINK2 inside the container:

```bash
docker run --rm \
  --platform linux/amd64 \
  ghcr.io/shiyuzhang0522/gel-gwas-regenie-plink2:alpha7-20260808 \
  plink2 --version
```

Push the image to GitHub Container Registry:

```bash
docker buildx build \
  --platform linux/amd64 \
  -t ghcr.io/shiyuzhang0522/gel-gwas-regenie-plink2:alpha7-20260808 \
  --push \
  containers/plink2
```

## Nextflow smoke test

The initial Nextflow workflow runs `plink2 --version` using the PLINK2 container.

```bash
nextflow run main.nf -profile docker
```

`nextflow` is not required to build the Docker image, but it is required to run the smoke workflow.
