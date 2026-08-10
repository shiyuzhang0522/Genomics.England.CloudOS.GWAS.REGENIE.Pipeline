# 01_shared_QC

Chromosome-parallel shared genotype quality control for preparing high-quality common SNP genotype datasets from Genomics England AggV3 WGS data for downstream REGENIE Step 1 analysis.

## Overview

This Nextflow pipeline performs the fixed genotype QC steps shared across downstream REGENIE Step 1 analyses.

The workflow runs independently across chromosomes 1–22 and generates chromosome-specific filtered PGEN datasets.

LD pruning is intentionally **not** performed in this pipeline. The resulting high-quality common genotype datasets are intended to be reused by a separate LD-pruning workflow so that different pruning parameters can be tested without repeating the upstream genotype QC.

## QC workflow

For each chromosome, the pipeline:

1. Restricts to the predefined GWAS analysis sample set.
2. Restricts to AggV3 variants passing basic site QC:
   - `FILTER = PASS`
   - `FILTER = LowMLSQ`
3. Retains A/C/G/T SNPs only.
4. Applies minor allele frequency filtering:
   - `MAF >= 0.01`
5. Applies variant missingness filtering:
   - missingness `<= 0.01`
6. Applies Hardy-Weinberg equilibrium filtering:
   - `P >= 1e-15`
7. Writes a new chromosome-specific PGEN/PVAR/PSAM genotype fileset.

The core PLINK2 command is:

```bash
plink2 \
  --pfile dragen \
  --keep GEL_CM_REGENIE.keep.txt \
  --extract basic_QC.variant_ids.txt \
  --snps-only just-acgt \
  --maf 0.01 \
  --geno 0.01 \
  --hwe 1e-15 0 \
  --make-pgen \
  --write-snplist \
  --out chrN.HQ_common
```

## Input data

The pipeline requires three runtime inputs.

### 1. AggV3 genotype root

A directory containing chromosome-specific Genomics England AggV3 PGEN datasets with the following structure:

```text
chrom-msvcf/
├── chrom-1/
│   └── postproc-pgen/
│       ├── dragen.pgen
│       ├── dragen.pvar
│       └── dragen.psam
├── chrom-2/
│   └── postproc-pgen/
│       └── ...
...
└── chrom-22/
    └── postproc-pgen/
        ├── dragen.pgen
        ├── dragen.pvar
        └── dragen.psam
```

Nextflow parameter:

```text
--pgen_root
```

### 2. Basic site-QC variant lists

A directory containing one variant-ID list per chromosome:

```text
chr1.PASS_or_LowMLSQ.variant_ids.txt
chr2.PASS_or_LowMLSQ.variant_ids.txt
...
chr22.PASS_or_LowMLSQ.variant_ids.txt
```

Each file contains one AggV3 variant ID per line.

Nextflow parameter:

```text
--variant_list_dir
```

### 3. Sample keep file

A PLINK2 two-column FID/IID keep file defining the analysis cohort.

The current GEL cutaneous melanoma analysis uses:

```text
GEL_CM_REGENIE.keep.txt
```

Nextflow parameter:

```text
--keep_file
```

## Outputs

For each chromosome, the workflow produces:

```text
chrN.HQ_common.pgen
chrN.HQ_common.pvar
chrN.HQ_common.psam
chrN.HQ_common.snplist
chrN.HQ_common.summary.tsv
chrN.HQ_common.log
```

The per-chromosome summary contains:

```text
chromosome
input_basic_QC_variants
HQ_common_variants
proportion_retained
input_samples
output_samples
```

CloudOS results are organized by chromosome:

```text
results/
├── chr1/
├── chr2/
...
└── chr22/
```

## Parallel execution

The workflow creates one independent `SHARED_QC` process for each autosome:

```text
chr1
chr2
...
chr22
```

These jobs can therefore be executed in parallel by the CloudOS AWS Batch executor.

## Container

The workflow uses the following PLINK2 container:

```text
ghcr.io/shiyuzhang0522/gel-gwas-regenie-plink2:alpha7-20260808
```

Container details:

- Base image: `ubuntu:22.04`
- Platform: `linux/amd64`
- PLINK2: Alpha 7.3
- PLINK2 build date: 2026-08-08

The container definition is maintained separately under:

```text
containers/plink2/
```

## Nextflow files

```text
01_shared_QC/
├── main.nf
├── nextflow.config
└── README.md
```

`main.nf` defines the scientific workflow and chromosome-parallel processing.

`nextflow.config` defines the PLINK2 container, QC thresholds, and process resource requirements.

## CloudOS execution

The workflow is designed for import into Genomics England CloudOS as a Nextflow pipeline.

Required CloudOS parameters:

```text
pgen_root
variant_list_dir
keep_file
```

The current implementation uses Nextflow DSL2 and has been tested with Nextflow `24.04.4` using the AWS Batch executor.

## Downstream workflow

The outputs from this pipeline are used as input to the subsequent LD-pruning stage:

```text
AggV3 WGS
   ↓
01_shared_QC
   ↓
HQ common SNP PGENs
   ↓
02_LD_pruning
   ↓
REGENIE Step 1 variant set
```

Separating genotype QC from LD pruning allows pruning parameters to be tested repeatedly without rerunning the more expensive upstream WGS QC.
