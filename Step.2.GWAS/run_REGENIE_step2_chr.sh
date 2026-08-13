#!/usr/bin/env bash
set -euo pipefail

################################################################################
#
# REGENIE Step 2 Genome-wide Association Testing
#
# Phenotype:
#   Cutaneous melanoma (CM)
#
# Dataset:
#   Genomics England (GEL) AggV3 WGS
#
# Execution:
#   Interactive CloudOS session
#
# Description:
#   Chromosome-specific REGENIE Step 2 GWAS
#
# Model:
#   Binary logistic regression
#   Additive model
#   LOCO prediction adjustment
#   Approximate Firth correction
#
# Author:
#   Shelley
#
# Date:
#   2026-08-13
#
################################################################################


################################################################################
# Arguments
################################################################################

if [[ $# -ne 1 ]]; then
    echo
    echo "Usage:"
    echo "  bash run_REGENIE_step2_chr.sh <chromosome>"
    echo
    echo "Example:"
    echo "  bash run_REGENIE_step2_chr.sh 1"
    echo
    exit 1
fi


CHR=$1


if ! [[ "${CHR}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "ERROR: chromosome must be 1-22"
    exit 1
fi



################################################################################
# Working directory
################################################################################

WORKDIR="/home/vscode/session_data/REGENIE.GWAS/GWAS.pipeline/REGENIE.Step2.GWAS"

cd "${WORKDIR}"



################################################################################
# Input paths
################################################################################


# WGS genotype

PGEN_ROOT="/home/vscode/session_data/filesystems/chrom-msvcf"


# QCed variant lists

VARIANT_LIST_DIR="/home/vscode/session_data/filesystems/QCed.WGS.Var.list"


# Phenotype

PHENO_FILE="/home/vscode/session_data/filesystems/GEL_CM_REGENIE.phenotype.tsv"


# Covariates

COVAR_FILE="/home/vscode/session_data/filesystems/GEL_CM_REGENIE.covariates.tsv"


# Step 1 LOCO prediction list
#
# Manually generated:
#
# CM /home/vscode/session_data/mounted-data-readonly/GEL_CM_REGENIE_step1_1.loco

PRED_FILE="${WORKDIR}/GEL_CM_REGENIE_step1_pred.list"



################################################################################
# Output
################################################################################


RESULT_DIR="${WORKDIR}/results/chr${CHR}"
LOG_DIR="${WORKDIR}/logs"

mkdir -p "${RESULT_DIR}"
mkdir -p "${LOG_DIR}"


LOG_FILE="${LOG_DIR}/chr${CHR}.REGENIE.step2.log"



################################################################################
# Logging
################################################################################


exec > >(tee -a "${LOG_FILE}") 2>&1


echo "============================================================"
echo "REGENIE Step 2 GWAS"
echo "============================================================"

date

echo
echo "Chromosome:"
echo "${CHR}"

echo
echo "Hostname:"
hostname

echo
echo "Working directory:"
pwd



################################################################################
# Environment
################################################################################


echo
echo "Environment information"
echo "------------------------------------------------------------"

conda info --envs || true

echo

regenie --version



################################################################################
# Hardware
################################################################################


echo
echo "CPU information"
echo "------------------------------------------------------------"

nproc


echo
echo "Memory information"
echo "------------------------------------------------------------"

free -h



################################################################################
# Input files
################################################################################


PGEN_PREFIX="${PGEN_ROOT}/chrom-${CHR}/postproc-pgen/dragen"

VARIANT_LIST="${VARIANT_LIST_DIR}/chr${CHR}.PASS_or_LowMLSQ.variant_ids.txt"



echo
echo "Input files"
echo "------------------------------------------------------------"

echo "PGEN:"
echo "${PGEN_PREFIX}"


echo "Variant list:"
echo "${VARIANT_LIST}"


echo "Phenotype:"
echo "${PHENO_FILE}"


echo "Covariates:"
echo "${COVAR_FILE}"


echo "Prediction list:"
echo "${PRED_FILE}"



################################################################################
# Input validation
################################################################################


for FILE in \
    "${PGEN_PREFIX}.pgen" \
    "${PGEN_PREFIX}.pvar" \
    "${PGEN_PREFIX}.psam" \
    "${VARIANT_LIST}" \
    "${PHENO_FILE}" \
    "${COVAR_FILE}" \
    "${PRED_FILE}"

do

    if [[ ! -s "${FILE}" ]]; then

        echo
        echo "ERROR: missing or empty file:"
        echo "${FILE}"

        exit 1

    fi

done



################################################################################
# Save command
################################################################################


COMMAND_FILE="${RESULT_DIR}/chr${CHR}.REGENIE.command.txt"



cat > "${COMMAND_FILE}" << EOF
regenie \\
    --step 2 \\
    --pgen ${PGEN_PREFIX} \\
    --extract ${VARIANT_LIST} \\
    --phenoFile ${PHENO_FILE} \\
    --phenoCol CM \\
    --covarFile ${COVAR_FILE} \\
    --catCovarList genetic_sex,study_source \\
    --maxCatLevels 30 \\
    --pred ${PRED_FILE} \\
    --bt \\
    --minMAC 20 \\
    --firth \\
    --approx \\
    --firth-se \\
    --pThresh 0.01 \\
    --bsize 400 \\
    --threads 16 \\
    --gz \\
    --out ${RESULT_DIR}/chr${CHR}.GEL_CM_REGENIE_step2
EOF



################################################################################
# Run REGENIE Step 2
################################################################################


cd "${RESULT_DIR}"


echo
echo "Starting REGENIE"
echo "------------------------------------------------------------"



regenie \
    --step 2 \
    --pgen "${PGEN_PREFIX}" \
    --extract "${VARIANT_LIST}" \
    --phenoFile "${PHENO_FILE}" \
    --phenoCol CM \
    --covarFile "${COVAR_FILE}" \
    --catCovarList genetic_sex,study_source \
    --maxCatLevels 30 \
    --pred "${PRED_FILE}" \
    --bt \
    --minMAC 20 \
    --firth \
    --approx \
    --firth-se \
    --pThresh 0.01 \
    --bsize 400 \
    --threads 16 \
    --gz \
    --out "${RESULT_DIR}/chr${CHR}.GEL_CM_REGENIE_step2"



################################################################################
# Output check
################################################################################


echo
echo "Checking output"
echo "------------------------------------------------------------"


GWAS_OUTPUT="${RESULT_DIR}/chr${CHR}.GEL_CM_REGENIE_step2_CM.regenie.gz"


if [[ ! -s "${GWAS_OUTPUT}" ]]; then

    echo
    echo "ERROR:"
    echo "Missing REGENIE output:"
    echo "${GWAS_OUTPUT}"

    exit 1

fi



echo
echo "SUCCESS"
echo "Output:"
ls -lh "${GWAS_OUTPUT}"


date

echo "============================================================"
echo "chr${CHR} completed successfully"
echo "============================================================"