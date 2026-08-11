#!/usr/bin/env bash

set -euo pipefail


################################################################################
# Extract LD-pruned variants from chromosome-specific HQ_common PLINK2 datasets
# for REGENIE Step 1 genotype preparation
#
# Input:
#   HQ_common PLINK2 datasets:
#     /home/vscode/session_data/filesystems/results/chrN/chrN.HQ_common.*
#
#   LD-pruned variant lists:
#     /home/vscode/session_data/filesystems/chrN/chrN.LD_pruned.prune.in
#
# Output:
#   Extracted.HQ.Pruned.Var.GT.Step.1/
#     chrN/
#       chrN.HQ_pruned.pgen
#       chrN.HQ_pruned.pvar
#       chrN.HQ_pruned.psam
#
# Author: Shelley
# Date: 2026-08-11
################################################################################


############################
# Paths
############################

FILESYSTEM_DIR="/home/vscode/session_data/filesystems"

INPUT_DIR="${FILESYSTEM_DIR}/results"

PRUNE_DIR="${FILESYSTEM_DIR}"

OUTPUT_DIR="/home/vscode/session_data/REGENIE.GWAS/GWAS.pipeline/Extracted.HQ.Pruned.Var.GT.Step.1"


mkdir -p "${OUTPUT_DIR}"


############################
# Environment check
############################

echo "Checking PLINK2..."

plink2 --version


echo ""
echo "=========================================="
echo "Extract LD-pruned variants"
echo "=========================================="


############################
# Extract variants chr1-22
############################

for chr in {1..22}; do

    echo ""
    echo "------------------------------------------"
    echo "Processing chromosome ${chr}"
    echo "------------------------------------------"


    CHR_OUT_DIR="${OUTPUT_DIR}/chr${chr}"

    mkdir -p "${CHR_OUT_DIR}"


    PFILE="${INPUT_DIR}/chr${chr}/chr${chr}.HQ_common"

    EXTRACT_FILE="${PRUNE_DIR}/chr${chr}/chr${chr}.LD_pruned.prune.in"

    OUT_PREFIX="${CHR_OUT_DIR}/chr${chr}.HQ_pruned"



    ############################
    # Sanity checks
    ############################

    if [[ ! -f "${PFILE}.pgen" ]]; then
        echo "ERROR: Missing ${PFILE}.pgen"
        exit 1
    fi


    if [[ ! -f "${PFILE}.pvar" ]]; then
        echo "ERROR: Missing ${PFILE}.pvar"
        exit 1
    fi


    if [[ ! -f "${PFILE}.psam" ]]; then
        echo "ERROR: Missing ${PFILE}.psam"
        exit 1
    fi


    if [[ ! -f "${EXTRACT_FILE}" ]]; then
        echo "ERROR: Missing ${EXTRACT_FILE}"
        exit 1
    fi



    echo "Input:"
    echo "  ${PFILE}"


    echo "Extract list:"
    echo "  ${EXTRACT_FILE}"


    echo "Output:"
    echo "  ${OUT_PREFIX}"



    ############################
    # PLINK2 extraction
    ############################

    plink2 \
        --pfile "${PFILE}" \
        --extract "${EXTRACT_FILE}" \
        --make-pgen \
        --out "${OUT_PREFIX}" \
        --threads 4



    ############################
    # Basic output check
    ############################

    if [[ ! -f "${OUT_PREFIX}.pgen" ]]; then
        echo "ERROR: Output pgen missing for chr${chr}"
        exit 1
    fi


    n_variants=$(grep -vc "^#" "${OUT_PREFIX}.pvar")


    echo "chr${chr}: ${n_variants} variants extracted"


done


echo ""
echo "=========================================="
echo "Extraction completed successfully"
echo "=========================================="

# Sanity check after finishing 22 chromosomes
cd /home/vscode/session_data/REGENIE.GWAS/GWAS.pipeline/Extracted.HQ.Pruned.Var.GT.Step.1

for chr in {1..22}; do
    echo -n "chr${chr}: "
    grep -vc "^#" chr${chr}/chr${chr}.HQ_pruned.pvar
done

