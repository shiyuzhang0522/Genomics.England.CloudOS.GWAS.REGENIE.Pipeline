#!/usr/bin/env bash

set -euo pipefail


################################################################################
# Merge chromosome-specific LD-pruned PLINK2 datasets
# for REGENIE Step 1
#
# Input:
#   Extracted.HQ.Pruned.Var.GT.Step.1/
#       chr1/chr1.HQ_pruned.*
#       ...
#       chr22/chr22.HQ_pruned.*
#
# Output:
#   REGENIE.Step1.HQ_pruned.GT.*
#
################################################################################


############################
# Paths
############################

BASE_DIR="/home/vscode/session_data/REGENIE.GWAS/GWAS.pipeline"


INPUT_DIR="${BASE_DIR}/Extracted.HQ.Pruned.Var.GT.Step.1"


OUTPUT_DIR="${BASE_DIR}/REGENIE.Step1.Genotype"


mkdir -p "${OUTPUT_DIR}"


MERGE_LIST="${OUTPUT_DIR}/chr_merge_list.txt"


OUT_PREFIX="${OUTPUT_DIR}/REGENIE.Step1.HQ_pruned"



############################
# Check PLINK2
############################

plink2 --version



############################
# Create merge list
############################

echo "Creating PLINK2 merge list..."

rm -f "${MERGE_LIST}"


for chr in {2..22}; do

    echo "${INPUT_DIR}/chr${chr}/chr${chr}.HQ_pruned" \
        >> "${MERGE_LIST}"

done


echo ""
echo "Merge list:"
cat "${MERGE_LIST}"



############################
# Check input datasets
############################

echo ""
echo "Checking input files..."


for chr in {1..22}; do

    PREFIX="${INPUT_DIR}/chr${chr}/chr${chr}.HQ_pruned"


    for ext in pgen pvar psam; do

        if [[ ! -f "${PREFIX}.${ext}" ]]; then
            echo "Missing ${PREFIX}.${ext}"
            exit 1
        fi

    done

done


echo "All chromosome files detected"



############################
# Merge
############################

echo ""
echo "Starting merge..."


plink2 \
    --pfile "${INPUT_DIR}/chr1/chr1.HQ_pruned" \
    --pmerge-list "${MERGE_LIST}" \
    --make-pgen \
    --out "${OUT_PREFIX}" \
    --threads 8



echo ""
echo "Merge completed"



############################
# Basic QC
############################

echo ""
echo "Final dataset summary:"


grep -vc "^#" "${OUT_PREFIX}.pvar"


echo "Samples:"
grep -v "^#" "${OUT_PREFIX}.psam" | wc -l


echo ""
echo "Output:"
echo "${OUT_PREFIX}.pgen"
echo "${OUT_PREFIX}.pvar"
echo "${OUT_PREFIX}.psam"