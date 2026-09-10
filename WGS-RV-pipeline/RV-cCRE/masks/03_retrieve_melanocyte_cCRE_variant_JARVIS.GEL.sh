#!/usr/bin/env bash

set -euo pipefail


################################################################################
# Retrieve variant-level JARVIS scores for GEL WGS variants located within
# epidermal melanocyte-specific cCREs
#
# Dataset:
#   Genomics England (GEL) AggV3 WGS
#
# Variant source:
#   Chromosome-specific DRAGEN PLINK2 .pvar files
#
#   /home/vscode/session_data/filesystems/chrom-msvcf/
#       chrom-${CHR}/postproc-pgen/dragen.pvar
#
# cCRE source:
#   ENCODE4 epidermal melanocyte-specific cCREs
#
# JARVIS source:
#   GRCh38 JARVIS BigWig
#
# Usage:
#   bash 03_retrieve_melanocyte_cCRE_variant_JARVIS.GEL.sh <CHR>
#
# Example:
#   bash 03_retrieve_melanocyte_cCRE_variant_JARVIS.GEL.sh 22
#
# Chromosomes:
#   1-22
#
# Coordinate conventions:
#   - GEL .pvar chromosome names:    1, 2, ..., 22
#   - GEL .pvar POS:                 1-based
#   - cCRE BED chromosome names:     chr1, chr2, ..., chr22
#   - cCRE BED coordinates:          0-based, half-open
#   - JARVIS BigWig chromosomes:     chr1, chr2, ..., chr22
#   - JARVIS BigWig coordinates:     0-based, half-open
#
# Variant inclusion:
#   - ALL variants in dragen.pvar are retained at the annotation stage.
#   - No FILTER-based variant exclusion is performed here.
#   - No JARVIS threshold is applied here.
#
# Missing JARVIS:
#   - Every variant-cCRE pair is retained.
#   - If no JARVIS score is available at the variant position,
#     JARVIS is set to NA.
#
# Output:
#   One variant-level JARVIS annotation table per chromosome.
#
# Output columns:
#   CHROM
#   POS
#   ID
#   REF
#   ALT
#   cCRE_ID
#   cCRE_start
#   cCRE_end
#   agnostic_class
#   DNase_Z
#   JARVIS
#
# Required software:
#   bedtools
#   bigWigToBedGraph
#   awk
#   wc
#
# Recommended environment:
#   GEL_RVAT_env
#
# Author: Shelley
# Date:   2026-09-09
################################################################################


################################################################################
# Step 0. Parse and validate chromosome
################################################################################

if [[ $# -ne 1 ]]; then
    echo "[ERROR] Exactly one chromosome argument is required." >&2
    echo "[USAGE] bash $0 <CHR>" >&2
    echo "[EXAMPLE] bash $0 22" >&2
    exit 1
fi

CHR="$1"

if ! [[ "${CHR}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "[ERROR] Invalid chromosome: ${CHR}" >&2
    echo "[ERROR] Expected an autosome number from 1 to 22." >&2
    exit 1
fi


echo "======================================================================"
echo "[INFO] Starting GEL melanocyte-cCRE JARVIS annotation"
echo "[INFO] Chromosome : chr${CHR}"
echo "[INFO] Started    : $(date)"
echo "======================================================================"


################################################################################
# Step 1. Check required software
################################################################################

for tool in bedtools bigWigToBedGraph awk wc; do

    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "[ERROR] Required tool not found: ${tool}" >&2
        exit 1
    fi

    echo "[INFO] ${tool}: $(command -v "${tool}")"

done


################################################################################
# Step 2. Define input and output paths
################################################################################

PVAR="/home/vscode/session_data/filesystems/chrom-msvcf/chrom-${CHR}/postproc-pgen/dragen.pvar"

CCRE_BED="/home/vscode/session_data/filesystems/ENCODE4.Epidermal.Melanocyte.cCREs.chr1-22.bed"

JARVIS_BW="/home/vscode/session_data/filesystems/jarvis.bw"

OUT_DIR="/home/vscode/session_data/JARVIS_raw/Melanocyte.cCRE.JARVIS.scores"

mkdir -p "${OUT_DIR}"


################################################################################
# Step 3. Check required input files
################################################################################

for file in \
    "${PVAR}" \
    "${CCRE_BED}" \
    "${JARVIS_BW}"
do

    if [[ ! -s "${file}" ]]; then
        echo "[ERROR] Missing or empty input file:" >&2
        echo "        ${file}" >&2
        exit 1
    fi

done


echo
echo "[INFO] PVAR      : ${PVAR}"
echo "[INFO] cCRE BED  : ${CCRE_BED}"
echo "[INFO] JARVIS BW : ${JARVIS_BW}"
echo "[INFO] Output    : ${OUT_DIR}"


################################################################################
# Step 4. Define intermediate and final files
################################################################################

CCRE_CHR="${OUT_DIR}/chr${CHR}.melanocyte.cCREs.bed"

VARIANT_BED="${OUT_DIR}/chr${CHR}.variants.bed"

VARIANT_CCRE="${OUT_DIR}/chr${CHR}.variants.in.melanocyte.cCREs.bed"

JARVIS_ANNOTATED="${OUT_DIR}/chr${CHR}.variants.in.melanocyte.cCREs.JARVIS.raw.bed"

OUT="${OUT_DIR}/chr${CHR}.Melanocyte.cCRE.variant.JARVIS.tsv"


################################################################################
# Step 5. Extract chromosome-specific melanocyte cCREs
#
# Input cCRE BED:
#
#   chr1  start  end  cCRE_ID  agnostic_class  DNase_Z
#
# Convert chromosome naming:
#
#   chr1 -> 1
#
# to match the GEL .pvar working representation.
################################################################################

awk -v chr="chr${CHR}" '
BEGIN {
    OFS="\t"
}

$1 == chr {
    sub(/^chr/, "", $1)
    print
}
' "${CCRE_BED}" > "${CCRE_CHR}"


N_CCRE=$(wc -l < "${CCRE_CHR}")


echo
echo "[INFO] chr${CHR}: melanocyte-specific cCREs = ${N_CCRE}"


if [[ "${N_CCRE}" -eq 0 ]]; then
    echo "[ERROR] No melanocyte-specific cCREs found for chr${CHR}." >&2
    exit 1
fi


################################################################################
# Step 6. Convert GEL dragen.pvar variants to 1-bp BED intervals
#
# GEL .pvar columns:
#
#   #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
#
# Example:
#
#   22  16559110  chr22:16559110:A:G  A  G  ...
#
# Convert to:
#
#   CHROM  POS-1  POS  ID  REF  ALT
#
# Important:
#   - ALL records are retained regardless of FILTER status.
################################################################################

awk -v expected_chr="${CHR}" '
BEGIN {
    OFS="\t"
}

!/^#/ {

    if (NF < 5) {
        print "[ERROR] Malformed .pvar record at line " NR > "/dev/stderr"
        exit 1
    }

    if ($1 != expected_chr) {
        print "[ERROR] Unexpected chromosome in .pvar at line " NR \
              ": observed=" $1 ", expected=" expected_chr > "/dev/stderr"
        exit 1
    }

    if ($2 !~ /^[0-9]+$/) {
        print "[ERROR] Non-numeric POS at line " NR ": " $2 > "/dev/stderr"
        exit 1
    }

    print \
        $1,      \
        $2 - 1,  \
        $2,      \
        $3,      \
        $4,      \
        $5
}
' "${PVAR}" > "${VARIANT_BED}"


N_VARIANTS=$(wc -l < "${VARIANT_BED}")


echo "[INFO] chr${CHR}: total variants in dragen.pvar = ${N_VARIANTS}"


if [[ "${N_VARIANTS}" -eq 0 ]]; then
    echo "[ERROR] No variants were extracted from:" >&2
    echo "        ${PVAR}" >&2
    exit 1
fi


################################################################################
# Step 7. Retain variants overlapping melanocyte-specific cCREs
#
# VARIANT_CCRE columns:
#
#    1  variant_CHROM
#    2  variant_start
#    3  variant_end
#    4  variant_ID
#    5  REF
#    6  ALT
#
#    7  cCRE_CHROM
#    8  cCRE_start
#    9  cCRE_end
#   10  cCRE_ID
#   11  agnostic_class
#   12  DNase_Z
#
# A variant may legitimately occur more than once if it overlaps more than
# one cCRE. Such variant-cCRE pairs are intentionally retained.
################################################################################

bedtools intersect \
    -a "${VARIANT_BED}" \
    -b "${CCRE_CHR}" \
    -wa \
    -wb \
    > "${VARIANT_CCRE}"


N_CCRE_VARIANTS=$(wc -l < "${VARIANT_CCRE}")


echo "[INFO] chr${CHR}: variant-cCRE pairs = ${N_CCRE_VARIANTS}"


if [[ "${N_CCRE_VARIANTS}" -eq 0 ]]; then
    echo "[ERROR] No GEL variants overlap melanocyte-specific cCREs on chr${CHR}." >&2
    exit 1
fi


################################################################################
# Step 8. Retrieve positional JARVIS scores
#
# JARVIS BigWig chromosome naming:
#
#   chr1, chr2, ..., chr22
#
# bigWigToBedGraph output:
#
#   chr22  start  end  JARVIS_score
#
# Convert:
#
#   chr22 -> 22
#
# before intersecting with VARIANT_CCRE.
#
# bedtools intersect -loj performs a left outer join:
#
#   - every variant-cCRE pair is retained;
#   - JARVIS information is appended when available;
#   - missing JARVIS annotations remain present and are later converted to NA.
#
# Final raw columns:
#
#   1-12   variant-cCRE
#   13     JARVIS_CHROM
#   14     JARVIS_start
#   15     JARVIS_end
#   16     JARVIS_score
################################################################################

bigWigToBedGraph \
    -chrom="chr${CHR}" \
    "${JARVIS_BW}" \
    stdout \
| awk '
BEGIN {
    OFS="\t"
}

{
    sub(/^chr/, "", $1)
    print
}
' \
| bedtools intersect \
    -a "${VARIANT_CCRE}" \
    -b stdin \
    -wa \
    -wb \
    -loj \
    > "${JARVIS_ANNOTATED}"


N_JARVIS_RAW=$(wc -l < "${JARVIS_ANNOTATED}")


echo "[INFO] chr${CHR}: rows after JARVIS left join = ${N_JARVIS_RAW}"


################################################################################
# Step 9. Write final variant-level JARVIS annotation table
#
# Every original variant-cCRE pair must be retained.
#
# Missing JARVIS values are written as:
#
#   NA
################################################################################

printf \
"CHROM\tPOS\tID\tREF\tALT\tcCRE_ID\tcCRE_start\tcCRE_end\tagnostic_class\tDNase_Z\tJARVIS\n" \
    > "${OUT}"


awk '
BEGIN {
    OFS="\t"
}

{
    jarvis = $16

    if (jarvis == "." || jarvis == "") {
        jarvis = "NA"
    }

    print \
        $1,  \
        $3,  \
        $4,  \
        $5,  \
        $6,  \
        $10, \
        $8,  \
        $9,  \
        $11, \
        $12, \
        jarvis
}
' "${JARVIS_ANNOTATED}" >> "${OUT}"


################################################################################
# Step 10. Final sanity checks
################################################################################

N_OUT=$(( $(wc -l < "${OUT}") - 1 ))


N_WITH_JARVIS=$(
    awk '
    NR > 1 && $11 != "NA" {
        n++
    }

    END {
        print n + 0
    }
    ' "${OUT}"
)


N_MISSING_JARVIS=$(
    awk '
    NR > 1 && $11 == "NA" {
        n++
    }

    END {
        print n + 0
    }
    ' "${OUT}"
)


echo
echo "[INFO] chr${CHR}: output variant-cCRE pairs = ${N_OUT}"
echo "[INFO] chr${CHR}: pairs with JARVIS score   = ${N_WITH_JARVIS}"
echo "[INFO] chr${CHR}: pairs missing JARVIS      = ${N_MISSING_JARVIS}"


################################################################################
# Step 10a. Require exact preservation of variant-cCRE pairs
################################################################################

if [[ "${N_OUT}" -ne "${N_CCRE_VARIANTS}" ]]; then

    echo "[ERROR] Final output does not preserve the original" >&2
    echo "        variant-cCRE pair count." >&2

    echo "[ERROR] Input variant-cCRE pairs : ${N_CCRE_VARIANTS}" >&2
    echo "[ERROR] Final output rows         : ${N_OUT}" >&2

    exit 1
fi


################################################################################
# Step 10b. Verify JARVIS accounting
################################################################################

if [[ $((N_WITH_JARVIS + N_MISSING_JARVIS)) -ne "${N_OUT}" ]]; then

    echo "[ERROR] JARVIS annotation accounting failed." >&2
    echo "[ERROR] Output rows       : ${N_OUT}" >&2
    echo "[ERROR] With JARVIS       : ${N_WITH_JARVIS}" >&2
    echo "[ERROR] Missing JARVIS    : ${N_MISSING_JARVIS}" >&2

    exit 1
fi


################################################################################
# Step 10c. Check for unexpected duplicate final rows
#
# Unique key:
#
#   CHROM + POS + ID + REF + ALT + cCRE_ID
#
# A variant overlapping two distinct cCREs is not considered a duplicate.
################################################################################

N_DUPLICATE_PAIRS=$(
    awk '
    NR > 1 {
        key = $1 FS $2 FS $3 FS $4 FS $5 FS $6
        count[key]++
    }

    END {
        n = 0

        for (key in count) {
            if (count[key] > 1) {
                n++
            }
        }

        print n + 0
    }
    ' "${OUT}"
)


echo "[INFO] chr${CHR}: duplicated variant-cCRE keys = ${N_DUPLICATE_PAIRS}"


if [[ "${N_DUPLICATE_PAIRS}" -ne 0 ]]; then

    echo "[ERROR] Unexpected duplicated variant-cCRE annotations detected." >&2
    echo "[ERROR] Number of duplicated keys: ${N_DUPLICATE_PAIRS}" >&2

    exit 1
fi


################################################################################
# Step 11. Completion summary
################################################################################

echo
echo "======================================================================"
echo "[PASS] GEL melanocyte-cCRE JARVIS annotation completed successfully"
echo
echo "Chromosome                     : chr${CHR}"
echo "Melanocyte-specific cCREs      : ${N_CCRE}"
echo "Total variants in dragen.pvar  : ${N_VARIANTS}"
echo "Variant-cCRE pairs             : ${N_CCRE_VARIANTS}"
echo "Pairs with JARVIS score        : ${N_WITH_JARVIS}"
echo "Pairs missing JARVIS score     : ${N_MISSING_JARVIS}"
echo "Duplicate variant-cCRE keys    : ${N_DUPLICATE_PAIRS}"
echo
echo "Final output:"
echo "${OUT}"
echo
echo "Finished: $(date)"
echo "======================================================================"