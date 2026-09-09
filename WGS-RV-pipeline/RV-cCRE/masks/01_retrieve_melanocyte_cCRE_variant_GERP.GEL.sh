#!/usr/bin/env bash

set -euo pipefail


################################################################################
# Retrieve variant-level GERP scores for GEL WGS variants located within
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
# GERP source:
#   GRCh38 GERP conservation score BigWig
#
# Usage:
#   bash 01_retrieve_melanocyte_cCRE_variant_GERP.GEL.sh <CHR>
#
# Example:
#   bash 01_retrieve_melanocyte_cCRE_variant_GERP.GEL.sh 21
#
# Chromosomes:
#   1-22
#
# Coordinate conventions:
#   - GEL .pvar chromosome names:   1, 2, ..., 22
#   - GEL .pvar POS:                1-based
#   - cCRE BED chromosome names:    chr1, chr2, ..., chr22
#   - cCRE BED coordinates:         0-based, half-open
#   - GERP BigWig chromosome names: 1, 2, ..., 22
#   - GERP coordinates:             0-based, half-open
#
# Variant inclusion:
#   - ALL variants in dragen.pvar are retained at the annotation stage.
#   - No FILTER-based variant exclusion is performed here.
#   - No GERP threshold is applied here.
#   - GERP >= 2 filtering will be performed later during mask construction.
#
# Missing GERP:
#   - Every variant-cCRE pair is retained.
#   - If no GERP score is available at the variant position, GERP is set to NA.
#
# Output:
#   One variant-level GERP annotation table per chromosome.
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
#   GERP
#
# Required software:
#   bedtools
#   bigWigToBedGraph
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
    echo "[EXAMPLE] bash $0 21" >&2
    exit 1
fi

CHR="$1"

if ! [[ "${CHR}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "[ERROR] Invalid chromosome: ${CHR}" >&2
    echo "[ERROR] Expected an autosome number from 1 to 22." >&2
    exit 1
fi


echo "======================================================================"
echo "[INFO] Starting GEL melanocyte-cCRE GERP annotation"
echo "[INFO] Chromosome : chr${CHR}"
echo "[INFO] Started    : $(date)"
echo "======================================================================"


################################################################################
# Step 1. Check required software
################################################################################

for tool in bedtools bigWigToBedGraph awk sort wc; do

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

GERP_BW="/home/vscode/session_data/filesystems/gerp_conservation_scores.homo_sapiens.GRCh38.bw"

OUT_DIR="/home/vscode/session_data/GERP_raw/Melanocyte.cCRE.GERP.scores"

mkdir -p "${OUT_DIR}"


################################################################################
# Step 3. Check required input files
################################################################################

for file in \
    "${PVAR}" \
    "${CCRE_BED}" \
    "${GERP_BW}"
do

    if [[ ! -s "${file}" ]]; then
        echo "[ERROR] Missing or empty input file:" >&2
        echo "        ${file}" >&2
        exit 1
    fi

done

echo
echo "[INFO] PVAR     : ${PVAR}"
echo "[INFO] cCRE BED : ${CCRE_BED}"
echo "[INFO] GERP BW  : ${GERP_BW}"
echo "[INFO] Output   : ${OUT_DIR}"


################################################################################
# Step 4. Define intermediate and final files
################################################################################

CCRE_CHR="${OUT_DIR}/chr${CHR}.melanocyte.cCREs.bed"

VARIANT_BED="${OUT_DIR}/chr${CHR}.variants.bed"

VARIANT_CCRE="${OUT_DIR}/chr${CHR}.variants.in.melanocyte.cCREs.bed"

GERP_ANNOTATED="${OUT_DIR}/chr${CHR}.variants.in.melanocyte.cCREs.GERP.raw.bed"

OUT="${OUT_DIR}/chr${CHR}.Melanocyte.cCRE.variant.GERP.tsv"


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
# This makes cCRE chromosome naming consistent with both:
#   - GEL dragen.pvar
#   - GERP BigWig
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
#   21  5010006  chr21:5010006:A:C  A  C  0  LowGTR  ...
#
# Convert:
#
#   CHROM  POS  ID  REF  ALT
#
# into:
#
#   CHROM  POS-1  POS  ID  REF  ALT
#
# Example:
#
#   21  5010005  5010006  chr21:5010006:A:C  A  C
#
# Important:
#   - ALL records are retained regardless of FILTER status.
#   - Header lines beginning with "#" are excluded.
################################################################################

awk '
BEGIN {
    OFS="\t"
}

!/^#/ {

    if (NF < 5) {
        print "[ERROR] Malformed .pvar record at line " NR > "/dev/stderr"
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
# A variant may appear more than once if it overlaps more than one cCRE.
# These are intentionally retained as distinct variant-cCRE pairs.
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
# Step 8. Retrieve positional GERP scores
#
# GERP BigWig chromosome naming:
#
#   1, 2, ..., 22
#
# bigWigToBedGraph output:
#
#   CHROM  start  end  GERP_score
#
# bedtools intersect -loj performs a left outer join:
#
#   - every variant-cCRE pair in VARIANT_CCRE is retained;
#   - overlapping GERP records are appended;
#   - if no GERP interval overlaps the variant, missing values are emitted
#     by bedtools and subsequently converted to GERP = NA.
#
# Output of this step:
#
#   columns 1-12  : variant-cCRE pair
#   columns 13-16 : GERP BedGraph annotation
#
#   column 16     : GERP score
################################################################################

bigWigToBedGraph \
    -chrom="${CHR}" \
    "${GERP_BW}" \
    stdout \
| bedtools intersect \
    -a "${VARIANT_CCRE}" \
    -b stdin \
    -wa \
    -wb \
    -loj \
    > "${GERP_ANNOTATED}"


N_GERP_RAW=$(wc -l < "${GERP_ANNOTATED}")


echo "[INFO] chr${CHR}: rows after GERP left join = ${N_GERP_RAW}"


################################################################################
# Step 9. Write final variant-level GERP annotation table
#
# Every original variant-cCRE pair must be retained.
#
# Missing GERP values are written as:
#
#   NA
################################################################################

printf \
"CHROM\tPOS\tID\tREF\tALT\tcCRE_ID\tcCRE_start\tcCRE_end\tagnostic_class\tDNase_Z\tGERP\n" \
    > "${OUT}"


awk '
BEGIN {
    OFS="\t"
}

{
    gerp = $16

    if (gerp == "." || gerp == "") {
        gerp = "NA"
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
        gerp
}
' "${GERP_ANNOTATED}" >> "${OUT}"

################################################################################
# Step 10. Final sanity checks
################################################################################

N_OUT=$(( $(wc -l < "${OUT}") - 1 ))


N_WITH_GERP=$(
    awk '
    NR > 1 && $11 != "NA" {
        n++
    }
    END {
        print n + 0
    }
    ' "${OUT}"
)


N_MISSING_GERP=$(
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
echo "[INFO] chr${CHR}: pairs with GERP score     = ${N_WITH_GERP}"
echo "[INFO] chr${CHR}: pairs missing GERP score  = ${N_MISSING_GERP}"


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
# Step 10b. Verify GERP accounting
################################################################################

if [[ $((N_WITH_GERP + N_MISSING_GERP)) -ne "${N_OUT}" ]]; then

    echo "[ERROR] GERP annotation accounting failed." >&2
    echo "[ERROR] Output rows          : ${N_OUT}" >&2
    echo "[ERROR] With GERP            : ${N_WITH_GERP}" >&2
    echo "[ERROR] Missing GERP         : ${N_MISSING_GERP}" >&2

    exit 1
fi


################################################################################
# Step 10c. Check for unexpected duplicate final rows
#
# The unique key here is:
#
#   CHROM + POS + ID + REF + ALT + cCRE_ID
#
# A variant legitimately occurring in two distinct cCREs is NOT considered
# a duplicate.
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
echo "[PASS] GEL melanocyte-cCRE GERP annotation completed successfully"
echo
echo "Chromosome                     : chr${CHR}"
echo "Melanocyte-specific cCREs      : ${N_CCRE}"
echo "Total variants in dragen.pvar  : ${N_VARIANTS}"
echo "Variant-cCRE pairs             : ${N_CCRE_VARIANTS}"
echo "Pairs with GERP score          : ${N_WITH_GERP}"
echo "Pairs missing GERP score       : ${N_MISSING_GERP}"
echo "Duplicate variant-cCRE keys    : ${N_DUPLICATE_PAIRS}"
echo
echo "Output:"
echo "${OUT}"
echo
echo "Finished: $(date)"
echo "======================================================================"