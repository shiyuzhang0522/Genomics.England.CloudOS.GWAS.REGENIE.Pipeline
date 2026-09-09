#!/usr/bin/env bash

set -euo pipefail


################################################################################
# Retrieve variant-level CADD PHRED scores for GEL WGS variants located within
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
# CADD source:
#   GEL chromosome-level VEP annotation tables
#
#   /home/vscode/session_data/filesystems/chromosome_level_VEP/
#       GEL.VEP.chr${CHR}.tsv
#
# Usage:
#   bash 02_retrieve_melanocyte_cCRE_variant_CADD.GEL.sh <CHR>
#
# Example:
#   bash 02_retrieve_melanocyte_cCRE_variant_CADD.GEL.sh 22
#
# Chromosomes:
#   1-22
#
# Coordinate conventions:
#   - GEL .pvar chromosome names:  1, 2, ..., 22
#   - GEL .pvar POS:               1-based
#   - cCRE BED chromosome names:   chr1, chr2, ..., chr22
#   - cCRE BED coordinates:        0-based, half-open
#
# Variant inclusion:
#   - ALL variants in dragen.pvar are retained at the annotation stage.
#   - No FILTER-based variant exclusion is performed here.
#   - No CADD threshold is applied here.
#   - CADD_PHRED >= 20 filtering will be performed later during mask building.
#
# CADD annotation:
#   - CADD_PHRED is obtained from the chromosome-level GEL VEP table.
#   - Matching is performed using the exact GEL variant ID.
#   - VEP tables are transcript-level and may contain multiple rows per variant.
#
# Transcript-level collapse:
#   - NA, ".", and empty CADD_PHRED values are ignored when collapsing.
#   - If a variant has no non-missing CADD_PHRED value:
#         final CADD_PHRED = NA
#   - If a variant has one distinct non-missing CADD_PHRED value:
#         that value is retained.
#   - If a variant has >1 distinct non-missing CADD_PHRED values:
#         all distinct values are recorded in a QC file and execution log,
#         and the maximum CADD_PHRED is assigned to the variant.
#
# Missing CADD:
#   - Every variant-cCRE pair is retained.
#   - Variants without a usable CADD_PHRED score are assigned CADD_PHRED = NA.
#
# Output:
#   One variant-level CADD annotation table per chromosome.
#
# Final output columns:
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
#   CADD_PHRED
#
# Required software:
#   bedtools
#   awk
#   sort
#   cut
#   wc
#   head
#
# Sorting:
#   Default sort threads: 16
#   Default sort memory:   8G
#
# These can be overridden at runtime, for example:
#
#   SORT_THREADS=8 SORT_MEM=4G \
#       bash 02_retrieve_melanocyte_cCRE_variant_CADD.GEL.sh 22
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


################################################################################
# Step 0a. Configure sorting resources
################################################################################

SORT_THREADS="${SORT_THREADS:-16}"
SORT_MEM="${SORT_MEM:-8G}"

if command -v nproc >/dev/null 2>&1; then
    N_CPU=$(nproc)

    if [[ "${SORT_THREADS}" -gt "${N_CPU}" ]]; then
        SORT_THREADS="${N_CPU}"
    fi
fi


echo "======================================================================"
echo "[INFO] Starting GEL melanocyte-cCRE CADD annotation"
echo "[INFO] Chromosome   : chr${CHR}"
echo "[INFO] Sort threads : ${SORT_THREADS}"
echo "[INFO] Sort memory  : ${SORT_MEM}"
echo "[INFO] Started      : $(date)"
echo "======================================================================"


################################################################################
# Step 1. Check required software
################################################################################

for tool in bedtools awk sort cut wc head; do

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

VEP="/home/vscode/session_data/filesystems/chromosome_level_VEP/GEL.VEP.chr${CHR}.tsv"

OUT_DIR="/home/vscode/session_data/CADD_raw/Melanocyte.cCRE.CADD.scores"

mkdir -p "${OUT_DIR}"


################################################################################
# Step 3. Check required input files
################################################################################

for file in \
    "${PVAR}" \
    "${CCRE_BED}" \
    "${VEP}"
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
echo "[INFO] VEP      : ${VEP}"
echo "[INFO] Output   : ${OUT_DIR}"


################################################################################
# Step 4. Define intermediate, QC, and final output files
################################################################################

CCRE_CHR="${OUT_DIR}/chr${CHR}.melanocyte.cCREs.bed"

VARIANT_BED="${OUT_DIR}/chr${CHR}.variants.bed"

VARIANT_CCRE="${OUT_DIR}/chr${CHR}.variants.in.melanocyte.cCREs.bed"

CCRE_VARIANT_IDS="${OUT_DIR}/chr${CHR}.cCRE.variant.ids.txt"

CADD_VALUES_RAW="${OUT_DIR}/chr${CHR}.CADD.nonmissing.transcript_values.tsv"

CADD_LOOKUP="${OUT_DIR}/chr${CHR}.CADD.variant.lookup.tsv"

MULTI_CADD="${OUT_DIR}/chr${CHR}.CADD.multiple_scores.tsv"

OUT="${OUT_DIR}/chr${CHR}.Melanocyte.cCRE.variant.CADD.tsv"


################################################################################
# Step 5. Extract chromosome-specific melanocyte cCREs
#
# Original cCRE BED:
#
#   chr1  start  end  cCRE_ID  agnostic_class  DNase_Z
#
# Convert:
#
#   chr1 -> 1
#
# so chromosome naming matches the GEL .pvar.
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
# GEL .pvar:
#
#   #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
#
# Example:
#
#   22  16050075  chr22:16050075:G:A  G  A  ...
#
# BED representation:
#
#   CHROM  POS-1  POS  ID  REF  ALT
#
# Important:
#   - ALL .pvar records are retained regardless of FILTER status.
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
# Step 7. Identify GEL variants overlapping melanocyte-specific cCREs
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
# A variant may occur more than once if it overlaps more than one cCRE.
# Such variant-cCRE pairs are intentionally retained.
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
# Step 8. Create unique cCRE variant-ID list
#
# CADD is allele-specific and will be matched using the exact GEL variant ID:
#
#   chr${CHR}:POS:REF:ALT
################################################################################

cut -f4 "${VARIANT_CCRE}" \
| LC_ALL=C sort \
    --parallel="${SORT_THREADS}" \
    -S "${SORT_MEM}" \
    -u \
> "${CCRE_VARIANT_IDS}"


N_UNIQUE_CCRE_VARIANTS=$(wc -l < "${CCRE_VARIANT_IDS}")


echo "[INFO] chr${CHR}: unique cCRE variants = ${N_UNIQUE_CCRE_VARIANTS}"


################################################################################
# Step 9. Identify ID and CADD_PHRED columns dynamically from the VEP header
#
# The current GEL VEP format has:
#
#   ID          = column 3
#   CADD_PHRED  = column 27
#
# However, column positions are detected by name rather than hard-coded.
################################################################################

read -r ID_COL CADD_COL < <(
    head -n 1 "${VEP}" \
    | awk -F'\t' '
    {
        id_col = 0
        cadd_col = 0

        for (i = 1; i <= NF; i++) {

            field = $i
            sub(/\r$/, "", field)

            if (field == "ID") {
                id_col = i
            }

            if (field == "CADD_PHRED") {
                cadd_col = i
            }
        }

        print id_col, cadd_col
    }
    '
)


if [[ "${ID_COL}" -eq 0 ]]; then
    echo "[ERROR] Could not identify ID column in VEP header." >&2
    exit 1
fi


if [[ "${CADD_COL}" -eq 0 ]]; then
    echo "[ERROR] Could not identify CADD_PHRED column in VEP header." >&2
    exit 1
fi


echo
echo "[INFO] VEP ID column         = ${ID_COL}"
echo "[INFO] VEP CADD_PHRED column = ${CADD_COL}"


################################################################################
# Step 10. Extract non-missing transcript-level CADD_PHRED values
#
# Only variants already identified within melanocyte-specific cCREs are
# extracted from the very large chromosome-level VEP table.
#
# Transcript-level VEP rows may repeat the same variant.
#
# Missing values:
#   NA
#   .
#   empty
#
# are ignored here.
#
# Non-missing CADD_PHRED values are required to be numeric.
#
# Output:
#
#   variant_ID    CADD_PHRED
################################################################################

echo
echo "[INFO] Scanning chromosome-level VEP table for cCRE CADD annotations..."


LC_ALL=C awk \
    -F'\t' \
    -v OFS="\t" \
    -v id_col="${ID_COL}" \
    -v cadd_col="${CADD_COL}" '
NR == FNR {
    keep[$1] = 1
    next
}

FNR == 1 {
    next
}

{
    id = $id_col
    val = $cadd_col

    sub(/\r$/, "", id)
    sub(/\r$/, "", val)

    if (!(id in keep)) {
        next
    }

    if (val == "NA" || val == "." || val == "") {
        next
    }

    if (val !~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/) {
        print "[ERROR] Non-numeric non-missing CADD_PHRED value:" \
              " ID=" id ", CADD_PHRED=" val > "/dev/stderr"
        bad = 1
        next
    }

    print id, val
}

END {
    if (bad) {
        exit 2
    }
}
' \
"${CCRE_VARIANT_IDS}" \
"${VEP}" \
> "${CADD_VALUES_RAW}"


N_CADD_TRANSCRIPT_ROWS=$(wc -l < "${CADD_VALUES_RAW}")


echo "[INFO] chr${CHR}: non-missing transcript-level CADD rows = ${N_CADD_TRANSCRIPT_ROWS}"


if [[ "${N_CADD_TRANSCRIPT_ROWS}" -eq 0 ]]; then
    echo "[ERROR] No non-missing CADD_PHRED values were found for cCRE variants." >&2
    echo "[ERROR] This may indicate an ID mismatch or unexpected VEP format." >&2
    exit 1
fi


################################################################################
# Step 11. Collapse transcript-level CADD annotations to one value per variant
#
# Rules:
#
#   1. Duplicate transcript rows carrying the same CADD_PHRED are collapsed.
#
#   2. If a variant has exactly one distinct non-missing score:
#          retain that score.
#
#   3. If a variant has >1 distinct non-missing score:
#          - record all distinct values;
#          - record their count;
#          - assign MAX(CADD_PHRED) as the final score.
#
# Numeric equality is used when defining distinct scores. Therefore values such
# as "2", "2.0", and "2.00" are treated as the same numerical score.
#
# CADD_LOOKUP:
#
#   ID    CADD_PHRED
#
# MULTI_CADD:
#
#   ID    N_distinct_CADD    CADD_values    MAX_CADD_PHRED
################################################################################

printf "ID\tCADD_PHRED\n" > "${CADD_LOOKUP}"

printf "ID\tN_distinct_CADD\tCADD_values\tMAX_CADD_PHRED\n" \
    > "${MULTI_CADD}"


LC_ALL=C sort \
    --parallel="${SORT_THREADS}" \
    -S "${SORT_MEM}" \
    -t $'\t' \
    -k1,1 \
    -k2,2n \
    "${CADD_VALUES_RAW}" \
| awk \
    -F'\t' \
    -v OFS="\t" \
    -v lookup="${CADD_LOOKUP}" \
    -v conflict="${MULTI_CADD}" '

function flush_variant(    i, values_string) {

    if (current_id == "") {
        return
    }

    print current_id, max_score_string >> lookup

    if (n_distinct > 1) {

        values_string = score_values[1]

        for (i = 2; i <= n_distinct; i++) {
            values_string = values_string "," score_values[i]
        }

        print \
            current_id, \
            n_distinct, \
            values_string, \
            max_score_string \
            >> conflict
    }
}


{
    id = $1
    score_string = $2
    score_numeric = $2 + 0

    if (id != current_id) {

        flush_variant()

        current_id = id
        n_distinct = 0
        have_previous_score = 0
        max_score_numeric = 0
        max_score_string = ""
    }

    if (!have_previous_score || score_numeric != previous_score_numeric) {

        n_distinct++
        score_values[n_distinct] = score_string

        previous_score_numeric = score_numeric
        have_previous_score = 1

        if (n_distinct == 1 || score_numeric > max_score_numeric) {
            max_score_numeric = score_numeric
            max_score_string = score_string
        }
    }
}


END {
    flush_variant()
}
'


N_CADD_LOOKUP=$(( $(wc -l < "${CADD_LOOKUP}") - 1 ))

N_MULTI_CADD=$(( $(wc -l < "${MULTI_CADD}") - 1 ))


echo
echo "[INFO] chr${CHR}: unique cCRE variants with CADD = ${N_CADD_LOOKUP}"
echo "[INFO] chr${CHR}: variants with >1 distinct CADD = ${N_MULTI_CADD}"


################################################################################
# Step 12. Record conflicting CADD annotations in the execution log
#
# Every variant with >1 distinct non-missing CADD_PHRED value is retained in:
#
#   chr${CHR}.CADD.multiple_scores.tsv
#
# and also printed to stdout/stderr so that the information is captured by the
# execution log.
################################################################################

if [[ "${N_MULTI_CADD}" -gt 0 ]]; then

    echo
    echo "[WARNING] ${N_MULTI_CADD} variants have multiple distinct CADD_PHRED values."
    echo "[WARNING] Maximum CADD_PHRED will be assigned for these variants."
    echo "[WARNING] Detailed QC file: ${MULTI_CADD}"

    awk -F'\t' '
    NR > 1 {
        print \
            "[CADD_CONFLICT]" \
            "\tID=" $1 \
            "\tN=" $2 \
            "\tVALUES=" $3 \
            "\tMAX=" $4
    }
    ' "${MULTI_CADD}"

else

    echo "[INFO] chr${CHR}: no variants have multiple distinct CADD_PHRED values."

fi


################################################################################
# Step 13. Verify that the collapsed lookup contains unique variant IDs
################################################################################

N_DUPLICATE_CADD_IDS=$(
    awk '
    NR > 1 {
        count[$1]++
    }

    END {
        n = 0

        for (id in count) {
            if (count[id] > 1) {
                n++
            }
        }

        print n + 0
    }
    ' "${CADD_LOOKUP}"
)


echo "[INFO] chr${CHR}: duplicate IDs in CADD lookup = ${N_DUPLICATE_CADD_IDS}"


if [[ "${N_DUPLICATE_CADD_IDS}" -ne 0 ]]; then
    echo "[ERROR] Duplicate variant IDs remain in the collapsed CADD lookup." >&2
    exit 1
fi


################################################################################
# Step 14. Annotate every variant-cCRE pair with CADD_PHRED
#
# All original variant-cCRE pairs are retained.
#
# If a variant is absent from CADD_LOOKUP:
#
#   CADD_PHRED = NA
################################################################################

printf \
"CHROM\tPOS\tID\tREF\tALT\tcCRE_ID\tcCRE_start\tcCRE_end\tagnostic_class\tDNase_Z\tCADD_PHRED\n" \
    > "${OUT}"


awk '
BEGIN {
    OFS="\t"
}

NR == FNR {

    if (FNR == 1) {
        next
    }

    cadd[$1] = $2
    next
}

{
    score = ($4 in cadd ? cadd[$4] : "NA")

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
        score
}
' \
"${CADD_LOOKUP}" \
"${VARIANT_CCRE}" \
>> "${OUT}"


################################################################################
# Step 15. Final annotation sanity checks
################################################################################

N_OUT=$(( $(wc -l < "${OUT}") - 1 ))


N_WITH_CADD=$(
    awk '
    NR > 1 && $11 != "NA" {
        n++
    }

    END {
        print n + 0
    }
    ' "${OUT}"
)


N_MISSING_CADD=$(
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
echo "[INFO] chr${CHR}: pairs with CADD score     = ${N_WITH_CADD}"
echo "[INFO] chr${CHR}: pairs missing CADD score  = ${N_MISSING_CADD}"


################################################################################
# Step 15a. Require exact preservation of variant-cCRE pairs
################################################################################

if [[ "${N_OUT}" -ne "${N_CCRE_VARIANTS}" ]]; then

    echo "[ERROR] Final output does not preserve the original" >&2
    echo "        variant-cCRE pair count." >&2

    echo "[ERROR] Input variant-cCRE pairs : ${N_CCRE_VARIANTS}" >&2
    echo "[ERROR] Final output rows         : ${N_OUT}" >&2

    exit 1
fi


################################################################################
# Step 15b. Verify CADD accounting
################################################################################

if [[ $((N_WITH_CADD + N_MISSING_CADD)) -ne "${N_OUT}" ]]; then

    echo "[ERROR] CADD annotation accounting failed." >&2
    echo "[ERROR] Output rows          : ${N_OUT}" >&2
    echo "[ERROR] With CADD            : ${N_WITH_CADD}" >&2
    echo "[ERROR] Missing CADD         : ${N_MISSING_CADD}" >&2

    exit 1
fi


################################################################################
# Step 15c. Check for unexpected duplicate final variant-cCRE rows
#
# Unique key:
#
#   CHROM + POS + ID + REF + ALT + cCRE_ID
#
# A variant legitimately overlapping two distinct cCREs is not considered a
# duplicate.
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
# Step 16. Completion summary
################################################################################

echo
echo "======================================================================"
echo "[PASS] GEL melanocyte-cCRE CADD annotation completed successfully"
echo
echo "Chromosome                           : chr${CHR}"
echo "Melanocyte-specific cCREs            : ${N_CCRE}"
echo "Total variants in dragen.pvar        : ${N_VARIANTS}"
echo "Variant-cCRE pairs                   : ${N_CCRE_VARIANTS}"
echo "Unique cCRE variants                 : ${N_UNIQUE_CCRE_VARIANTS}"
echo "Non-missing transcript CADD rows     : ${N_CADD_TRANSCRIPT_ROWS}"
echo "Unique cCRE variants with CADD       : ${N_CADD_LOOKUP}"
echo "Variants with >1 distinct CADD score : ${N_MULTI_CADD}"
echo "Pairs with CADD score                : ${N_WITH_CADD}"
echo "Pairs missing CADD score             : ${N_MISSING_CADD}"
echo "Duplicate variant-cCRE keys          : ${N_DUPLICATE_PAIRS}"
echo
echo "Final output:"
echo "${OUT}"
echo
echo "Multiple-CADD QC file:"
echo "${MULTI_CADD}"
echo
echo "Finished: $(date)"
echo "======================================================================"