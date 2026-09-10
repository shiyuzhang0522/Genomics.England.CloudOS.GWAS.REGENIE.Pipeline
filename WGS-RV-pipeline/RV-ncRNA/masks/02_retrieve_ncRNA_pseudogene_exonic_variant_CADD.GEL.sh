#!/usr/bin/env bash

################################################################################
# Retrieve variant-level CADD PHRED scores for Genomics England (GEL) WGS
# variants located within annotated exons of autosomal ncRNA/pseudogene genes.
#
# Dataset:
#   Genomics England (GEL) AggV3 WGS
#
# Variant source:
#   Chromosome-specific DRAGEN PVAR:
#
#     /home/vscode/session_data/filesystems/
#       chrom-msvcf/chrom-${CHR}/postproc-pgen/dragen.pvar
#
# ncRNA/pseudogene exon annotation:
#   Ensembl release 116, GRCh38
#
#     Ensembl116.GRCh38.autosomal.ncRNA_pseudogene.
#       gene_union_exons.merged.tsv
#
#   Columns:
#       1. ensembl_gene_id
#       2. external_gene_name
#       3. gene_biotype
#       4. chromosome_name
#       5. exon_start
#       6. exon_end
#       7. strand
#
# CADD source:
#   GEL chromosome-level VEP tables:
#
#     /home/vscode/session_data/filesystems/chromosome_level_VEP/
#       GEL.VEP.chr${CHR}.tsv
#
#   Required VEP columns are detected dynamically by header name:
#
#       ID
#       CADD_PHRED
#
# Coordinate conventions:
#   - GEL PVAR chromosome names:    1, 2, ..., 22
#   - Ensembl chromosome names:     1, 2, ..., 22
#   - Ensembl exon coordinates:     1-based, inclusive
#   - BED coordinates:              0-based, half-open
#   - GEL PVAR POS:                 1-based
#
# Variant-to-gene mapping:
#   - Each WGS variant is represented by its PVAR POS as:
#
#         [POS-1, POS)
#
#   - Exonic intervals are gene-level unions of all annotated exons across
#     transcripts for each ncRNA/pseudogene Ensembl gene.
#
#   - A variant may legitimately map to more than one gene when exon intervals
#     belonging to different genes overlap.
#
# CADD annotation:
#   - CADD is allele-specific.
#   - Exact GEL variant IDs are used for matching:
#
#         chrCHR:POS:REF:ALT
#
#   - GEL VEP is transcript-level, so the same variant may occur on multiple
#     transcript rows.
#
#   - Missing CADD values ("NA", ".", empty) are ignored when constructing the
#     variant-level lookup.
#
#   - For each exact variant ID:
#
#       no non-missing CADD values
#           -> final CADD_PHRED = NA
#
#       one distinct non-missing CADD value
#           -> use that value
#
#       >1 distinct non-missing CADD values
#           -> record all distinct values in a QC table
#           -> use the maximum CADD_PHRED
#
#   - All original variant-gene pairs are retained.
#   - No CADD threshold is applied here.
#   - CADD_PHRED >= 20 will be applied later during mask construction.
#
# FILTER policy:
#   - All variants present in dragen.pvar are considered irrespective of
#     FILTER status, consistent with the GEL cCRE annotation workflow.
#
# Final output:
#
#   CHROM
#   POS
#   ID
#   REF
#   ALT
#   ensembl_gene_id
#   external_gene_name
#   gene_biotype
#   exon_start
#   exon_end
#   strand
#   CADD_PHRED
#
# Usage:
#
#   bash 02_retrieve_ncRNA_pseudogene_exonic_variant_CADD.GEL.sh 22
#
# Required software:
#   awk
#   bedtools
#   sort
#   cut
#   wc
#
# Author: Shelley
# Date:   2026-09-10
################################################################################


set -euo pipefail


################################################################################
# Step 1. Parse chromosome argument
################################################################################

if [[ $# -ne 1 ]]; then
    echo "[ERROR] Usage: $0 <chromosome 1-22>" >&2
    exit 1
fi

CHR="$1"

if ! [[ "${CHR}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "[ERROR] Chromosome must be an integer from 1 to 22." >&2
    exit 1
fi


################################################################################
# Step 2. Check required software
################################################################################

for tool in \
    awk \
    bedtools \
    sort \
    cut \
    wc
do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "[ERROR] Required tool not found: ${tool}" >&2
        exit 1
    fi
done


################################################################################
# Step 3. Sorting resources
#
# These can be overridden before running, for example:
#
#   SORT_THREADS=8 SORT_MEM=6G bash script.sh 22
################################################################################

SORT_THREADS="${SORT_THREADS:-8}"
SORT_MEM="${SORT_MEM:-8G}"

AVAILABLE_CPUS=$(nproc 2>/dev/null || echo 1)

if [[ "${SORT_THREADS}" -gt "${AVAILABLE_CPUS}" ]]; then
    SORT_THREADS="${AVAILABLE_CPUS}"
fi


################################################################################
# Step 4. Define input and output paths
################################################################################

PVAR="/home/vscode/session_data/filesystems/chrom-msvcf/chrom-${CHR}/postproc-pgen/dragen.pvar"

NCRNA_EXONS="/home/vscode/session_data/filesystems/Ensembl116.GRCh38.autosomal.ncRNA_pseudogene.gene_union_exons.merged.tsv"

VEP="/home/vscode/session_data/filesystems/chromosome_level_VEP/GEL.VEP.chr${CHR}.tsv"

OUT_DIR="/home/vscode/session_data/CADD_ncRNA"

LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${OUT_DIR}"
mkdir -p "${LOG_DIR}"


################################################################################
# Step 5. Define intermediate, QC and final files
################################################################################

EXON_BED="${OUT_DIR}/chr${CHR}.ncRNA_pseudogene.exons.bed"

VARIANT_BED="${OUT_DIR}/chr${CHR}.variants.bed"

VARIANT_GENE="${OUT_DIR}/chr${CHR}.variants.in.ncRNA_pseudogene.exons.bed"

EXONIC_VARIANT_IDS="${OUT_DIR}/chr${CHR}.ncRNA_exonic.variant.ids.txt"

CADD_TRANSCRIPT_VALUES="${OUT_DIR}/chr${CHR}.CADD.nonmissing.transcript_values.tsv"

CADD_SORTED_VALUES="${OUT_DIR}/chr${CHR}.CADD.nonmissing.transcript_values.sorted.tsv"

CADD_LOOKUP="${OUT_DIR}/chr${CHR}.CADD.variant.lookup.tsv"

MULTI_CADD="${OUT_DIR}/chr${CHR}.CADD.multiple_scores.tsv"

OUT="${OUT_DIR}/chr${CHR}.ncRNA_pseudogene.exonic_variants.CADD.tsv"

SORT_TMP="${OUT_DIR}/tmp_sort_chr${CHR}"

mkdir -p "${SORT_TMP}"

trap 'rm -rf "${SORT_TMP}"' EXIT


################################################################################
# Step 6. Initial report
################################################################################

echo
echo "======================================================================"
echo "[INFO] GEL ncRNA/pseudogene exon CADD annotation"
echo "[INFO] Chromosome   : chr${CHR}"
echo "[INFO] Started      : $(date)"
echo "======================================================================"

echo
echo "[INFO] Input PVAR:"
echo "       ${PVAR}"

echo "[INFO] ncRNA/pseudogene exon annotation:"
echo "       ${NCRNA_EXONS}"

echo "[INFO] GEL chromosome-level VEP:"
echo "       ${VEP}"

echo "[INFO] Output directory:"
echo "       ${OUT_DIR}"

echo
echo "[INFO] Sort threads : ${SORT_THREADS}"
echo "[INFO] Sort memory  : ${SORT_MEM}"


################################################################################
# Step 7. Check required input files
################################################################################

for file in \
    "${PVAR}" \
    "${NCRNA_EXONS}" \
    "${VEP}"
do
    if [[ ! -s "${file}" ]]; then
        echo "[ERROR] Missing or empty input file:" >&2
        echo "        ${file}" >&2
        exit 1
    fi
done

echo
echo "[PASS] Required input files found."


################################################################################
# Step 8. Validate PVAR header
################################################################################

PVAR_HEADER=$(
    awk -F '\t' '
    $1 == "#CHROM" {
        print
        exit
    }
    ' "${PVAR}"
)

if [[ -z "${PVAR_HEADER}" ]]; then
    echo "[ERROR] #CHROM header not found in PVAR:" >&2
    echo "        ${PVAR}" >&2
    exit 1
fi

PVAR_NCOL=$(
    awk -F '\t' '
    $1 == "#CHROM" {
        print NF
        exit
    }
    ' "${PVAR}"
)

if [[ "${PVAR_NCOL}" -lt 5 ]]; then
    echo "[ERROR] PVAR contains fewer than five required columns." >&2
    exit 1
fi

echo "[PASS] PVAR header detected."


################################################################################
# Step 9. Detect VEP ID and CADD_PHRED columns
#
# We deliberately detect these columns by exact header name rather than
# hard-coding column numbers.
################################################################################

ID_COL=$(
    awk -F '\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "ID") {
                print i
                exit
            }
        }
    }
    ' "${VEP}"
)

CADD_COL=$(
    awk -F '\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "CADD_PHRED") {
                print i
                exit
            }
        }
    }
    ' "${VEP}"
)

if [[ -z "${ID_COL}" ]]; then
    echo "[ERROR] ID column not found in VEP header." >&2
    exit 1
fi

if [[ -z "${CADD_COL}" ]]; then
    echo "[ERROR] CADD_PHRED column not found in VEP header." >&2
    exit 1
fi

echo
echo "[PASS] Required VEP columns detected."
echo "[INFO] VEP ID column         : ${ID_COL}"
echo "[INFO] VEP CADD_PHRED column : ${CADD_COL}"


################################################################################
# Step 10. Extract chromosome-specific ncRNA/pseudogene merged exons
#
# Input TSV:
#
#   1  ensembl_gene_id
#   2  external_gene_name
#   3  gene_biotype
#   4  chromosome_name
#   5  exon_start             [1-based inclusive]
#   6  exon_end               [1-based inclusive]
#   7  strand
#
# BED output:
#
#   1  CHROM
#   2  exon_start - 1
#   3  exon_end
#   4  ensembl_gene_id
#   5  external_gene_name
#   6  gene_biotype
#   7  strand
#
# Empty gene names are replaced with NA so that the seven-field BED structure
# is preserved.
################################################################################

awk \
    -F '\t' \
    -v OFS='\t' \
    -v chr="${CHR}" '

NR == 1 {
    next
}

$4 == chr {

    gene_name = $2

    if (gene_name == "" || gene_name == ".") {
        gene_name = "NA"
    }

    if ($1 == "" || $3 == "" || $5 !~ /^[0-9]+$/ || $6 !~ /^[0-9]+$/ || $5 > $6) {
        print "[ERROR] Invalid ncRNA exon record:", $0 > "/dev/stderr"
        exit 1
    }

    print \
        $4,        \
        $5 - 1,    \
        $6,        \
        $1,        \
        gene_name, \
        $3,        \
        $7
}

' "${NCRNA_EXONS}" > "${EXON_BED}"

N_EXONS=$(wc -l < "${EXON_BED}")

echo
echo "[INFO] chr${CHR}: merged ncRNA/pseudogene exon intervals = ${N_EXONS}"

if [[ "${N_EXONS}" -eq 0 ]]; then
    echo "[ERROR] No ncRNA/pseudogene exon intervals found on chr${CHR}." >&2
    exit 1
fi


################################################################################
# Step 11. Convert all GEL PVAR variants to 1-bp BED intervals
#
# No FILTER-based exclusion is applied.
#
# BED:
#
#   1  CHROM
#   2  POS - 1
#   3  POS
#   4  ID
#   5  REF
#   6  ALT
################################################################################

awk \
    -F '\t' \
    -v OFS='\t' \
    -v chr="${CHR}" '

BEGIN {
    header_seen = 0
}

$1 == "#CHROM" {
    header_seen = 1
    next
}

/^##/ {
    next
}

header_seen {

    if ($1 != chr) {
        print "[ERROR] Unexpected chromosome in chr" chr " PVAR: " $1 > "/dev/stderr"
        exit 1
    }

    if ($2 !~ /^[0-9]+$/ || $2 < 1) {
        print "[ERROR] Invalid POS in PVAR: " $2 > "/dev/stderr"
        exit 1
    }

    if ($3 == "" || $3 == "." || $4 == "" || $4 == "." || $5 == "" || $5 == ".") {
        print "[ERROR] Missing ID/REF/ALT in PVAR record: " $0 > "/dev/stderr"
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

echo "[INFO] chr${CHR}: total GEL WGS variants = ${N_VARIANTS}"

if [[ "${N_VARIANTS}" -eq 0 ]]; then
    echo "[ERROR] No GEL WGS variants extracted from ${PVAR}." >&2
    exit 1
fi


################################################################################
# Step 12. Identify variants overlapping ncRNA/pseudogene exons
#
# VARIANT_GENE columns:
#
# Variant:
#    1  variant_CHROM
#    2  variant_start
#    3  variant_end
#    4  variant_ID
#    5  REF
#    6  ALT
#
# Gene/exon:
#    7  exon_CHROM
#    8  exon_start             [0-based BED]
#    9  exon_end
#   10  ensembl_gene_id
#   11  external_gene_name
#   12  gene_biotype
#   13  strand
################################################################################

bedtools intersect \
    -a "${VARIANT_BED}" \
    -b "${EXON_BED}" \
    -wa \
    -wb \
    > "${VARIANT_GENE}"

N_VARIANT_GENE=$(wc -l < "${VARIANT_GENE}")

if [[ "${N_VARIANT_GENE}" -eq 0 ]]; then
    echo "[ERROR] No GEL WGS variants overlap ncRNA/pseudogene exons on chr${CHR}." >&2
    exit 1
fi

N_UNIQUE_EXONIC_VARIANTS=$(
    cut -f4 "${VARIANT_GENE}" \
    | LC_ALL=C sort \
        --parallel="${SORT_THREADS}" \
        -S "${SORT_MEM}" \
        -T "${SORT_TMP}" \
        -u \
    | wc -l
)

N_GENES_WITH_VARIANTS=$(
    cut -f10 "${VARIANT_GENE}" \
    | LC_ALL=C sort \
        --parallel="${SORT_THREADS}" \
        -S "${SORT_MEM}" \
        -T "${SORT_TMP}" \
        -u \
    | wc -l
)

echo
echo "[INFO] chr${CHR}: variant-gene pairs                 = ${N_VARIANT_GENE}"
echo "[INFO] chr${CHR}: unique exonic GEL WGS variants    = ${N_UNIQUE_EXONIC_VARIANTS}"
echo "[INFO] chr${CHR}: genes containing >=1 WGS variant = ${N_GENES_WITH_VARIANTS}"


################################################################################
# Step 13. Verify variant-gene key uniqueness
#
# Multiple different genes per variant are legitimate.
# Duplicate ID + ensembl_gene_id pairs are not expected because exon intervals
# have already been merged within each Ensembl gene.
################################################################################

N_DUPLICATE_VARIANT_GENE=$(
    awk -F '\t' '
    {
        key = $4 "\t" $10
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
    ' "${VARIANT_GENE}"
)

echo "[INFO] chr${CHR}: duplicated variant-gene keys      = ${N_DUPLICATE_VARIANT_GENE}"

if [[ "${N_DUPLICATE_VARIANT_GENE}" -ne 0 ]]; then
    echo "[ERROR] Duplicate variant_ID + ensembl_gene_id pairs detected." >&2
    exit 1
fi


################################################################################
# Step 14. Create unique exonic variant-ID list
################################################################################

cut -f4 "${VARIANT_GENE}" \
| LC_ALL=C sort \
    --parallel="${SORT_THREADS}" \
    -S "${SORT_MEM}" \
    -T "${SORT_TMP}" \
    -u \
> "${EXONIC_VARIANT_IDS}"

N_EXONIC_IDS=$(wc -l < "${EXONIC_VARIANT_IDS}")

if [[ "${N_EXONIC_IDS}" -ne "${N_UNIQUE_EXONIC_VARIANTS}" ]]; then
    echo "[ERROR] Unique exonic variant-ID count is inconsistent." >&2
    exit 1
fi

echo "[PASS] Unique exonic variant-ID list generated."


################################################################################
# Step 15. Extract non-missing transcript-level CADD_PHRED values from GEL VEP
#
# Only variants in EXONIC_VARIANT_IDS are retained.
#
# Missing values:
#
#   NA
#   .
#   empty
#
# are ignored.
#
# Output:
#
#   ID    CADD_PHRED
#
# There may be multiple rows per variant because GEL VEP is transcript-level.
################################################################################

awk \
    -F '\t' \
    -v OFS='\t' \
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
    score = $cadd_col

    if (!(id in keep)) {
        next
    }

    if (score == "" || score == "." || toupper(score) == "NA") {
        next
    }

    if (score !~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/) {
        print "[ERROR] Non-numeric non-missing CADD_PHRED for " id ": " score > "/dev/stderr"
        exit 1
    }

    print id, score
}

' \
"${EXONIC_VARIANT_IDS}" \
"${VEP}" \
> "${CADD_TRANSCRIPT_VALUES}"

N_CADD_TRANSCRIPT_ROWS=$(wc -l < "${CADD_TRANSCRIPT_VALUES}")

echo
echo "[INFO] chr${CHR}: non-missing transcript-level CADD rows = ${N_CADD_TRANSCRIPT_ROWS}"

if [[ "${N_CADD_TRANSCRIPT_ROWS}" -eq 0 ]]; then
    echo "[ERROR] No non-missing CADD values were retrieved for exonic variants." >&2
    exit 1
fi


################################################################################
# Step 16. Sort transcript-level CADD values by ID and numeric score
#
# Numeric sorting allows representations such as:
#
#   2
#   2.0
#
# to be treated as the same underlying numerical value during collapse.
################################################################################

LC_ALL=C sort \
    --parallel="${SORT_THREADS}" \
    -S "${SORT_MEM}" \
    -T "${SORT_TMP}" \
    -t $'\t' \
    -k1,1 \
    -k2,2n \
    "${CADD_TRANSCRIPT_VALUES}" \
    > "${CADD_SORTED_VALUES}"


################################################################################
# Step 17. Collapse transcript-level CADD values to one score per variant
#
# Rules:
#
#   1 distinct non-missing numeric value
#       -> retain that value
#
#   >1 distinct non-missing numeric values
#       -> retain maximum
#       -> write conflict information to MULTI_CADD
#
# Numeric-equivalent strings, e.g. 2 and 2.0, are treated as the same value.
################################################################################

printf "ID\tCADD_PHRED\n" > "${CADD_LOOKUP}"

printf "ID\tN_distinct_CADD\tCADD_values\tMAX_CADD_PHRED\n" > "${MULTI_CADD}"

awk \
    -F '\t' \
    -v OFS='\t' \
    -v lookup="${CADD_LOOKUP}" \
    -v conflicts="${MULTI_CADD}" '

function flush_variant() {

    if (current_id == "") {
        return
    }

    print current_id, max_string >> lookup

    if (n_distinct > 1) {
        print current_id, n_distinct, values, max_string >> conflicts
        print "[CADD_CONFLICT] " current_id \
              " n_distinct=" n_distinct \
              " values=" values \
              " max=" max_string \
              > "/dev/stderr"
    }
}

{
    id = $1
    score_string = $2
    score_numeric = $2 + 0

    if (id != current_id) {

        flush_variant()

        current_id = id
        previous_numeric = score_numeric
        max_numeric = score_numeric
        max_string = score_string
        values = score_string
        n_distinct = 1

        next
    }

    if (score_numeric != previous_numeric) {
        n_distinct++
        values = values "," score_string
        previous_numeric = score_numeric
    }

    if (score_numeric >= max_numeric) {
        max_numeric = score_numeric
        max_string = score_string
    }
}

END {
    flush_variant()
}

' "${CADD_SORTED_VALUES}"


N_CADD_LOOKUP=$(
    awk 'END {print NR - 1}' "${CADD_LOOKUP}"
)

N_MULTI_CADD=$(
    awk 'END {print NR - 1}' "${MULTI_CADD}"
)

echo
echo "[INFO] chr${CHR}: unique exonic variants with CADD     = ${N_CADD_LOOKUP}"
echo "[INFO] chr${CHR}: variants with >1 distinct CADD score = ${N_MULTI_CADD}"


################################################################################
# Step 18. Check CADD lookup uniqueness
################################################################################

N_DUPLICATE_LOOKUP_IDS=$(
    awk -F '\t' '
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

echo "[INFO] chr${CHR}: duplicate CADD lookup IDs = ${N_DUPLICATE_LOOKUP_IDS}"

if [[ "${N_DUPLICATE_LOOKUP_IDS}" -ne 0 ]]; then
    echo "[ERROR] Duplicate variant IDs detected in CADD lookup." >&2
    exit 1
fi


################################################################################
# Step 19. Ensure every CADD lookup ID belongs to the exonic variant universe
################################################################################

N_LOOKUP_OUTSIDE_EXONIC=$(
    awk -F '\t' '
    NR == FNR {
        keep[$1] = 1
        next
    }

    FNR == 1 {
        next
    }

    !($1 in keep) {
        n++
    }

    END {
        print n + 0
    }
    ' \
    "${EXONIC_VARIANT_IDS}" \
    "${CADD_LOOKUP}"
)

if [[ "${N_LOOKUP_OUTSIDE_EXONIC}" -ne 0 ]]; then
    echo "[ERROR] ${N_LOOKUP_OUTSIDE_EXONIC} CADD lookup IDs are absent from the exonic variant universe." >&2
    exit 1
fi

echo "[PASS] CADD lookup is concordant with exonic variant IDs."


################################################################################
# Step 20. Annotate all original variant-gene pairs with CADD_PHRED
#
# Every original variant-gene pair is retained.
#
# Variants absent from CADD_LOOKUP receive:
#
#     CADD_PHRED = NA
################################################################################

printf \
"CHROM\tPOS\tID\tREF\tALT\tensembl_gene_id\texternal_gene_name\tgene_biotype\texon_start\texon_end\tstrand\tCADD_PHRED\n" \
    > "${OUT}"

awk \
    -F '\t' \
    -v OFS='\t' '

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
        $1,      \
        $3,      \
        $4,      \
        $5,      \
        $6,      \
        $10,     \
        $11,     \
        $12,     \
        $8 + 1,  \
        $9,      \
        $13,     \
        score
}

' \
"${CADD_LOOKUP}" \
"${VARIANT_GENE}" \
>> "${OUT}"


################################################################################
# Step 21. Final row-count and annotation sanity checks
################################################################################

N_OUTPUT=$(
    awk 'END {print NR - 1}' "${OUT}"
)

N_WITH_CADD=$(
    awk -F '\t' '
    NR > 1 && $12 != "NA" {
        n++
    }

    END {
        print n + 0
    }
    ' "${OUT}"
)

N_MISSING_CADD=$(
    awk -F '\t' '
    NR > 1 && $12 == "NA" {
        n++
    }

    END {
        print n + 0
    }
    ' "${OUT}"
)

echo
echo "[INFO] chr${CHR}: output variant-gene pairs = ${N_OUTPUT}"
echo "[INFO] chr${CHR}: pairs with CADD score      = ${N_WITH_CADD}"
echo "[INFO] chr${CHR}: pairs with CADD = NA       = ${N_MISSING_CADD}"

if [[ "${N_OUTPUT}" -ne "${N_VARIANT_GENE}" ]]; then
    echo "[ERROR] Final output row count does not equal input variant-gene count." >&2
    echo "[ERROR] Input variant-gene pairs : ${N_VARIANT_GENE}" >&2
    echo "[ERROR] Final output rows        : ${N_OUTPUT}" >&2
    exit 1
fi

if [[ $((N_WITH_CADD + N_MISSING_CADD)) -ne "${N_OUTPUT}" ]]; then
    echo "[ERROR] CADD scored + missing rows do not equal total output rows." >&2
    exit 1
fi


################################################################################
# Step 22. Check final variant-gene key uniqueness
################################################################################

N_FINAL_DUPLICATES=$(
    awk -F '\t' '
    NR > 1 {
        key = $3 "\t" $6
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

echo "[INFO] chr${CHR}: duplicate final variant-gene keys = ${N_FINAL_DUPLICATES}"

if [[ "${N_FINAL_DUPLICATES}" -ne 0 ]]; then
    echo "[ERROR] Final CADD table contains duplicate variant-gene keys." >&2
    exit 1
fi


################################################################################
# Step 23. Additional accounting
#
# Because each variant can map to >1 gene, N_WITH_CADD does not necessarily
# equal N_CADD_LOOKUP. Therefore lookup-level and pair-level counts are reported
# separately rather than directly compared.
################################################################################

N_UNIQUE_FINAL_CADD_VARIANTS=$(
    awk -F '\t' '
    NR > 1 && $12 != "NA" {
        seen[$3] = 1
    }

    END {
        for (id in seen) {
            n++
        }

        print n + 0
    }
    ' "${OUT}"
)

echo "[INFO] chr${CHR}: unique final variants with CADD = ${N_UNIQUE_FINAL_CADD_VARIANTS}"

if [[ "${N_UNIQUE_FINAL_CADD_VARIANTS}" -ne "${N_CADD_LOOKUP}" ]]; then
    echo "[ERROR] Number of unique final variants with CADD does not match CADD lookup." >&2
    echo "[ERROR] Lookup variants with CADD : ${N_CADD_LOOKUP}" >&2
    echo "[ERROR] Final variants with CADD  : ${N_UNIQUE_FINAL_CADD_VARIANTS}" >&2
    exit 1
fi


echo "[PASS] CADD annotation sanity checks passed."


################################################################################
# Step 24. Completion summary
################################################################################

echo
echo "======================================================================"
echo "[PASS] GEL ncRNA/pseudogene exon CADD annotation completed"
echo
echo "Chromosome                         : chr${CHR}"
echo "Merged exon intervals              : ${N_EXONS}"
echo "Total GEL WGS variants             : ${N_VARIANTS}"
echo "Unique exonic GEL WGS variants     : ${N_UNIQUE_EXONIC_VARIANTS}"
echo "Genes containing >=1 variant       : ${N_GENES_WITH_VARIANTS}"
echo "Variant-gene pairs                 : ${N_VARIANT_GENE}"
echo "Non-missing transcript CADD rows   : ${N_CADD_TRANSCRIPT_ROWS}"
echo "Unique variants with CADD          : ${N_CADD_LOOKUP}"
echo "Variants with >1 distinct CADD     : ${N_MULTI_CADD}"
echo "Pairs with CADD score              : ${N_WITH_CADD}"
echo "Pairs with CADD = NA               : ${N_MISSING_CADD}"
echo "Duplicate variant-gene keys        : ${N_FINAL_DUPLICATES}"
echo
echo "CADD conflict QC:"
echo "${MULTI_CADD}"
echo
echo "Output:"
echo "${OUT}"
echo
echo "Finished: $(date)"
echo "======================================================================"