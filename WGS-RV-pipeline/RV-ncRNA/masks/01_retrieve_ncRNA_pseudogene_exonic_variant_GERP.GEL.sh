#!/usr/bin/env bash

################################################################################
# Retrieve variant-level GERP conservation scores for Genomics England (GEL)
# WGS variants located within annotated exons of autosomal ncRNA/pseudogene
# genes.
#
# Dataset:
#   Genomics England (GEL) AggV3 WGS
#
# Variant source:
#   Chromosome-specific DRAGEN PVAR files:
#
#     /home/vscode/session_data/filesystems/
#       chrom-msvcf/chrom-${CHR}/postproc-pgen/dragen.pvar
#
# ncRNA/pseudogene exon annotation:
#   Ensembl release 116, GRCh38
#
#   Ensembl116.GRCh38.autosomal.ncRNA_pseudogene.
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
# Coordinate conventions:
#   - GEL PVAR chromosome names:      1, 2, ..., 22
#   - GERP BigWig chromosome names:  1, 2, ..., 22
#   - Ensembl exon coordinates:       1-based, inclusive
#   - BED / BigWig coordinates:       0-based, half-open
#   - GEL PVAR POS:                   1-based
#
# Variant-to-gene mapping:
#   - Each GEL WGS variant is represented by its POS as a 1-bp BED interval:
#
#         [POS-1, POS)
#
#   - Exonic intervals are the union of all annotated exons across transcripts
#     within each Ensembl ncRNA/pseudogene gene.
#
#   - A variant may legitimately map to more than one gene when exon intervals
#     from different genes overlap.
#
# GERP annotation:
#   - GERP is position-specific rather than allele-specific.
#   - Multiple ALT alleles at the same genomic position therefore receive the
#     same positional GERP score.
#   - All original variant-gene pairs are retained.
#   - Missing positional GERP annotations are reported as GERP = NA.
#   - No GERP threshold is applied in this script.
#   - GERP >= 2 will be applied later during mask construction.
#
# FILTER policy:
#   - All variants present in dragen.pvar are considered, irrespective of the
#     FILTER column. This is consistent with the GEL cCRE annotation workflow.
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
#   GERP
#
# Usage:
#
#   bash 01_retrieve_ncRNA_pseudogene_exonic_variant_GERP.GEL.sh 22
#
# Required software:
#   awk
#   bedtools
#   bigWigToBedGraph
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
    bigWigToBedGraph \
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
# Step 3. Define input and output paths
################################################################################

PVAR="/home/vscode/session_data/filesystems/chrom-msvcf/chrom-${CHR}/postproc-pgen/dragen.pvar"


NCRNA_EXONS="/home/vscode/session_data/filesystems/Ensembl116.GRCh38.autosomal.ncRNA_pseudogene.gene_union_exons.merged.tsv"


GERP_BW="/home/vscode/session_data/filesystems/gerp_conservation_scores.homo_sapiens.GRCh38.bw"


OUT_DIR="/home/vscode/session_data/GERP_ncRNA"


LOG_DIR="${OUT_DIR}/logs"


mkdir -p "${OUT_DIR}"
mkdir -p "${LOG_DIR}"



################################################################################
# Step 4. Define intermediate and final files
################################################################################

EXON_BED="${OUT_DIR}/chr${CHR}.ncRNA_pseudogene.exons.bed"


VARIANT_BED="${OUT_DIR}/chr${CHR}.variants.bed"


VARIANT_GENE="${OUT_DIR}/chr${CHR}.variants.in.ncRNA_pseudogene.exons.bed"


GERP_ANNOTATED="${OUT_DIR}/chr${CHR}.variants.in.ncRNA_pseudogene.exons.GERP.raw.bed"


OUT="${OUT_DIR}/chr${CHR}.ncRNA_pseudogene.exonic_variants.GERP.tsv"



################################################################################
# Step 5. Initial report
################################################################################

echo
echo "======================================================================"
echo "[INFO] GEL ncRNA/pseudogene exon GERP annotation"
echo "[INFO] Chromosome : chr${CHR}"
echo "[INFO] Started    : $(date)"
echo "======================================================================"

echo
echo "[INFO] Input PVAR:"
echo "       ${PVAR}"

echo "[INFO] ncRNA/pseudogene exon annotation:"
echo "       ${NCRNA_EXONS}"

echo "[INFO] GERP BigWig:"
echo "       ${GERP_BW}"

echo "[INFO] Output directory:"
echo "       ${OUT_DIR}"



################################################################################
# Step 6. Check required input files
################################################################################

for file in \
    "${PVAR}" \
    "${NCRNA_EXONS}" \
    "${GERP_BW}"

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
# Step 7. Validate PVAR header
#
# Expected:
#
#   #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
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
echo "[INFO] PVAR header:"
echo "       ${PVAR_HEADER}"



################################################################################
# Step 8. Extract chromosome-specific ncRNA/pseudogene merged exons
#
# Input TSV:
#
#   1  ensembl_gene_id
#   2  external_gene_name
#   3  gene_biotype
#   4  chromosome_name
#   5  exon_start            [1-based inclusive]
#   6  exon_end              [1-based inclusive]
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
# Empty gene names are replaced with NA to preserve the field structure.
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

    if (
        gene_name == "" ||
        gene_name == "."
    ) {
        gene_name = "NA"
    }


    if (
        $1 == "" ||
        $3 == "" ||
        $5 !~ /^[0-9]+$/ ||
        $6 !~ /^[0-9]+$/ ||
        $5 > $6
    ) {

        print \
            "[ERROR] Invalid ncRNA exon record:", \
            $0 \
            > "/dev/stderr"

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
# Step 9. Convert GEL PVAR variants to 1-bp BED intervals
#
# PVAR:
#
#   #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
#
# BED:
#
#   1  CHROM
#   2  POS - 1
#   3  POS
#   4  ID
#   5  REF
#   6  ALT
#
# Important:
#   No FILTER-based exclusion is applied.
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

        print \
            "[ERROR] Unexpected chromosome in chr" chr \
            " PVAR: " $1 \
            > "/dev/stderr"

        exit 1
    }


    if ($2 !~ /^[0-9]+$/ || $2 < 1) {

        print \
            "[ERROR] Invalid POS in PVAR: " $2 \
            > "/dev/stderr"

        exit 1
    }


    if (
        $3 == "" ||
        $3 == "." ||
        $4 == "" ||
        $4 == "." ||
        $5 == "" ||
        $5 == "."
    ) {

        print \
            "[ERROR] Missing ID/REF/ALT in PVAR record: " $0 \
            > "/dev/stderr"

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

    echo "[ERROR] No GEL WGS variants extracted from:" >&2
    echo "        ${PVAR}" >&2
    exit 1

fi



################################################################################
# Step 10. Identify GEL WGS variants overlapping ncRNA/pseudogene exons
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
#    8  exon_start            [0-based BED]
#    9  exon_end
#   10  ensembl_gene_id
#   11  external_gene_name
#   12  gene_biotype
#   13  strand
#
# Because exon intervals were previously merged within each Ensembl gene,
# transcript-level exon redundancy should not create duplicate variant-gene
# mappings.
################################################################################

bedtools intersect \
    -a "${VARIANT_BED}" \
    -b "${EXON_BED}" \
    -wa \
    -wb \
    > "${VARIANT_GENE}"


N_VARIANT_GENE=$(wc -l < "${VARIANT_GENE}")


N_UNIQUE_EXONIC_VARIANTS=$(
    cut -f4 "${VARIANT_GENE}" \
    | LC_ALL=C sort -u \
    | wc -l
)


N_GENES_WITH_VARIANTS=$(
    cut -f10 "${VARIANT_GENE}" \
    | LC_ALL=C sort -u \
    | wc -l
)


echo
echo "[INFO] chr${CHR}: variant-gene pairs                  = ${N_VARIANT_GENE}"
echo "[INFO] chr${CHR}: unique exonic GEL WGS variants     = ${N_UNIQUE_EXONIC_VARIANTS}"
echo "[INFO] chr${CHR}: genes containing >=1 WGS variant  = ${N_GENES_WITH_VARIANTS}"


if [[ "${N_VARIANT_GENE}" -eq 0 ]]; then

    echo "[ERROR] No GEL WGS variants overlap ncRNA/pseudogene exons on chr${CHR}." >&2
    exit 1

fi



################################################################################
# Step 11. Check variant-gene pair uniqueness
#
# A variant is allowed to map to multiple DIFFERENT genes.
#
# What is not expected is more than one row for the same:
#
#       variant_ID + ensembl_gene_id
#
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


echo "[INFO] chr${CHR}: duplicated variant-gene keys       = ${N_DUPLICATE_VARIANT_GENE}"


if [[ "${N_DUPLICATE_VARIANT_GENE}" -ne 0 ]]; then

    echo "[ERROR] Duplicate variant_ID + ensembl_gene_id pairs detected." >&2
    echo "[ERROR] This suggests unexpected overlapping exon intervals within a gene." >&2
    exit 1

fi



################################################################################
# Step 12. Retrieve positional GERP scores
#
# GERP BigWig chromosome convention:
#
#       1, 2, ..., 22
#
# This directly matches VARIANT_GENE; therefore no chromosome-name conversion
# is required.
#
# bigWigToBedGraph output:
#
#   1  GERP_CHROM
#   2  GERP_start
#   3  GERP_end
#   4  GERP_score
#
# VARIANT_GENE contains 13 columns.
#
# After bedtools intersect -loj:
#
#   14  GERP_CHROM
#   15  GERP_start
#   16  GERP_end
#   17  GERP_score
#
# -loj guarantees that all original variant-gene pairs are retained.
################################################################################

bigWigToBedGraph \
    -chrom="${CHR}" \
    "${GERP_BW}" \
    stdout \
| bedtools intersect \
    -a "${VARIANT_GENE}" \
    -b stdin \
    -wa \
    -wb \
    -loj \
    > "${GERP_ANNOTATED}"



################################################################################
# Step 13. Generate final variant-gene GERP table
#
# exon_start is converted back from BED:
#
#       0-based exon_start -> 1-based Ensembl exon_start
################################################################################

printf \
"CHROM\tPOS\tID\tREF\tALT\tensembl_gene_id\texternal_gene_name\tgene_biotype\texon_start\texon_end\tstrand\tGERP\n" \
    > "${OUT}"


awk \
    -F '\t' \
    -v OFS='\t' '

{

    gerp = $17


    if (
        $14 == "." ||
        gerp == "." ||
        gerp == ""
    ) {
        gerp = "NA"
    }


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
        gerp
}

' "${GERP_ANNOTATED}" >> "${OUT}"



################################################################################
# Step 14. Final row-count and annotation sanity checks
################################################################################

N_OUTPUT=$(
    awk 'END {print NR - 1}' "${OUT}"
)


N_GERP_SCORED=$(
    awk \
        -F '\t' '

    NR > 1 && $12 != "NA" {
        n++
    }

    END {
        print n + 0
    }

    ' "${OUT}"
)


N_GERP_MISSING=$(
    awk \
        -F '\t' '

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
echo "[INFO] chr${CHR}: pairs with GERP score      = ${N_GERP_SCORED}"
echo "[INFO] chr${CHR}: pairs with GERP = NA       = ${N_GERP_MISSING}"



#------------------------------------------------------------------------------
# Every original variant-gene pair must be retained.
#------------------------------------------------------------------------------

if [[ "${N_OUTPUT}" -ne "${N_VARIANT_GENE}" ]]; then

    echo "[ERROR] Final output row count does not equal input variant-gene count." >&2
    echo "[ERROR] Input variant-gene pairs : ${N_VARIANT_GENE}" >&2
    echo "[ERROR] Final output rows        : ${N_OUTPUT}" >&2

    exit 1

fi



#------------------------------------------------------------------------------
# Scored + missing must account for every final row.
#------------------------------------------------------------------------------

if [[ $((N_GERP_SCORED + N_GERP_MISSING)) -ne "${N_OUTPUT}" ]]; then

    echo "[ERROR] GERP scored + missing rows do not equal total output rows." >&2

    exit 1

fi



################################################################################
# Step 15. Check final variant-gene key uniqueness
################################################################################

N_FINAL_DUPLICATES=$(
    awk \
        -F '\t' '

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

    echo "[ERROR] Final GERP table contains duplicate variant-gene keys." >&2
    exit 1

fi


echo "[PASS] GERP annotation sanity checks passed."



################################################################################
# Step 16. Completion summary
################################################################################

echo
echo "======================================================================"
echo "[PASS] GEL ncRNA/pseudogene exon GERP annotation completed"
echo
echo "Chromosome                         : chr${CHR}"
echo "Merged exon intervals              : ${N_EXONS}"
echo "Total GEL WGS variants             : ${N_VARIANTS}"
echo "Unique exonic GEL WGS variants     : ${N_UNIQUE_EXONIC_VARIANTS}"
echo "Genes containing >=1 variant       : ${N_GENES_WITH_VARIANTS}"
echo "Variant-gene pairs                 : ${N_VARIANT_GENE}"
echo "Pairs with GERP score              : ${N_GERP_SCORED}"
echo "Pairs with GERP = NA               : ${N_GERP_MISSING}"
echo "Duplicate variant-gene keys        : ${N_FINAL_DUPLICATES}"
echo
echo "Output:"
echo "${OUT}"
echo
echo "Finished: $(date)"
echo "======================================================================"