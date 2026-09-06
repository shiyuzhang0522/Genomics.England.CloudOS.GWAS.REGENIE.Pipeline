#!/usr/bin/env bash

################################################################################
# Merge GEL AggV3 VEP annotation subshards into chromosome-level TSV files
#
# Description:
#   For autosomes chr1-chr22:
#     1. Identify all VEP annotation subshard TSV files for each chromosome.
#     2. Print the number of subshards detected.
#     3. Parse genomic START and END coordinates from each filename.
#     4. Sort subshards numerically by genomic START coordinate.
#     5. Verify that all subshards have an identical TSV header.
#     6. Concatenate subshards into one chromosome-level TSV file,
#        retaining the header only once.
#     7. Report the first and last genomic regions included.
#     8. Produce a chromosome-level merge summary table.
#
# Input filename format:
#   GEL.VEP.shard{SHARD}.subshard{SUBSHARD}.chr{CHR}:{START}-{END}.tsv
#
# Example:
#   GEL.VEP.shard51.subshard24.chr9:66299445-67401967.tsv
#
# Input directory:
#   /home/vscode/session_data/filesystems/VEP_annotation
#
# Output directory:
#   /home/vscode/session_data/chromosome_level_VEP
#
# Output files:
#   GEL.VEP.chr1.tsv
#   GEL.VEP.chr2.tsv
#   ...
#   GEL.VEP.chr22.tsv
#
#   chromosome_level_VEP.merge_summary.tsv
#
# Author:
#   Shiyu Zhang
#
# Date:
#   2026-09-06
################################################################################

set -euo pipefail


################################################################################
# Configuration
################################################################################

INPUT_DIR="/home/vscode/session_data/filesystems/VEP_annotation"

OUTPUT_DIR="/home/vscode/session_data/chromosome_level_VEP"

SUMMARY_FILE="${OUTPUT_DIR}/chromosome_level_VEP.merge_summary.tsv"


################################################################################
# Initial checks
################################################################################

echo "============================================================"
echo "GEL AggV3 VEP chromosome-level merger"
echo "============================================================"
echo
echo "Input directory : ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo

if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "ERROR: Input directory does not exist:"
    echo "       ${INPUT_DIR}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"


################################################################################
# Initialize summary table
################################################################################

printf "chromosome\tn_subshards\tfirst_region\tlast_region\toutput_file\tstatus\n" \
    > "${SUMMARY_FILE}"


################################################################################
# Counters
################################################################################

TOTAL_SUBSHARDS=0


################################################################################
# Process chromosomes 1-22
################################################################################

for CHR in $(seq 1 22); do

    echo
    echo "============================================================"
    echo "Chromosome ${CHR}"
    echo "============================================================"


    ########################################################################
    # Temporary file containing:
    #
    # START    END    FILE
    #
    # Use mktemp so that filenames are not dependent on concurrent sessions.
    ########################################################################

    FILE_LIST=$(mktemp "${OUTPUT_DIR}/.chr${CHR}.filelist.XXXXXX")

    TMP_OUTPUT=$(mktemp "${OUTPUT_DIR}/.chr${CHR}.merge.XXXXXX")

    cleanup_chr() {
        rm -f "${FILE_LIST:-}" "${TMP_OUTPUT:-}"
    }

    trap cleanup_chr EXIT


    ########################################################################
    # Identify chromosome-specific subshard files and parse coordinates
    ########################################################################

    while IFS= read -r -d '' FILE; do

        BASENAME=$(basename "${FILE}")

        # Remove everything through ".chrN:"
        REGION_PART="${BASENAME##*.chr${CHR}:}"

        # Remove ".tsv"
        REGION_PART="${REGION_PART%.tsv}"

        START="${REGION_PART%%-*}"
        END="${REGION_PART##*-}"


        ####################################################################
        # Validate parsed coordinates
        ####################################################################

        if [[ ! "${START}" =~ ^[0-9]+$ ]] || [[ ! "${END}" =~ ^[0-9]+$ ]]; then
            echo "ERROR: Could not parse genomic coordinates from:"
            echo "       ${BASENAME}" >&2
            exit 1
        fi

        printf "%s\t%s\t%s\n" \
            "${START}" \
            "${END}" \
            "${FILE}" \
            >> "${FILE_LIST}"

    done < <(
        find "${INPUT_DIR}" \
            -maxdepth 1 \
            -type f \
            -name "GEL.VEP.shard*.subshard*.chr${CHR}:*.tsv" \
            -print0
    )


    ########################################################################
    # Count subshards
    ########################################################################

    N_SUBSHARDS=$(wc -l < "${FILE_LIST}")
    N_SUBSHARDS="${N_SUBSHARDS//[[:space:]]/}"

    echo "Number of subshards: ${N_SUBSHARDS}"


    if [[ "${N_SUBSHARDS}" -eq 0 ]]; then
        echo "ERROR: No subshard files found for chromosome ${CHR}." >&2
        exit 1
    fi

    TOTAL_SUBSHARDS=$((TOTAL_SUBSHARDS + N_SUBSHARDS))


    ########################################################################
    # Sort numerically by genomic START coordinate, then END coordinate
    ########################################################################

    sort \
        -t $'\t' \
        -k1,1n \
        -k2,2n \
        "${FILE_LIST}" \
        -o "${FILE_LIST}"


    ########################################################################
    # Extract first and last genomic regions
    ########################################################################

    FIRST_START=$(awk -F'\t' 'NR==1 {print $1}' "${FILE_LIST}")
    FIRST_END=$(awk -F'\t' 'NR==1 {print $2}' "${FILE_LIST}")

    LAST_START=$(awk -F'\t' 'END {print $1}' "${FILE_LIST}")
    LAST_END=$(awk -F'\t' 'END {print $2}' "${FILE_LIST}")

    FIRST_REGION="chr${CHR}:${FIRST_START}-${FIRST_END}"
    LAST_REGION="chr${CHR}:${LAST_START}-${LAST_END}"

    echo "First region        : ${FIRST_REGION}"
    echo "Last region         : ${LAST_REGION}"


    ########################################################################
    # Check for overlapping or incorrectly ordered genomic intervals
    #
    # Gaps are reported as warnings rather than treated as fatal because
    # genomic partitioning does not necessarily require every base to be
    # represented.
    ########################################################################

    PREV_END=""
    N_GAPS=0

    while IFS=$'\t' read -r START END FILE; do

        if [[ -n "${PREV_END}" ]]; then

            if (( START <= PREV_END )); then
                echo
                echo "ERROR: Overlapping or non-increasing regions detected:"
                echo "       previous END = ${PREV_END}"
                echo "       current region = chr${CHR}:${START}-${END}"
                exit 1
            fi

            if (( START > PREV_END + 1 )); then
                N_GAPS=$((N_GAPS + 1))
            fi

        fi

        PREV_END="${END}"

    done < "${FILE_LIST}"

    echo "Interval gaps       : ${N_GAPS}"


    ########################################################################
    # Obtain reference header from first genomic subshard
    ########################################################################

    FIRST_FILE=$(awk -F'\t' 'NR==1 {print $3}' "${FILE_LIST}")

    REFERENCE_HEADER=$(head -n 1 "${FIRST_FILE}")

    if [[ -z "${REFERENCE_HEADER}" ]]; then
        echo "ERROR: Empty header in:"
        echo "       ${FIRST_FILE}" >&2
        exit 1
    fi


    ########################################################################
    # Check expected number of columns
    ########################################################################

    HEADER_COLUMNS=$(
        printf '%s\n' "${REFERENCE_HEADER}" |
        awk -F'\t' '{print NF}'
    )

    echo "Header columns      : ${HEADER_COLUMNS}"

    if [[ "${HEADER_COLUMNS}" -ne 45 ]]; then
        echo "ERROR: Expected 45 columns but detected ${HEADER_COLUMNS}."
        echo "       File: ${FIRST_FILE}" >&2
        exit 1
    fi


    ########################################################################
    # Write header once
    ########################################################################

    printf '%s\n' "${REFERENCE_HEADER}" > "${TMP_OUTPUT}"


    ########################################################################
    # Merge subshards in genomic order
    ########################################################################

    MERGED_COUNT=0

    while IFS=$'\t' read -r START END FILE; do

        BASENAME=$(basename "${FILE}")

        ####################################################################
        # Ensure input is non-empty
        ####################################################################

        if [[ ! -s "${FILE}" ]]; then
            echo "ERROR: Empty input file:"
            echo "       ${FILE}" >&2
            exit 1
        fi


        ####################################################################
        # Confirm header is identical across all subshards
        ####################################################################

        CURRENT_HEADER=$(head -n 1 "${FILE}")

        if [[ "${CURRENT_HEADER}" != "${REFERENCE_HEADER}" ]]; then
            echo
            echo "ERROR: Header mismatch detected:"
            echo "       ${FILE}" >&2
            exit 1
        fi


        ####################################################################
        # Append data rows only.
        #
        # tail -n +2 removes the repeated header but otherwise leaves the
        # original TSV records unchanged.
        ####################################################################

        tail -n +2 "${FILE}" >> "${TMP_OUTPUT}"

        MERGED_COUNT=$((MERGED_COUNT + 1))

    done < "${FILE_LIST}"


    ########################################################################
    # Verify all identified subshards were processed
    ########################################################################

    if [[ "${MERGED_COUNT}" -ne "${N_SUBSHARDS}" ]]; then
        echo "ERROR: Merge count mismatch for chromosome ${CHR}."
        echo "       Expected : ${N_SUBSHARDS}"
        echo "       Processed: ${MERGED_COUNT}" >&2
        exit 1
    fi


    ########################################################################
    # Finalize chromosome-level output atomically
    ########################################################################

    OUTPUT_FILE="${OUTPUT_DIR}/GEL.VEP.chr${CHR}.tsv"

    mv "${TMP_OUTPUT}" "${OUTPUT_FILE}"

    # TMP_OUTPUT no longer exists after mv.
    TMP_OUTPUT=""


    ########################################################################
    # Record summary
    ########################################################################

    printf "chr%s\t%s\t%s\t%s\t%s\tPASS\n" \
        "${CHR}" \
        "${N_SUBSHARDS}" \
        "${FIRST_REGION}" \
        "${LAST_REGION}" \
        "$(basename "${OUTPUT_FILE}")" \
        >> "${SUMMARY_FILE}"


    ########################################################################
    # Report completion
    ########################################################################

    OUTPUT_SIZE=$(du -h "${OUTPUT_FILE}" | cut -f1)

    echo "Subshards merged    : ${MERGED_COUNT}"
    echo "Output size         : ${OUTPUT_SIZE}"
    echo "Output file         : ${OUTPUT_FILE}"
    echo "Status              : PASS"


    ########################################################################
    # Remove chromosome-specific temporary file list
    ########################################################################

    rm -f "${FILE_LIST}"
    FILE_LIST=""

    trap - EXIT

done


################################################################################
# Final summary
################################################################################

echo
echo "============================================================"
echo "Merge completed successfully"
echo "============================================================"
echo
echo "Chromosomes processed : 22"
echo "Total subshards merged: ${TOTAL_SUBSHARDS}"
echo
echo "Chromosome-level files:"
echo "  ${OUTPUT_DIR}/GEL.VEP.chr1.tsv"
echo "  ..."
echo "  ${OUTPUT_DIR}/GEL.VEP.chr22.tsv"
echo
echo "Merge summary:"
echo "  ${SUMMARY_FILE}"
echo
echo "All chromosome-level merges: PASS"
echo "============================================================"