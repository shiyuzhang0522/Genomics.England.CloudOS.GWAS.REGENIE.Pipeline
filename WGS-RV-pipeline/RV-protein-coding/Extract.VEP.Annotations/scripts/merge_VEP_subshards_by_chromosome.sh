#!/usr/bin/env bash

################################################################################
# Parallel merger of GEL AggV3 VEP annotation subshards
#
# Description:
#   Merge GEL AggV3 VEP annotation subshards into chromosome-level TSV files
#   for autosomes chr1-chr22.
#
#   Optimizations:
#     1. Scan the input directory only once.
#     2. Construct chromosome-specific ordered manifests.
#     3. Process multiple chromosomes concurrently.
#     4. Preserve genomic order within each chromosome.
#     5. Validate all headers before performing the large merge.
#     6. Write to temporary files and atomically rename only after success.
#     7. Generate independent per-chromosome logs.
#     8. Generate a final chromosome-level merge summary.
#
# Input filename format:
#   GEL.VEP.shard{SHARD}.subshard{SUBSHARD}.chr{CHR}:{START}-{END}.tsv
#
# Example:
#   GEL.VEP.shard51.subshard24.chr9:66299445-67401967.tsv
#
# Input:
#   /home/vscode/session_data/filesystems/VEP_annotation
#
# Output:
#   /home/vscode/session_data/chromosome_level_VEP
#
# Outputs:
#   GEL.VEP.chr1.tsv
#   ...
#   GEL.VEP.chr22.tsv
#
#   chromosome_level_VEP.merge_summary.tsv
#
# Parallelism:
#   Default: 4 chromosomes concurrently.
#
#   Can be overridden at runtime:
#
#       MAX_PARALLEL=4 ./merge_VEP_subshards_by_chromosome.parallel.sh
#
# Author:
#   Shiyu Zhang
#
# Date:
#   2026-09-06
################################################################################

set -Eeuo pipefail

export LC_ALL=C


################################################################################
# Configuration
################################################################################

INPUT_DIR="/home/vscode/session_data/filesystems/VEP_annotation"

OUTPUT_DIR="/home/vscode/session_data/chromosome_level_VEP"

SUMMARY_FILE="${OUTPUT_DIR}/chromosome_level_VEP.merge_summary.tsv"

LOG_DIR="${OUTPUT_DIR}/logs"

MAX_PARALLEL="${MAX_PARALLEL:-4}"


################################################################################
# Basic validation
################################################################################

if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "ERROR: Input directory does not exist:" >&2
    echo "       ${INPUT_DIR}" >&2
    exit 1
fi

if [[ ! "${MAX_PARALLEL}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: MAX_PARALLEL must be a positive integer." >&2
    exit 1
fi

# wait -n requires Bash >= 4.3
if (( BASH_VERSINFO[0] < 4 ||
      (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 3) )); then

    echo "ERROR: Bash >= 4.3 is required." >&2
    echo "Current Bash: ${BASH_VERSION}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"


################################################################################
# Temporary metadata directory
################################################################################

META_DIR=$(mktemp -d "${OUTPUT_DIR}/.merge_metadata.XXXXXX")

cleanup_global() {
    rm -rf "${META_DIR}"
}

trap cleanup_global EXIT


################################################################################
# Header
################################################################################

echo "============================================================"
echo "GEL AggV3 VEP chromosome-level parallel merger"
echo "============================================================"
echo
echo "Input directory      : ${INPUT_DIR}"
echo "Output directory     : ${OUTPUT_DIR}"
echo "Parallel chromosomes : ${MAX_PARALLEL}"
echo "Bash version         : ${BASH_VERSION}"
echo


################################################################################
# Initialize chromosome-specific metadata files
################################################################################

for CHR in $(seq 1 22); do
    : > "${META_DIR}/chr${CHR}.list"
done


################################################################################
# Discover input files ONCE
#
# Metadata columns:
#
#   START    END    FILE
################################################################################

echo "Scanning input directory..."

N_DISCOVERED=0
N_AUTOSOMAL=0

while IFS= read -r -d '' FILE; do

    BASENAME="${FILE##*/}"

    ########################################################################
    # Parse:
    #
    # GEL.VEP.shard51.subshard24.chr9:66299445-67401967.tsv
    ########################################################################

    if [[ "${BASENAME}" =~ ^GEL\.VEP\.shard[0-9]+\.subshard[0-9]+\.chr([0-9]+):([0-9]+)-([0-9]+)\.tsv$ ]]; then

        CHR="${BASH_REMATCH[1]}"
        START="${BASH_REMATCH[2]}"
        END="${BASH_REMATCH[3]}"

        N_DISCOVERED=$((N_DISCOVERED + 1))

        ####################################################################
        # Autosomes only
        ####################################################################

        if (( CHR >= 1 && CHR <= 22 )); then

            printf "%s\t%s\t%s\n" \
                "${START}" \
                "${END}" \
                "${FILE}" \
                >> "${META_DIR}/chr${CHR}.list"

            N_AUTOSOMAL=$((N_AUTOSOMAL + 1))
        fi

    fi

done < <(
    find "${INPUT_DIR}" \
        -maxdepth 1 \
        -type f \
        -name 'GEL.VEP.shard*.subshard*.chr*.tsv' \
        -print0
)

echo "Recognized VEP files : ${N_DISCOVERED}"
echo "Autosomal subshards  : ${N_AUTOSOMAL}"
echo


################################################################################
# Sort chromosome manifests and validate counts
################################################################################

echo "Preparing chromosome manifests..."
echo

TOTAL_SUBSHARDS=0

for CHR in $(seq 1 22); do

    LIST="${META_DIR}/chr${CHR}.list"

    sort \
        -t $'\t' \
        -k1,1n \
        -k2,2n \
        "${LIST}" \
        -o "${LIST}"

    N_SUBSHARDS=$(wc -l < "${LIST}")
    N_SUBSHARDS="${N_SUBSHARDS//[[:space:]]/}"

    if [[ "${N_SUBSHARDS}" -eq 0 ]]; then
        echo "ERROR: No input subshards found for chromosome ${CHR}." >&2
        exit 1
    fi

    TOTAL_SUBSHARDS=$((TOTAL_SUBSHARDS + N_SUBSHARDS))

    printf "chr%-2s : %s subshards\n" \
        "${CHR}" \
        "${N_SUBSHARDS}"

done

echo
echo "Total autosomal subshards: ${TOTAL_SUBSHARDS}"
echo


################################################################################
# Chromosome processing function
#
# Executed in a subshell so each chromosome has independent:
#   - variables
#   - temporary output
#   - cleanup trap
################################################################################

process_chromosome() (

    set -Eeuo pipefail

    CHR="$1"

    LIST="${META_DIR}/chr${CHR}.list"

    OUTPUT_FILE="${OUTPUT_DIR}/GEL.VEP.chr${CHR}.tsv"

    CHR_SUMMARY="${META_DIR}/chr${CHR}.summary"

    TMP_OUTPUT=""


    ########################################################################
    # Cleanup incomplete output
    ########################################################################

    cleanup_chr() {

        if [[ -n "${TMP_OUTPUT}" && -f "${TMP_OUTPUT}" ]]; then
            rm -f "${TMP_OUTPUT}"
        fi
    }

    trap cleanup_chr EXIT INT TERM


    ########################################################################
    # Count subshards
    ########################################################################

    N_SUBSHARDS=$(wc -l < "${LIST}")
    N_SUBSHARDS="${N_SUBSHARDS//[[:space:]]/}"


    echo "============================================================"
    echo "Chromosome ${CHR}"
    echo "============================================================"
    echo "Number of subshards : ${N_SUBSHARDS}"


    ########################################################################
    # First / last interval
    ########################################################################

    FIRST_START=$(awk -F'\t' 'NR==1 {print $1}' "${LIST}")
    FIRST_END=$(awk -F'\t' 'NR==1 {print $2}' "${LIST}")

    LAST_START=$(awk -F'\t' 'END {print $1}' "${LIST}")
    LAST_END=$(awk -F'\t' 'END {print $2}' "${LIST}")

    FIRST_REGION="chr${CHR}:${FIRST_START}-${FIRST_END}"
    LAST_REGION="chr${CHR}:${LAST_START}-${LAST_END}"

    echo "First region        : ${FIRST_REGION}"
    echo "Last region         : ${LAST_REGION}"


    ########################################################################
    # Check genomic ordering and overlapping intervals
    ########################################################################

    PREV_END=""
    N_GAPS=0

    while IFS=$'\t' read -r START END FILE; do

        if [[ -n "${PREV_END}" ]]; then

            if (( START <= PREV_END )); then

                echo "ERROR: Overlapping or non-increasing intervals." >&2
                echo "Previous END : ${PREV_END}" >&2
                echo "Current      : chr${CHR}:${START}-${END}" >&2

                exit 1
            fi

            if (( START > PREV_END + 1 )); then
                N_GAPS=$((N_GAPS + 1))
            fi
        fi

        PREV_END="${END}"

    done < "${LIST}"

    echo "Interval gaps       : ${N_GAPS}"


    ########################################################################
    # Determine reference header
    ########################################################################

    FIRST_FILE=$(awk -F'\t' 'NR==1 {print $3}' "${LIST}")

    if [[ ! -s "${FIRST_FILE}" ]]; then
        echo "ERROR: First input file is missing or empty:" >&2
        echo "       ${FIRST_FILE}" >&2
        exit 1
    fi

    REFERENCE_HEADER=$(head -n 1 -- "${FIRST_FILE}")

    if [[ -z "${REFERENCE_HEADER}" ]]; then
        echo "ERROR: Empty header:" >&2
        echo "       ${FIRST_FILE}" >&2
        exit 1
    fi

    HEADER_COLUMNS=$(
        printf '%s\n' "${REFERENCE_HEADER}" |
        awk -F'\t' '{print NF}'
    )

    echo "Header columns      : ${HEADER_COLUMNS}"

    if [[ "${HEADER_COLUMNS}" -ne 45 ]]; then

        echo "ERROR: Expected 45 columns, detected ${HEADER_COLUMNS}." >&2
        echo "       ${FIRST_FILE}" >&2

        exit 1
    fi


    ########################################################################
    # PRE-MERGE validation
    #
    # Validate every input before writing a potentially ~100-GB chromosome
    # output.
    ########################################################################

    echo "Validating input headers..."

    VALIDATED_COUNT=0

    while IFS=$'\t' read -r START END FILE; do

        if [[ ! -s "${FILE}" ]]; then

            echo "ERROR: Missing or empty input file:" >&2
            echo "       ${FILE}" >&2

            exit 1
        fi

        CURRENT_HEADER=$(head -n 1 -- "${FILE}")

        if [[ "${CURRENT_HEADER}" != "${REFERENCE_HEADER}" ]]; then

            echo "ERROR: Header mismatch:" >&2
            echo "       ${FILE}" >&2

            exit 1
        fi

        VALIDATED_COUNT=$((VALIDATED_COUNT + 1))

    done < "${LIST}"

    if [[ "${VALIDATED_COUNT}" -ne "${N_SUBSHARDS}" ]]; then

        echo "ERROR: Validation count mismatch." >&2
        echo "Expected : ${N_SUBSHARDS}" >&2
        echo "Validated: ${VALIDATED_COUNT}" >&2

        exit 1
    fi

    echo "Headers validated   : ${VALIDATED_COUNT}"


    ########################################################################
    # Protect existing completed output
    ########################################################################

    if [[ -e "${OUTPUT_FILE}" ]]; then

        echo "ERROR: Final output already exists:" >&2
        echo "       ${OUTPUT_FILE}" >&2
        echo "Remove it explicitly before rerunning chromosome ${CHR}." >&2

        exit 1
    fi


    ########################################################################
    # Create chromosome-specific temporary output
    ########################################################################

    TMP_OUTPUT=$(
        mktemp "${OUTPUT_DIR}/.GEL.VEP.chr${CHR}.merge.XXXXXX"
    )


    ########################################################################
    # Write header exactly once
    ########################################################################

    printf '%s\n' "${REFERENCE_HEADER}" > "${TMP_OUTPUT}"


    ########################################################################
    # Merge input files in genomic order
    ########################################################################

    echo "Starting merge..."

    MERGED_COUNT=0

    while IFS=$'\t' read -r START END FILE; do

        tail -n +2 -- "${FILE}" >> "${TMP_OUTPUT}"

        MERGED_COUNT=$((MERGED_COUNT + 1))

    done < "${LIST}"


    ########################################################################
    # Validate merge count
    ########################################################################

    if [[ "${MERGED_COUNT}" -ne "${N_SUBSHARDS}" ]]; then

        echo "ERROR: Merge count mismatch." >&2
        echo "Expected: ${N_SUBSHARDS}" >&2
        echo "Merged  : ${MERGED_COUNT}" >&2

        exit 1
    fi


    ########################################################################
    # Atomic finalization
    ########################################################################

    mv "${TMP_OUTPUT}" "${OUTPUT_FILE}"

    TMP_OUTPUT=""


    ########################################################################
    # Output metadata
    #
    # stat is metadata-only; unlike wc -l it does not reread the entire
    # potentially 100-GB output.
    ########################################################################

    OUTPUT_BYTES=$(stat -c '%s' "${OUTPUT_FILE}")

    OUTPUT_HUMAN=$(du -h "${OUTPUT_FILE}" | cut -f1)


    ########################################################################
    # Per-chromosome summary
    ########################################################################

    printf \
        "chr%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n" \
        "${CHR}" \
        "${N_SUBSHARDS}" \
        "${FIRST_REGION}" \
        "${LAST_REGION}" \
        "${N_GAPS}" \
        "$(basename "${OUTPUT_FILE}")" \
        "${OUTPUT_BYTES}" \
        > "${CHR_SUMMARY}"


    ########################################################################
    # Done
    ########################################################################

    echo "Subshards merged    : ${MERGED_COUNT}"
    echo "Output size         : ${OUTPUT_HUMAN}"
    echo "Output file         : ${OUTPUT_FILE}"
    echo "Status              : PASS"
    echo

)


################################################################################
# Parallel scheduler
################################################################################

echo "============================================================"
echo "Launching chromosome merges"
echo "Maximum parallel jobs: ${MAX_PARALLEL}"
echo "============================================================"
echo

RUNNING=0

for CHR in $(seq 1 22); do

    CHR_LOG="${LOG_DIR}/GEL.VEP.chr${CHR}.merge.log"

    echo "Launching chr${CHR}"
    echo "  log: ${CHR_LOG}"

    process_chromosome "${CHR}" \
        > "${CHR_LOG}" \
        2>&1 &

    RUNNING=$((RUNNING + 1))


    ########################################################################
    # Keep at most MAX_PARALLEL chromosome processes active
    ########################################################################

    if (( RUNNING >= MAX_PARALLEL )); then

        if ! wait -n; then

            echo
            echo "ERROR: A chromosome merge failed." >&2
            echo "Check logs in:" >&2
            echo "       ${LOG_DIR}" >&2

            ################################################################
            # Stop remaining chromosome merges
            ################################################################

            for PID in $(jobs -pr); do
                kill "${PID}" 2>/dev/null || true
            done

            wait 2>/dev/null || true

            exit 1
        fi

        RUNNING=$((RUNNING - 1))
    fi

done


################################################################################
# Wait for remaining chromosomes
################################################################################

while (( RUNNING > 0 )); do

    if ! wait -n; then

        echo
        echo "ERROR: A chromosome merge failed." >&2
        echo "Check logs in:" >&2
        echo "       ${LOG_DIR}" >&2

        for PID in $(jobs -pr); do
            kill "${PID}" 2>/dev/null || true
        done

        wait 2>/dev/null || true

        exit 1
    fi

    RUNNING=$((RUNNING - 1))

done


################################################################################
# Verify all 22 chromosome summaries exist
################################################################################

for CHR in $(seq 1 22); do

    if [[ ! -s "${META_DIR}/chr${CHR}.summary" ]]; then

        echo "ERROR: Missing successful summary for chromosome ${CHR}." >&2
        exit 1
    fi

done


################################################################################
# Construct final master summary atomically
################################################################################

TMP_SUMMARY=$(mktemp "${OUTPUT_DIR}/.merge_summary.XXXXXX")

printf \
    "chromosome\tn_subshards\tfirst_region\tlast_region\tinterval_gaps\toutput_file\toutput_bytes\tstatus\n" \
    > "${TMP_SUMMARY}"

for CHR in $(seq 1 22); do
    cat "${META_DIR}/chr${CHR}.summary" >> "${TMP_SUMMARY}"
done

mv "${TMP_SUMMARY}" "${SUMMARY_FILE}"


################################################################################
# Final validation
################################################################################

N_FINAL_OUTPUTS=$(
    find "${OUTPUT_DIR}" \
        -maxdepth 1 \
        -type f \
        -name 'GEL.VEP.chr*.tsv' \
        | wc -l
)

N_FINAL_OUTPUTS="${N_FINAL_OUTPUTS//[[:space:]]/}"

if [[ "${N_FINAL_OUTPUTS}" -ne 22 ]]; then

    echo "ERROR: Expected 22 final chromosome TSVs." >&2
    echo "Detected: ${N_FINAL_OUTPUTS}" >&2

    exit 1
fi


################################################################################
# Final report
################################################################################

echo
echo "============================================================"
echo "Parallel merge completed successfully"
echo "============================================================"
echo
echo "Chromosomes completed : 22"
echo "Autosomal subshards   : ${TOTAL_SUBSHARDS}"
echo "Parallel jobs         : ${MAX_PARALLEL}"
echo
echo "Output directory:"
echo "  ${OUTPUT_DIR}"
echo
echo "Master summary:"
echo "  ${SUMMARY_FILE}"
echo
echo "Per-chromosome logs:"
echo "  ${LOG_DIR}"
echo
echo "All chromosome-level merges: PASS"
echo "============================================================"