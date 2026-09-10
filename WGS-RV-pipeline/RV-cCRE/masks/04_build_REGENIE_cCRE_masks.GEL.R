#!/usr/bin/env Rscript


################################################################################
# Build REGENIE set-based testing input files for epidermal
# melanocyte-specific cCREs using Genomics England (GEL) WGS annotations
#
# Dataset:
#   Genomics England (GEL) AggV3 WGS
#
# Variant-level annotation inputs:
#   1. CADD
#   2. GERP
#   3. JARVIS
#
# Biological masks:
#
#   CADD:
#       CADD_PHRED >= 20
#
#   GERP:
#       GERP >= 2
#
#   JARVIS:
#       JARVIS >= 0.99
#
#   FUNC_ALL:
#       union of variants passing CADD, GERP, or JARVIS
#
#   ALL:
#       all GEL WGS variants overlapping melanocyte-specific cCREs,
#       irrespective of functional score
#
# For each chromosome, independent REGENIE input files are generated for:
#
#   CADD/
#   GERP/
#   JARVIS/
#   FUNC_ALL/
#   ALL/
#
# Each mask directory contains:
#
#   chrN.<MASK>.annotation.txt
#   chrN.<MASK>.setlist.txt
#   <MASK>.mask.def
#
# REGENIE annotation format:
#
#   variant_id    cCRE_ID    annotation
#
# REGENIE set-list format:
#
#   cCRE_ID       chromosome      cCRE_start      variant_list
#
# REGENIE mask-definition format:
#
#   mask_name     annotation
#
# Important:
#   - Score-specific missing values (NA) are retained in the input tables but
#     do not enter score-thresholded masks.
#   - The complete variant-cCRE universe is taken from the final GERP table,
#     because the GEL GERP annotation pipeline explicitly retains all original
#     variant-cCRE pairs and assigns GERP=NA when the score is unavailable.
#   - Before mask construction, the CADD, GERP, and JARVIS inputs are required
#     to contain exactly the same variant-cCRE universe.
#
# Usage:
#
#   Rscript 04_build_REGENIE_cCRE_masks.GEL.R --chr 22
#
# Required R packages:
#
#   data.table
#   optparse
#
# Author: Shelley
# Date:   2026-09-10
################################################################################



################################################################################
# Step 0. Load packages
################################################################################

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})



################################################################################
# Step 1. Parse arguments
################################################################################

option_list <- list(

    make_option(
        c("--chr"),
        type = "integer",
        help = "Autosome number: 1-22"
    )

)


opt <- parse_args(
    OptionParser(option_list = option_list)
)


if (is.null(opt$chr)) {
    stop("[ERROR] --chr is required")
}


chr_num <- opt$chr


if (
    length(chr_num) != 1L ||
    is.na(chr_num) ||
    chr_num < 1L ||
    chr_num > 22L
) {
    stop("[ERROR] --chr must be an integer from 1 to 22")
}



################################################################################
# Step 2. Define input and output paths
################################################################################

CADD_DIR <- paste0(
    "/home/vscode/session_data/filesystems/",
    "CADD_raw/Melanocyte.cCRE.CADD.scores"
)


GERP_DIR <- paste0(
    "/home/vscode/session_data/filesystems/",
    "GERP_raw/Melanocyte.cCRE.GERP.scores"
)


JARVIS_DIR <- paste0(
    "/home/vscode/session_data/filesystems/",
    "JARVIS_raw/Melanocyte.cCRE.JARVIS.scores"
)


CADD_FILE <- file.path(
    CADD_DIR,
    paste0(
        "chr", chr_num,
        ".Melanocyte.cCRE.variant.CADD.tsv"
    )
)


GERP_FILE <- file.path(
    GERP_DIR,
    paste0(
        "chr", chr_num,
        ".Melanocyte.cCRE.variant.GERP.tsv"
    )
)


JARVIS_FILE <- file.path(
    JARVIS_DIR,
    paste0(
        "chr", chr_num,
        ".Melanocyte.cCRE.variant.JARVIS.tsv"
    )
)


OUT_ROOT <- "/home/vscode/session_data/REGENIE_cCRE_masks"


OUT_BASE <- file.path(
    OUT_ROOT,
    paste0("chr", chr_num)
)


dir.create(
    OUT_BASE,
    recursive = TRUE,
    showWarnings = FALSE
)



################################################################################
# Step 3. Helper functions
################################################################################


#------------------------------------------------------------------------------
# Check that a required input exists and is non-empty
#------------------------------------------------------------------------------

check_file <- function(path) {

    if (
        !file.exists(path) ||
        is.na(file.size(path)) ||
        file.size(path) == 0
    ) {

        stop(
            paste0(
                "[ERROR] Missing or empty file:\n",
                path
            )
        )
    }
}



#------------------------------------------------------------------------------
# Check required columns
#------------------------------------------------------------------------------

check_columns <- function(dt, required, object_name) {

    missing_columns <- setdiff(
        required,
        names(dt)
    )


    if (length(missing_columns) > 0L) {

        stop(
            paste0(
                "[ERROR] ",
                object_name,
                " is missing required columns: ",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            )
        )
    }
}



#------------------------------------------------------------------------------
# Check score column can be interpreted numerically
#
# fread() will normally infer these columns as numeric with NA values.
# This function also safely handles a character representation if necessary.
#------------------------------------------------------------------------------

convert_numeric_score <- function(dt, column_name, object_name) {

    x <- dt[[column_name]]


    if (!is.numeric(x)) {

        suppressWarnings(
            x_numeric <- as.numeric(x)
        )


        bad <- (
            !is.na(x) &
            trimws(as.character(x)) != "" &
            toupper(trimws(as.character(x))) != "NA" &
            trimws(as.character(x)) != "." &
            is.na(x_numeric)
        )


        if (any(bad)) {

            example_values <- unique(
                as.character(x[bad])
            )


            stop(
                paste0(
                    "[ERROR] ",
                    object_name,
                    " contains non-numeric non-missing values in ",
                    column_name,
                    ". Examples: ",
                    paste(
                        head(example_values, 10L),
                        collapse = ", "
                    )
                )
            )
        }


        set(
            dt,
            j = column_name,
            value = x_numeric
        )
    }


    invisible(dt)
}



#------------------------------------------------------------------------------
# Canonical representation of the complete variant-cCRE universe
#
# We use all identifying and cCRE-coordinate fields that should be identical
# among the CADD, GERP, and JARVIS final annotation outputs.
#------------------------------------------------------------------------------

get_variant_ccre_pairs <- function(dt) {

    unique(
        dt[
            ,
            .(
                CHROM,
                POS,
                ID,
                REF,
                ALT,
                cCRE_ID,
                cCRE_start,
                cCRE_end
            )
        ]
    )
}



#------------------------------------------------------------------------------
# Check for duplicated variant-cCRE keys within an annotation table
#
# A variant overlapping two DIFFERENT cCREs is legitimate.
# The duplicate definition therefore includes cCRE_ID.
#------------------------------------------------------------------------------

check_duplicate_pairs <- function(dt, object_name) {

    dup <- dt[
        ,
        .N,
        by = .(
            CHROM,
            POS,
            ID,
            REF,
            ALT,
            cCRE_ID
        )
    ][
        N > 1L
    ]


    if (nrow(dup) > 0L) {

        stop(
            paste0(
                "[ERROR] ",
                object_name,
                " contains ",
                nrow(dup),
                " duplicated variant-cCRE keys"
            )
        )
    }


    invisible(TRUE)
}



#------------------------------------------------------------------------------
# Require two annotation tables to contain exactly the same variant-cCRE
# universe
#------------------------------------------------------------------------------

check_pair_universe <- function(reference_pairs,
                                query_pairs,
                                reference_name,
                                query_name) {

    missing_from_query <- fsetdiff(
        reference_pairs,
        query_pairs
    )


    extra_in_query <- fsetdiff(
        query_pairs,
        reference_pairs
    )


    if (
        nrow(missing_from_query) > 0L ||
        nrow(extra_in_query) > 0L
    ) {

        stop(
            paste0(
                "[ERROR] Variant-cCRE universe mismatch between ",
                reference_name,
                " and ",
                query_name,
                "\n",
                "Pairs present in ",
                reference_name,
                " but absent from ",
                query_name,
                ": ",
                nrow(missing_from_query),
                "\n",
                "Pairs present in ",
                query_name,
                " but absent from ",
                reference_name,
                ": ",
                nrow(extra_in_query)
            )
        )
    }


    invisible(TRUE)
}



#------------------------------------------------------------------------------
# Write REGENIE mask definition
#
# Example:
#
#   CADD    CADD
#------------------------------------------------------------------------------

write_mask_def <- function(mask_dir, mask_name) {

    fwrite(
        data.table(
            mask_name,
            mask_name
        ),
        file = file.path(
            mask_dir,
            paste0(
                mask_name,
                ".mask.def"
            )
        ),
        sep = "\t",
        col.names = FALSE,
        quote = FALSE
    )
}



#------------------------------------------------------------------------------
# Expand a REGENIE set-list back to variant-cCRE pairs
#
# Used only for annotation/set-list concordance QC.
#------------------------------------------------------------------------------

expand_setlist <- function(setlist) {

    setlist[
        ,
        .(
            variant_id = unlist(
                strsplit(
                    variants,
                    ",",
                    fixed = TRUE
                )
            )
        ),
        by = .(
            cCRE_ID,
            CHR,
            cCRE_start
        )
    ]
}



#------------------------------------------------------------------------------
# Write annotation, set-list and mask-definition files for one mask
#------------------------------------------------------------------------------

write_regenie_files <- function(df, mask_name) {


    cat(
        "[INFO] Writing mask:",
        mask_name,
        "\n"
    )


    if (nrow(df) == 0L) {

        stop(
            paste0(
                "[ERROR] Mask ",
                mask_name,
                " contains zero variant-cCRE pairs on chr",
                chr_num
            )
        )
    }


    #--------------------------------------------------------------------------
    # Ensure exact unique variant-cCRE membership
    #--------------------------------------------------------------------------

    df <- unique(
        df[
            ,
            .(
                variant_id,
                cCRE_ID,
                cCRE_start,
                CHR
            )
        ]
    )


    #--------------------------------------------------------------------------
    # Create mask-specific directory
    #--------------------------------------------------------------------------

    mask_dir <- file.path(
        OUT_BASE,
        mask_name
    )


    dir.create(
        mask_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )


    #--------------------------------------------------------------------------
    # Annotation file
    #
    # variant_id    cCRE_ID    annotation
    #--------------------------------------------------------------------------

    annotation <- unique(
        df[
            ,
            .(
                variant_id,
                cCRE_ID,
                annotation = mask_name
            )
        ]
    )


    annotation_file <- file.path(
        mask_dir,
        paste0(
            "chr",
            chr_num,
            ".",
            mask_name,
            ".annotation.txt"
        )
    )


    fwrite(
        annotation,
        file = annotation_file,
        sep = "\t",
        col.names = FALSE,
        quote = FALSE
    )


    #--------------------------------------------------------------------------
    # Set-list
    #
    # cCRE_ID    chromosome    cCRE_start    comma-separated variants
    #--------------------------------------------------------------------------

    tmp <- unique(
        df[
            ,
            .(
                cCRE_ID,
                CHR,
                cCRE_start,
                variant_id
            )
        ]
    )


    #--------------------------------------------------------------------------
    # Each cCRE must map to exactly one chromosome and one start coordinate
    #--------------------------------------------------------------------------

    bad_coordinates <- tmp[
        ,
        .(
            n_chr = uniqueN(CHR),
            n_start = uniqueN(cCRE_start)
        ),
        by = cCRE_ID
    ][
        n_chr > 1L |
        n_start > 1L
    ]


    if (nrow(bad_coordinates) > 0L) {

        stop(
            paste0(
                "[ERROR] ",
                nrow(bad_coordinates),
                " cCREs have inconsistent chromosome/start coordinates in ",
                mask_name
            )
        )
    }


    setlist <- tmp[
        ,
        .(
            CHR = CHR[1L],
            cCRE_start = cCRE_start[1L],
            variants = paste(
                unique(variant_id),
                collapse = ","
            )
        ),
        by = cCRE_ID
    ]


    setlist_file <- file.path(
        mask_dir,
        paste0(
            "chr",
            chr_num,
            ".",
            mask_name,
            ".setlist.txt"
        )
    )


    fwrite(
        setlist[
            ,
            .(
                cCRE_ID,
                CHR,
                cCRE_start,
                variants
            )
        ],
        file = setlist_file,
        sep = "\t",
        col.names = FALSE,
        quote = FALSE
    )


    #--------------------------------------------------------------------------
    # Mask definition
    #--------------------------------------------------------------------------

    write_mask_def(
        mask_dir,
        mask_name
    )


    #--------------------------------------------------------------------------
    # Annotation/set-list concordance QC
    #--------------------------------------------------------------------------

    annotation_pairs <- unique(
        annotation[
            ,
            .(
                variant_id,
                cCRE_ID
            )
        ]
    )


    setlist_pairs <- unique(
        expand_setlist(setlist)[
            ,
            .(
                variant_id,
                cCRE_ID
            )
        ]
    )


    missing_from_setlist <- fsetdiff(
        annotation_pairs,
        setlist_pairs
    )


    extra_in_setlist <- fsetdiff(
        setlist_pairs,
        annotation_pairs
    )


    if (
        nrow(missing_from_setlist) > 0L ||
        nrow(extra_in_setlist) > 0L
    ) {

        stop(
            paste0(
                "[ERROR] Annotation/set-list mismatch for mask ",
                mask_name,
                "\n",
                "Annotation pairs absent from set-list: ",
                nrow(missing_from_setlist),
                "\n",
                "Set-list pairs absent from annotation: ",
                nrow(extra_in_setlist)
            )
        )
    }


    #--------------------------------------------------------------------------
    # Return mask summary
    #--------------------------------------------------------------------------

    data.table(
        mask = mask_name,
        n_pairs = nrow(annotation),
        n_variants = uniqueN(annotation$variant_id),
        n_cCRE = uniqueN(annotation$cCRE_ID)
    )
}



################################################################################
# Step 4. Initial report
################################################################################

cat("\n")
cat("======================================================================\n")
cat("[INFO] Building GEL REGENIE melanocyte-cCRE masks\n")
cat("[INFO] Chromosome :", chr_num, "\n")
cat("[INFO] Started    :", date(), "\n")
cat("======================================================================\n")

cat("\n[INFO] Input files:\n")
cat("CADD   :", CADD_FILE, "\n")
cat("GERP   :", GERP_FILE, "\n")
cat("JARVIS :", JARVIS_FILE, "\n")
cat("Output :", OUT_BASE, "\n\n")



################################################################################
# Step 5. Check required input files
################################################################################

for (
    f in c(
        CADD_FILE,
        GERP_FILE,
        JARVIS_FILE
    )
) {

    check_file(f)
}



################################################################################
# Step 6. Read final GEL annotation tables
################################################################################

cat("[INFO] Reading CADD annotation table\n")

CADD <- fread(
    CADD_FILE,
    na.strings = c(
        "NA",
        ".",
        ""
    )
)


check_columns(
    CADD,
    c(
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "cCRE_ID",
        "cCRE_start",
        "cCRE_end",
        "CADD_PHRED"
    ),
    "CADD"
)


convert_numeric_score(
    CADD,
    "CADD_PHRED",
    "CADD"
)



cat("[INFO] Reading GERP annotation table\n")

GERP <- fread(
    GERP_FILE,
    na.strings = c(
        "NA",
        ".",
        ""
    )
)


check_columns(
    GERP,
    c(
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "cCRE_ID",
        "cCRE_start",
        "cCRE_end",
        "GERP"
    ),
    "GERP"
)


convert_numeric_score(
    GERP,
    "GERP",
    "GERP"
)



cat("[INFO] Reading JARVIS annotation table\n")

JARVIS <- fread(
    JARVIS_FILE,
    na.strings = c(
        "NA",
        ".",
        ""
    )
)


check_columns(
    JARVIS,
    c(
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "cCRE_ID",
        "cCRE_start",
        "cCRE_end",
        "JARVIS"
    ),
    "JARVIS"
)


convert_numeric_score(
    JARVIS,
    "JARVIS",
    "JARVIS"
)



################################################################################
# Step 7. Input-level chromosome and duplicate QC
################################################################################

cat("\n[INFO] Performing input-level QC\n")


for (
    object_name in c(
        "CADD",
        "GERP",
        "JARVIS"
    )
) {

    dt <- get(object_name)


    observed_chr <- unique(
        dt$CHROM[
            !is.na(dt$CHROM)
        ]
    )


    if (
        length(observed_chr) != 1L ||
        observed_chr[1L] != chr_num
    ) {

        stop(
            paste0(
                "[ERROR] Unexpected chromosome value(s) in ",
                object_name,
                ". Expected chr",
                chr_num,
                "; observed: ",
                paste(
                    observed_chr,
                    collapse = ","
                )
            )
        )
    }


    check_duplicate_pairs(
        dt,
        object_name
    )


    cat(
        "[INFO] ",
        object_name,
        " rows = ",
        nrow(dt),
        "; unique IDs = ",
        uniqueN(dt$ID),
        "; unique cCREs = ",
        uniqueN(dt$cCRE_ID),
        "\n",
        sep = ""
    )
}



################################################################################
# Step 8. Verify CADD/GERP/JARVIS contain the same variant-cCRE universe
################################################################################

cat("\n")
cat("[INFO] Checking cross-annotation variant-cCRE concordance\n")


CADD_PAIRS <- get_variant_ccre_pairs(CADD)

GERP_PAIRS <- get_variant_ccre_pairs(GERP)

JARVIS_PAIRS <- get_variant_ccre_pairs(JARVIS)


check_pair_universe(
    GERP_PAIRS,
    CADD_PAIRS,
    "GERP",
    "CADD"
)


check_pair_universe(
    GERP_PAIRS,
    JARVIS_PAIRS,
    "GERP",
    "JARVIS"
)


cat(
    "[PASS] CADD, GERP, and JARVIS contain identical variant-cCRE universes\n"
)


N_ALL_PAIRS <- nrow(GERP_PAIRS)

N_ALL_VARIANTS <- uniqueN(GERP_PAIRS$ID)

N_ALL_CCRE <- uniqueN(GERP_PAIRS$cCRE_ID)


cat(
    "[INFO] Complete variant-cCRE pairs :",
    N_ALL_PAIRS,
    "\n"
)

cat(
    "[INFO] Complete unique variants    :",
    N_ALL_VARIANTS,
    "\n"
)

cat(
    "[INFO] Complete unique cCREs       :",
    N_ALL_CCRE,
    "\n"
)



################################################################################
# Step 9. Construct biological masks
################################################################################

cat("\n")
cat("[INFO] Constructing biological masks\n")



#------------------------------------------------------------------------------
# CADD mask
#
# CADD_PHRED >= 20
#------------------------------------------------------------------------------

CADD_MASK <- unique(
    CADD[
        !is.na(CADD_PHRED) &
        CADD_PHRED >= 20,
        .(
            variant_id = ID,
            cCRE_ID,
            cCRE_start,
            CHR = CHROM
        )
    ]
)



#------------------------------------------------------------------------------
# GERP mask
#
# GERP >= 2
#------------------------------------------------------------------------------

GERP_MASK <- unique(
    GERP[
        !is.na(GERP) &
        GERP >= 2,
        .(
            variant_id = ID,
            cCRE_ID,
            cCRE_start,
            CHR = CHROM
        )
    ]
)



#------------------------------------------------------------------------------
# JARVIS mask
#
# JARVIS >= 0.99
#------------------------------------------------------------------------------

JARVIS_MASK <- unique(
    JARVIS[
        !is.na(JARVIS) &
        JARVIS >= 0.99,
        .(
            variant_id = ID,
            cCRE_ID,
            cCRE_start,
            CHR = CHROM
        )
    ]
)



#------------------------------------------------------------------------------
# FUNC_ALL
#
# Union:
#
#   CADD OR GERP OR JARVIS
#------------------------------------------------------------------------------

FUNC_ALL <- unique(
    rbindlist(
        list(
            CADD_MASK,
            GERP_MASK,
            JARVIS_MASK
        ),
        use.names = TRUE
    )
)



#------------------------------------------------------------------------------
# ALL
#
# Every variant-cCRE pair, regardless of CADD/GERP/JARVIS score.
#
# The complete GERP output is used because the GEL GERP pipeline explicitly
# retains all variant-cCRE pairs and uses NA for missing GERP scores.
#------------------------------------------------------------------------------

ALL_MASK <- unique(
    GERP[
        ,
        .(
            variant_id = ID,
            cCRE_ID,
            cCRE_start,
            CHR = CHROM
        )
    ]
)



################################################################################
# Step 10. Mask-level QC
################################################################################

cat("\n")
cat("[INFO] Performing mask-level QC\n")


#------------------------------------------------------------------------------
# Recalculate the expected functional union independently
#------------------------------------------------------------------------------

EXPECTED_FUNC_ALL <- unique(
    rbindlist(
        list(
            CADD_MASK,
            GERP_MASK,
            JARVIS_MASK
        ),
        use.names = TRUE
    )
)


MISSING_FROM_FUNC_ALL <- fsetdiff(
    EXPECTED_FUNC_ALL,
    FUNC_ALL
)


EXTRA_IN_FUNC_ALL <- fsetdiff(
    FUNC_ALL,
    EXPECTED_FUNC_ALL
)


if (
    nrow(MISSING_FROM_FUNC_ALL) > 0L ||
    nrow(EXTRA_IN_FUNC_ALL) > 0L
) {

    stop(
        "[ERROR] FUNC_ALL union check failed"
    )
}



#------------------------------------------------------------------------------
# Every score-based mask must be a subset of ALL
#------------------------------------------------------------------------------

ALL_MEMBERSHIP <- unique(
    ALL_MASK[
        ,
        .(
            variant_id,
            cCRE_ID,
            cCRE_start,
            CHR
        )
    ]
)


for (
    mask_name in c(
        "CADD_MASK",
        "GERP_MASK",
        "JARVIS_MASK",
        "FUNC_ALL"
    )
) {

    mask_dt <- get(mask_name)


    outside_all <- fsetdiff(
        unique(mask_dt),
        ALL_MEMBERSHIP
    )


    if (nrow(outside_all) > 0L) {

        stop(
            paste0(
                "[ERROR] ",
                mask_name,
                " contains ",
                nrow(outside_all),
                " variant-cCRE pairs absent from ALL"
            )
        )
    }
}


cat("[PASS] All score-based masks are subsets of ALL\n")



################################################################################
# Step 11. Report mask sizes before writing
################################################################################

cat("\n")
cat("[INFO] Mask sizes before writing\n")


prewrite_summary <- rbindlist(
    list(

        data.table(
            mask = "CADD",
            n_pairs = nrow(CADD_MASK),
            n_variants = uniqueN(CADD_MASK$variant_id),
            n_cCRE = uniqueN(CADD_MASK$cCRE_ID)
        ),

        data.table(
            mask = "GERP",
            n_pairs = nrow(GERP_MASK),
            n_variants = uniqueN(GERP_MASK$variant_id),
            n_cCRE = uniqueN(GERP_MASK$cCRE_ID)
        ),

        data.table(
            mask = "JARVIS",
            n_pairs = nrow(JARVIS_MASK),
            n_variants = uniqueN(JARVIS_MASK$variant_id),
            n_cCRE = uniqueN(JARVIS_MASK$cCRE_ID)
        ),

        data.table(
            mask = "FUNC_ALL",
            n_pairs = nrow(FUNC_ALL),
            n_variants = uniqueN(FUNC_ALL$variant_id),
            n_cCRE = uniqueN(FUNC_ALL$cCRE_ID)
        ),

        data.table(
            mask = "ALL",
            n_pairs = nrow(ALL_MASK),
            n_variants = uniqueN(ALL_MASK$variant_id),
            n_cCRE = uniqueN(ALL_MASK$cCRE_ID)
        )

    ),
    use.names = TRUE
)


print(prewrite_summary)



################################################################################
# Step 12. Write REGENIE files
################################################################################

cat("\n")
cat("[INFO] Writing REGENIE input files\n")


summary <- rbindlist(
    list(
        write_regenie_files(
            CADD_MASK,
            "CADD"
        ),

        write_regenie_files(
            GERP_MASK,
            "GERP"
        ),

        write_regenie_files(
            JARVIS_MASK,
            "JARVIS"
        ),

        write_regenie_files(
            FUNC_ALL,
            "FUNC_ALL"
        ),

        write_regenie_files(
            ALL_MASK,
            "ALL"
        )
    ),
    use.names = TRUE
)



################################################################################
# Step 13. Verify written summary against pre-write summary
################################################################################

setorder(
    summary,
    mask
)


setorder(
    prewrite_summary,
    mask
)


if (!identical(summary, prewrite_summary)) {

    stop(
        "[ERROR] Written mask summary does not match pre-write mask summary"
    )
}



################################################################################
# Step 14. Save chromosome-level mask summary
################################################################################

SUMMARY_FILE <- file.path(
    OUT_BASE,
    paste0(
        "chr",
        chr_num,
        ".mask_summary.tsv"
    )
)


fwrite(
    summary,
    file = SUMMARY_FILE,
    sep = "\t",
    quote = FALSE
)



################################################################################
# Step 15. Save chromosome-level input/QC summary
################################################################################

QC_SUMMARY <- data.table(

    chromosome = chr_num,

    input_variant_cCRE_pairs = N_ALL_PAIRS,

    input_unique_variants = N_ALL_VARIANTS,

    input_unique_cCREs = N_ALL_CCRE,

    CADD_nonmissing_pairs = CADD[
        !is.na(CADD_PHRED),
        .N
    ],

    GERP_nonmissing_pairs = GERP[
        !is.na(GERP),
        .N
    ],

    JARVIS_nonmissing_pairs = JARVIS[
        !is.na(JARVIS),
        .N
    ],

    CADD_mask_pairs = nrow(CADD_MASK),

    GERP_mask_pairs = nrow(GERP_MASK),

    JARVIS_mask_pairs = nrow(JARVIS_MASK),

    FUNC_ALL_pairs = nrow(FUNC_ALL),

    ALL_pairs = nrow(ALL_MASK)

)


QC_FILE <- file.path(
    OUT_BASE,
    paste0(
        "chr",
        chr_num,
        ".input_QC_summary.tsv"
    )
)


fwrite(
    QC_SUMMARY,
    file = QC_FILE,
    sep = "\t",
    quote = FALSE
)



################################################################################
# Step 16. Final report
################################################################################

cat("\n")
cat("======================================================================\n")
cat("[PASS] GEL REGENIE melanocyte-cCRE mask construction completed\n")
cat("\n")

cat(
    "Chromosome                  : chr",
    chr_num,
    "\n",
    sep = ""
)

cat(
    "Complete variant-cCRE pairs : ",
    N_ALL_PAIRS,
    "\n",
    sep = ""
)

cat(
    "Complete unique variants    : ",
    N_ALL_VARIANTS,
    "\n",
    sep = ""
)

cat(
    "Complete unique cCREs       : ",
    N_ALL_CCRE,
    "\n",
    sep = ""
)

cat("\nMask summary:\n")
print(summary)

cat("\nOutput directory:\n")
cat(OUT_BASE, "\n")

cat("\nMask summary:\n")
cat(SUMMARY_FILE, "\n")

cat("\nInput/QC summary:\n")
cat(QC_FILE, "\n")

cat("\nFinished:", date(), "\n")
cat("======================================================================\n")