#!/usr/bin/env Rscript


################################################################################
# Build REGENIE gene-based testing input files for autosomal ncRNA/pseudogene
# exons using Genomics England (GEL) WGS annotations
#
# Dataset:
#   Genomics England (GEL) AggV3 WGS
#
# Reference annotation:
#   Ensembl release 116, GRCh38, autosomal gene-union merged exons
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
#       union of variant-gene pairs passing CADD, GERP, or JARVIS
#
#   ALL:
#       all ncRNA/pseudogene exonic variant-gene pairs, irrespective of score
#
# UKBB-compatible output structure beneath --out-root:
#
#   REGENIE_inputs/chrN/CADD/
#   REGENIE_inputs/chrN/GERP/
#   REGENIE_inputs/chrN/JARVIS/
#   REGENIE_inputs/chrN/FUNC_ALL/
#   REGENIE_inputs/chrN/ALL/
#
# Each mask directory contains:
#
#   chrN.<MASK>.annotation.txt
#   chrN.<MASK>.setlist.txt
#   <MASK>.mask.def
#
# REGENIE annotation format:
#
#   variant_id    ensembl_gene_id    annotation
#
# REGENIE set-list format:
#
#   ensembl_gene_id    chromosome    gene_start    variant_list
#
# REGENIE mask-definition format:
#
#   mask_name    annotation
#
# Gene anchor:
#   Minimum 1-based exon_start across the ORIGINAL Ensembl gene annotation,
#   independent of strand and mask membership. No coordinate conversion.
#
# Important:
#   - Preserve original GEL variant IDs and mappings to multiple genes.
#   - Missing scores do not enter the corresponding score-thresholded mask.
#   - Before mask construction, require identical ID x ensembl_gene_id sets
#     across CADD, GERP, and JARVIS, including rows with missing scores.
#   - Reject duplicate pairs and inconsistent shared annotation metadata.
#   - Derive ALL from the complete GERP table only after concordance passes.
#   - No allele-frequency filtering is performed during mask creation.
#   - Empty masks are errors, consistent with the supplied UKBB builder.
#
# Output handling:
#   Stage and validate files before publishing the chromosome directory.
#   Save diagnostic reports on failure. Preserve existing chromosome outputs;
#   use another --out-root for reruns. Different chromosomes may run in parallel.
#   Do not run the same chromosome concurrently against the same output root.
#
# Usage:
#
#   Rscript 05_build_REGENIE_ncRNA_pseudogene_exonic_masks.GEL.R --chr 22
#
# Required R packages:
#
#   data.table
#   optparse
#
# Author: Shelley
# Date:   2026-09-11
# Version: 1.0.0 (formatting revision; analysis logic unchanged)
################################################################################


################################################################################
# Step 0. Load packages
################################################################################

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})
options(warn = 2) # Malformed reads must not be silently accepted.
main <- function() {


    ################################################################################
    # Step 1. Parse arguments
    ################################################################################

    opt <- parse_args(
        OptionParser(
            option_list = list(
                make_option(
                    "--chr",
                    type = "integer",
                    help = "Autosome, 1-22 [required]"
                ),
                make_option(
                    "--input-root",
                    dest = "input_root",
                    type = "character",
                    default = "/home/vscode/session_data/filesystems",
                    help = "Directory containing CADD_ncRNA, GERP_ncRNA, JARVIS_ncRNA [%default]"
                ),
                make_option(
                    "--exon-file",
                    dest = "exon_file",
                    type = "character",
                    default = NULL,
                    help = "Merged exon TSV; default: Ensembl116 TSV under input-root"
                ),
                make_option(
                    "--out-root",
                    dest = "out_root",
                    type = "character",
                    default = "/home/vscode/session_data/REGENIE_ncRNA_masks",
                    help = "Parent of UKBB-style REGENIE_inputs directory [%default]"
                ),
                make_option(
                    "--threads",
                    type = "integer",
                    default = 2L,
                    help = "data.table threads [%default]"
                )
            )
        )
    )
    fail <- function(...) stop(
        paste0(...),
        call. = FALSE
    )
    log <- function(...) cat(
        format(
            Sys.time(),
            "[%Y-%m-%d %H:%M:%S] "
        ),
        paste0(...),
        "\n",
        sep = ""
    )
    if (
        is.null(opt$chr) || is.na(opt$chr) || !opt$chr %in% 1:22
    )
    fail("--chr must be an integer from 1 to 22.")
    if (
        is.na(opt$threads) || opt$threads < 1L
    ) fail("--threads must be positive.")
    setDTthreads(opt$threads)


    ################################################################################
    # Step 2. Define input and output paths
    ################################################################################

    chr <- opt$chr
    prefix <- paste0(
        "chr",
        chr
    )
    scores <- c(
        CADD = "CADD_PHRED",
        GERP = "GERP",
        JARVIS = "JARVIS"
    )
    thresholds <- c(
        CADD = 20,
        GERP = 2,
        JARVIS = 0.99
    )
    mask_names <- c(
        names(scores),
        "FUNC_ALL",
        "ALL"
    )
    paths <- setNames(
        vapply(
            names(scores),
            function(s) file.path(
                opt$input_root,
                paste0(
                    s,
                    "_ncRNA"
                ),
                paste0(
                    prefix,
                    ".ncRNA_pseudogene.exonic_variants.",
                    s,
                    ".tsv"
                )
            ),
            character(1)
        ),
        names(scores)
    )
    exon_path <- opt$exon_file
    if (
        is.null(exon_path)
    ) exon_path <- file.path(
        opt$input_root,
        "Ensembl116.GRCh38.autosomal.ncRNA_pseudogene.gene_union_exons.merged.tsv"
    )
    parent <- file.path(
        opt$out_root,
        "REGENIE_inputs"
    )
    final <- file.path(
        parent,
        prefix
    )
    if (
        file.exists(final)
    ) fail(
        "Output already exists: ",
        final,
        ". Choose another --out-root; existing results have been preserved."
    )
    dir.create(
        parent,
        recursive = TRUE,
        showWarnings = FALSE
    )
    stage <- tempfile(
        paste0(
            prefix,
            ".building."
        ),
        tmpdir = parent
    )
    if (
        !dir.create(stage)
    ) fail(
        "Cannot create staging directory: ",
        stage
    )
    published <- FALSE
    on.exit(
        {
            if (
                !published
            ) message(
                "[FAILED] Diagnostic files retained at: ",
                stage
            )
        },
        add = TRUE
    )


    ################################################################################
    # Step 3. Helper functions
    ################################################################################

    report <- function(
        dt,
        suffix
    ) fwrite(
        dt,
        file.path(
            stage,
            paste0(
                prefix,
                ".",
                suffix
            )
        ),
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )


    #------------------------------------------------------------------------------
    # Read a complete tab-delimited input and check required columns
    #------------------------------------------------------------------------------

    read_tsv <- function(
        path,
        required
    ) {
        if (
            !file.exists(path) || is.na(
                file.size(path)
            ) || file.size(path) == 0
        )
        fail(
            "Missing or empty input: ",
            path
        )
        dt <- fread(
            path,
            sep = "\t",
            header = TRUE,
            quote = "",
            na.strings = c(
                "NA",
                ".",
                ""
            ),
            colClasses = "character",
            showProgress = interactive()
        )
        if (
            anyDuplicated(
                names(dt)
            )
        ) fail(
            "Duplicate header names: ",
            path
        )
        missing <- setdiff(
            required,
            names(dt)
        )
        if (
            length(missing)
        ) fail(
            "Missing columns in ",
            path,
            ": ",
            paste(
                missing,
                collapse = ", "
            )
        )
        if (
            !nrow(dt)
        ) fail(
            "No data rows: ",
            path
        )
        dt
    }


    #------------------------------------------------------------------------------
    # Validate numeric columns while preserving allowed missing scores
    #------------------------------------------------------------------------------

    numeric_column <- function(
        dt,
        col,
        label,
        missing_ok = FALSE
    ) {
        raw <- dt[[col]]
        value <- suppressWarnings(
            as.numeric(raw)
        )
        bad <- (
            !is.na(raw) & (
                is.na(value) | !is.finite(value)
            )
        ) |
        (
            !missing_ok & is.na(raw)
        )
        if (
            any(bad)
        ) fail(
            label,
            ": invalid ",
            col,
            " at data row ",
            which(bad)[
                1
            ]
        )
        set(
            dt,
            j = col,
            value = value
        )
    }


    #------------------------------------------------------------------------------
    # Validate chromosome and coordinate columns
    #------------------------------------------------------------------------------

    integers <- function(
        dt,
        cols,
        label
    ) {
        for (
            col in cols
        ) {
            numeric_column(
                dt,
                col,
                label
            )
            x <- dt[[col]]
            if (
                any(
                    x < 1 | x != floor(x) | x > .Machine$integer.max
                )
            )
            fail(
                label,
                ": invalid positive integer in ",
                col
            )
            set(
                dt,
                j = col,
                value = as.integer(x)
            )
        }
    }


    #------------------------------------------------------------------------------
    # Validate identifiers for headerless REGENIE output
    #------------------------------------------------------------------------------

    identifiers <- function(
        dt,
        cols,
        label
    ) {
        for (
            col in cols
        ) {
            x <- dt[[col]]
            if (
                anyNA(x) || any(
                    grepl(
                        "[[:space:],]",
                        x
                    )
                ) || any(x == "")
            )
            fail(
                label,
                ": missing/unsafe identifier in ",
                col
            )
        }
    }


    #------------------------------------------------------------------------------
    # Canonical representation of unique variant-gene pairs
    #------------------------------------------------------------------------------

    pairs <- function(dt) unique(
        dt[
            ,
            .(
                ID,
                ensembl_gene_id
            )
        ]
    )


    ################################################################################
    # Step 4. Initial report
    ################################################################################

    log(
        "[START] GEL ncRNA/pseudogene masks: ",
        prefix
    )
    log(
        "[INFO] Output: ",
        final
    )
    ################################################################################
    # Step 5. Read complete score tables and check duplicate pairs
    ################################################################################
    common <- c(
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "ensembl_gene_id",
        "external_gene_name",
        "gene_biotype",
        "exon_start",
        "exon_end",
        "strand"
    )
    inputs <- list()
    input_qc <- list()
    duplicate_failure <- FALSE
    for (
        s in names(scores)
    ) {
        log(
            "[READ] ",
            s,
            ": ",
            paths[[s]]
        )
        dt <- read_tsv(
            paths[[s]],
            c(
                common,
                scores[[s]]
            )
        )
        identifiers(
            dt,
            c(
                "ID",
                "ensembl_gene_id",
                "REF",
                "ALT"
            ),
            s
        )
        integers(
            dt,
            c(
                "CHROM",
                "POS",
                "exon_start",
                "exon_end"
            ),
            s
        )
        if (
            any(dt$CHROM != chr)
        ) fail(
            s,
            ": unexpected chromosome."
        )
        if (
            anyNA(dt$strand) || any(
                !dt$strand %in% c(
                    "+",
                    "-"
                )
            )
        )
        fail(
            s,
            ": strand must be + or -."
        )
        if (
            any(dt$exon_start > dt$exon_end)
        ) fail(
            s,
            ": exon_start exceeds exon_end."
        )
        numeric_column(
            dt,
            scores[[s]],
            s,
            missing_ok = TRUE
        )
        duplicates <- dt[
            ,
            .N,
            by = .(
                ID,
                ensembl_gene_id
            )
        ][
            N > 1L
        ]
        report(
            duplicates,
            paste0(
                s,
                ".duplicate_pairs.tsv"
            )
        )
        duplicate_failure <- duplicate_failure || nrow(duplicates) > 0L
        input_qc[[s]] <- data.table(
            source = s,
            n_rows = nrow(dt),
            n_pairs = nrow(
                pairs(dt)
            ),
            n_variants = uniqueN(dt$ID),
            n_genes = uniqueN(dt$ensembl_gene_id),
            n_nonmissing_score_pairs = sum(
                !is.na(
                    dt[[
                        scores[[s]]
                    ]]
                )
            ),
            n_duplicate_keys = nrow(duplicates)
        )
        inputs[[s]] <- dt
        log(
            "[INFO] ",
            s,
            ": ",
            nrow(dt),
            " rows; ",
            uniqueN(dt$ID),
            " variants; ",
            uniqueN(dt$ensembl_gene_id),
            " genes"
        )
    }


    report(
        rbindlist(input_qc),
        "input_QC_summary.tsv"
    )
    ################################################################################
    # Step 6. Verify cross-annotation variant-gene concordance
    ################################################################################
    universe <- pairs(inputs$GERP)
    concordance <- list()
    for (
        s in c(
            "CADD",
            "JARVIS"
        )
    ) {
        p <- pairs(
            inputs[[s]]
        )
        missing <- fsetdiff(
            universe,
            p
        )
        extra <- fsetdiff(
            p,
            universe
        )
        report(
            missing,
            paste0(
                s,
                ".missing_vs_GERP.tsv"
            )
        )
        report(
            extra,
            paste0(
                s,
                ".extra_vs_GERP.tsv"
            )
        )
        concordance[[s]] <- data.table(
            reference = "GERP",
            source = s,
            missing_pairs = nrow(missing),
            extra_pairs = nrow(extra)
        )
    }


    concordance <- rbindlist(concordance)
    report(
        concordance,
        "pair_concordance.tsv"
    )
    if (
        duplicate_failure || any(
            concordance$missing_pairs + concordance$extra_pairs > 0L
        )
    )
    fail(
        "Input pair QC failed. See duplicate and pair-concordance TSV reports."
    )
    log(
        "[PASS] All three variant-gene sets are identical and duplicate-free."
    )
    ################################################################################
    # Step 7. Verify shared annotation metadata
    ################################################################################
    for (
        s in names(scores)
    ) setkeyv(
        inputs[[s]],
        c(
            "ID",
            "ensembl_gene_id"
        )
    )
    base <- inputs$GERP
    metadata_failure <- FALSE
    for (
        s in c(
            "CADD",
            "JARVIS"
        )
    ) {
        bad <- rep(
            FALSE,
            nrow(base)
        )
        for (
            col in common
        ) {
            a <- base[[col]];
            b <- inputs[[s]][[col]]
            equal <- (
                is.na(a) & is.na(b)
            ) | (
                !is.na(a) & !is.na(b) & a == b
            )
            bad <- bad | !equal
        }
        report(
            base[
                bad,
                .(
                    ID,
                    ensembl_gene_id
                )
            ],
            paste0(
                s,
                ".metadata_mismatch_pairs.tsv"
            )
        )
        metadata_failure <- metadata_failure || any(bad)
    }


    if (
        metadata_failure
    ) fail(
        "Shared annotation metadata differs between score tables."
    )
    variant_coordinates <- unique(
        base[
            ,
            .(
                ID,
                CHROM,
                POS,
                REF,
                ALT
            )
        ]
    )
    if (
        anyDuplicated(variant_coordinates$ID)
    ) fail(
        "An ID maps to inconsistent variant coordinates/alleles."
    )
    ################################################################################
    # Step 8. Read Ensembl exons and construct stable gene anchors
    ################################################################################
    log(
        "[READ] Ensembl 116 exon reference: ",
        exon_path
    )
    exons <- read_tsv(
        exon_path,
        c(
            "ensembl_gene_id",
            "external_gene_name",
            "gene_biotype",
            "chromosome_name",
            "exon_start",
            "exon_end",
            "strand"
        )
    )
    identifiers(
        exons,
        "ensembl_gene_id",
        "Exons"
    )
    integers(
        exons,
        c(
            "chromosome_name",
            "exon_start",
            "exon_end"
        ),
        "Exons"
    )
    if (
        any(!exons$chromosome_name %in% 1:22) || any(exons$exon_start > exons$exon_end)
    )
    fail("Invalid autosomal exon coordinates.")
    if (
        anyNA(exons$strand) || any(
            !exons$strand %in% c(
                "+",
                "-"
            )
        )
    )
    fail("Exon strand must be + or -.")
    gene_qc <- exons[
        ,
        .(
            n_chr = uniqueN(chromosome_name),
            n_strand = uniqueN(strand)
        ),
        by = ensembl_gene_id
    ]
    if (
        nrow(
            gene_qc[
                n_chr != 1L | n_strand != 1L
            ]
        )
    ) fail(
        "Inconsistent gene chromosome/strand in exon reference."
    )
    exons <- exons[
        chromosome_name == chr
    ]
    if (
        !nrow(exons)
    ) fail(
        "No exon reference rows for ",
        prefix
    )
    # Input interval and strand must belong to the supplied reference.
    observed <- unique(
        base[
            ,
            .(
                ensembl_gene_id,
                exon_start,
                exon_end,
                strand
            )
        ]
    )
    reference <- unique(
        exons[
            ,
            .(
                ensembl_gene_id,
                exon_start,
                exon_end,
                strand
            )
        ]
    )
    absent <- fsetdiff(
        observed,
        reference
    )
    report(
        absent,
        "unmatched_exon_intervals.tsv"
    )
    if (
        nrow(absent)
    ) fail(
        "Input exon intervals/strands absent from exon reference."
    )
    anchors <- exons[
        ,
        .(
            CHR = chromosome_name[
                1L
            ],
            gene_start = min(exon_start)
        ),
        by = ensembl_gene_id
    ]
    ################################################################################
    # Step 9. Construct biological masks
    ################################################################################
    log(
        "[INFO] Constructing five masks; no frequency filtering."
    )
    masks <- lapply(
        names(scores),
        function(s) {
            value <- inputs[[s]][[
                scores[[s]]
            ]]
            which(
                !is.na(value) & value >= thresholds[[s]]
            )
        }
    )
    names(masks) <- names(scores)
    masks$FUNC_ALL <- sort(
        unique(
            unlist(
                masks,
                use.names = FALSE
            )
        )
    )
    masks$ALL <- seq_len(
        nrow(base)
    )


    ################################################################################
    # Step 10. Summarize masks and check for empty masks
    ################################################################################

    summary_dt <- rbindlist(
        lapply(
            mask_names,
            function(s) {
                dt <- base[
                    masks[[s]]
                ]
                data.table(
                    mask = s,
                    n_pairs = nrow(dt),
                    n_variants = uniqueN(dt$ID),
                    n_genes = uniqueN(dt$ensembl_gene_id)
                )
            }
        )
    )
    summary_dt[
        ,
        proportion_of_ALL_pairs := n_pairs / nrow(base)
    ]
    summary_dt[
        ,
        proportion_of_ALL_variants := n_variants / uniqueN(base$ID)
    ]
    summary_dt[
        ,
        proportion_of_ALL_genes := n_genes / uniqueN(base$ensembl_gene_id)
    ]
    report(
        summary_dt,
        "mask_summary.tsv"
    )
    print(summary_dt)
    if (
        any(summary_dt$n_pairs == 0L)
    ) fail(
        "Empty mask(s): ",
        paste(
            summary_dt[
                n_pairs == 0L,
                mask
            ],
            collapse = ", "
        ),
        ". No chromosome outputs published (same policy as UKBB)."
    )
    ################################################################################
    # Step 11. Write REGENIE files and verify written pair membership
    ################################################################################
    for (
        s in mask_names
    ) {
        log(
            "[WRITE] ",
            s
        )
        dt <- merge(
            base[
                masks[[s]],
                .(
                    ID,
                    ensembl_gene_id,
                    POS
                )
            ],
            anchors,
            by = "ensembl_gene_id",
            all.x = TRUE,
            sort = FALSE
        )
        if (
            anyNA(dt$gene_start)
        ) fail(
            "Missing gene anchor in ",
            s
        )
        setorder(
            dt,
            ensembl_gene_id,
            POS,
            ID
        )
        annotation <- dt[
            ,
            .(
                ID,
                ensembl_gene_id,
                annotation = s
            )
        ]
        setlist <- dt[
            ,
            .(
                CHR = CHR[
                    1L
                ],
                gene_start = gene_start[
                    1L
                ],
                variants = paste(
                    ID,
                    collapse = ","
                )
            ),
            by = ensembl_gene_id
        ]
        directory <- file.path(
            stage,
            s
        )
        dir.create(directory)
        af <- file.path(
            directory,
            paste0(
                prefix,
                ".",
                s,
                ".annotation.txt"
            )
        )
        sf <- file.path(
            directory,
            paste0(
                prefix,
                ".",
                s,
                ".setlist.txt"
            )
        )
        mf <- file.path(
            directory,
            paste0(
                s,
                ".mask.def"
            )
        )
        fwrite(
            annotation,
            af,
            sep = "\t",
            col.names = FALSE,
            quote = FALSE
        )
        fwrite(
            setlist,
            sf,
            sep = "\t",
            col.names = FALSE,
            quote = FALSE
        )
        fwrite(
            data.table(
                mask = s,
                annotation = s
            ),
            mf,
            sep = "\t",
            col.names = FALSE,
            quote = FALSE
        )
        aa <- fread(
            af,
            header = FALSE,
            sep = "\t",
            colClasses = "character",
            col.names = c(
                "ID",
                "ensembl_gene_id",
                "annotation"
            )
        )
        ss <- fread(
            sf,
            header = FALSE,
            sep = "\t",
            col.names = c(
                "ensembl_gene_id",
                "CHR",
                "gene_start",
                "variants"
            )
        )
        mm <- fread(
            mf,
            header = FALSE,
            sep = "\t",
            colClasses = "character"
        )
        expanded <- ss[
            ,
            .(
                ID = unlist(
                    strsplit(
                        variants,
                        ",",
                        fixed = TRUE
                    )
                )
            ),
            by = ensembl_gene_id
        ][
            ,
            .(
                ID,
                ensembl_gene_id
            )
        ]
        expected <- dt[
            ,
            .(
                ID,
                ensembl_gene_id
            )
        ]
        if (
            nrow(aa) != nrow(expected) || any(aa$annotation != s) ||
            nrow(expanded) != nrow(expected) || anyDuplicated(ss$ensembl_gene_id) ||
            nrow(
                fsetdiff(
                    pairs(aa),
                    expected
                )
            ) || nrow(
                fsetdiff(
                    expected,
                    pairs(aa)
                )
            ) ||
            nrow(
                fsetdiff(
                    expanded,
                    expected
                )
            ) || nrow(
                fsetdiff(
                    expected,
                    expanded
                )
            ) ||
            !identical(
                as.character(
                    unlist(
                        mm,
                        use.names = FALSE
                    )
                ),
                c(
                    s,
                    s
                )
            )
        )
        fail(
            "Written annotation/set-list/mask-definition QC failed: ",
            s
        )
        if (
            !isTRUE(
                all.equal(
                    ss[
                        ,
                        .(
                            ensembl_gene_id,
                            CHR,
                            gene_start
                        )
                    ],
                    setlist[
                        ,
                        .(
                            ensembl_gene_id,
                            CHR,
                            gene_start
                        )
                    ],
                    check.attributes = FALSE
                )
            )
        ) fail(
            "Written anchor QC failed: ",
            s
        )
    }


    ################################################################################
    # Step 12. Save input provenance and session information
    ################################################################################

    # Record resolved paths, file sizes and modification times without a costly
    # second full read of each large input. Session details capture R/packages.
    all_paths <- c(
        paths,
        EXONS = exon_path
    )
    info <- file.info(all_paths)
    report(
        data.table(
            source = names(all_paths),
            path = normalizePath(all_paths),
            size_bytes = info$size,
            modified = as.character(info$mtime)
        ),
        "input_manifest.tsv"
    )
    writeLines(
        c(
            paste("Script version: 1.0.0"),
            paste(
                "Command:",
                paste(
                    commandArgs(),
                    collapse = " "
                )
            ),
            capture.output(
                sessionInfo()
            )
        ),
        file.path(
            stage,
            paste0(
                prefix,
                ".sessionInfo.txt"
            )
        )
    )


    ################################################################################
    # Step 13. Publish validated chromosome outputs and report completion
    ################################################################################

    writeLines(
        paste(
            "[PASS]",
            prefix,
            "completed",
            format(
                Sys.time(),
                tz = "UTC"
            )
        ),
        file.path(
            stage,
            paste0(
                prefix,
                ".SUCCESS"
            )
        )
    )
    if (
        file.exists(final) || !file.rename(
            stage,
            final
        )
    )
    fail(
        "Cannot publish chromosome directory; validated files remain at ",
        stage
    )
    published <- TRUE
    log(
        "[PASS] All five masks written and verified: ",
        final
    )
}


################################################################################
# Run the workflow and return a nonzero exit status on failure
################################################################################

tryCatch(
    main(),
    error = function(e) {
        message(
            "[ERROR] ",
            conditionMessage(e)
        )
        quit(
            save = "no",
            status = 1L
        )
    }
)
