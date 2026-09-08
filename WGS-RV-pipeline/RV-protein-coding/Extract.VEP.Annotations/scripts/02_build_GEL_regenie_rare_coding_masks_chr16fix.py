#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Construct GEL rare-coding masks for REGENIE
# CloudOS interactive session | UKBB-compatible annotation and mask rules
################################################################################

"""Build GEL rare protein-coding REGENIE inputs, preserving the supplied UKBB rules.

Author: Shelley
Adapted: 2026-09-07
Reference: supplied 02_Create.RareCoding.Regenie.input.masks(2).py
Dependencies: Python >=3.10, numpy, pandas, pyarrow.

Input: chromosome-level GEL VEP TSV (or gzip TSV), one annotation per row;
       Ensembl 114 GRCh38 gene-coordinate TSV used for UKBB.
Output: 3 files per mask (pLoF_only, pLoF_Dmis, Missense, Synonymous),
        one audit Parquet and one summary TSV: 14 files per chromosome.

Scientific rules preserved from the UKBB implementation:
* Keep protein_coding rows with MANE_SELECT present OR, on that same row,
  MANE_SELECT missing and CANONICAL == YES. This is not a gene-wide fallback.
* pLoF: LoF == HC; damaging: a MISSENSE_CSQS consequence with REVEL >=0.773
  or CADD_PHRED >=28.1, OR row-specific SpliceAI DSmax >=0.20, OR LoF == LC.
* Collapse to variant x gene by highest annotation priority. Broad Missense
  uses any retained missense/protein-altering consequence; Synonymous uses
  the final hierarchical class. Masks can overlap.
* REVEL/CADD: first valid ampersand-separated number, as in UKBB.
* No cohort AF, gnomAD AF, or VCF FILTER restriction is applied here.

GEL adaptation:
* Preserve original ID exactly in REGENIE outputs; do not add DRAGEN:.
* Use CHROM/POS/REF/ALT only for a separate internal matching/sorting key.
  This key standardizes chromosome naming; it does not left-align indels.
* Use the four gene/transcript-specific SpliceAI columns supplied in the TSV.
  DSmax uses available numeric scores; four missing scores remain missing.
* Use CADD_PHRED directly; no external CADD or SpliceAI files or tabix needed.
* Read in chunks, then retain selected rows in memory for global collapse.
  Chunking therefore does not place a fixed bound on total memory use.
* Build outputs in a temporary directory; validate written annotation/set-list
  pairs before publishing. Rerunning replaces this chromosome's named outputs.

Example (after placing this script in the output directory):
  python 02_build_GEL_regenie_rare_coding_masks.py --chrom 22
"""

# ==============================================================================
# Step 1. Load packages
# ==============================================================================

from __future__ import annotations

import argparse
import platform
import time
from datetime import datetime
import re
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

# ==============================================================================
# Step 2. Define BRaVa annotation criteria
# ==============================================================================

REVEL_CUTOFF = 0.773

CADD_CUTOFF = 28.1

SPLICEAI_CUTOFF = 0.20

# Exact consequence groups used by the BRaVa annotation pipeline.
# PLOF_CSQS is retained for reference, as in the UKBB script.
# Classification below uses LOFTEE HC directly; it does not test this set.

PLOF_CSQS = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
}


MISSENSE_CSQS = {
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
}

SYNONYMOUS_CSQS = {
    "stop_retained_variant",
    "synonymous_variant",
}

OTHER_CSQS = {
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
}

INFRAME_CSQS = {
    "inframe_deletion",
    "inframe_insertion",
}

MASK_NAMES = [
    "pLoF_only",
    "pLoF_Dmis",
    "Missense",
    "Synonymous",
]

ANNOTATION_PRIORITY = {
    "": 0,
    "synonymous": 1,
    "non_coding": 2,
    "other_missense_or_protein_altering": 3,
    "damaging_missense_or_protein_altering": 4,
    "pLoF": 5,
}

# ==============================================================================
# Step 3. Logging utilities
# ==============================================================================

def log(message: str) -> None:
    print(
        f"[INFO] {message}",
        file=sys.stderr,
        flush=True,
    )

def warning(message: str) -> None:
    print(
        f"[WARNING] {message}",
        file=sys.stderr,
        flush=True,
    )

def fail(message: str) -> None:
    raise RuntimeError(message)

# ==============================================================================
# Step 4. General helper functions
# ==============================================================================

def strip_ensembl_version(value) -> str | None:
    """
    Remove Ensembl version suffix.

    Example:
        ENSG00000123456.8
            ->
        ENSG00000123456
    """

    if pd.isna(value):
        return None

    value = str(value).strip()

    if value in {"", "."}:
        return None

    return value.split(".")[0]

def normalize_chromosome(chrom) -> str:
    """
    Normalize chromosome naming.

    Examples:
        chr10 -> 10
        10    -> 10
    """

    chrom = str(chrom).strip()

    return re.sub(
        r"^chr",
        "",
        chrom,
        flags=re.IGNORECASE,
    )

def get_variant_position(
    value: str,
) -> int:
    """
    Extract POS from:
        chrCHROM:POS:REF:ALT
    """

    parts = value.split(":")

    if len(parts) != 4:
        fail(
            f"Malformed normalized variant ID: {value}"
        )

    return int(
        parts[1]
    )

def consequence_set(
    value,
) -> set[str]:
    """
    Convert '&'-separated VEP consequence field to a set.
    """

    if pd.isna(value):
        return set()

    return {
        x.strip()
        for x in str(value).split("&")
        if x.strip() not in {"", "."}
    }

# ==============================================================================
# Step 5. Parse REVEL and CADD
# ==============================================================================

def parse_revel(
    value,
) -> float:
    """
    Follow the BRaVa REVEL parsing implementation.

    REVEL may contain '&'-separated values.

    BRaVa retains the first non-missing value.
    """

    if pd.isna(value):
        return np.nan

    values = str(value).split("&")

    for x in values:

        x = x.strip()

        if x in {"", "."}:
            continue

        try:
            return float(x)

        except ValueError:
            continue

    return np.nan

def parse_cadd(
    value,
) -> float:
    """
    Parse CADD PHRED annotation.

    If an '&'-separated field is encountered, retain the first valid
    non-missing numeric value.
    """

    if pd.isna(value):
        return np.nan

    values = str(value).split("&")

    for x in values:

        x = x.strip()

        if x in {"", "."}:
            continue

        try:
            return float(x)

        except ValueError:
            continue

    return np.nan

# ==============================================================================
# Step 6. Read Ensembl 114 gene-coordinate reference
# ==============================================================================

def read_gene_coordinates(
    path: Path,
    chromosome: str,
) -> pd.DataFrame:
    """
    Read Ensembl 114 GRCh38 gene-coordinate reference.

    Required columns:
        ensembl_gene_id
        hgnc_symbol
        chromosome_name
        start_position
        end_position

    Ensembl gene ID is the final REGENIE set identifier.

    HGNC symbol is used only as an intermediate mapping key for raw
    Illumina SpliceAI annotations.
    """

    log(
        f"Reading Ensembl gene coordinates: {path}"
    )

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
    )

    required = {
        "ensembl_gene_id",
        "hgnc_symbol",
        "chromosome_name",
        "start_position",
        "end_position",
    }

    missing = (
        required
        - set(df.columns)
    )

    if missing:

        fail(
            "Gene-coordinate file is missing required columns: "
            + ", ".join(
                sorted(missing)
            )
        )

    df[
        "ensembl_gene_id"
    ] = (
        df[
            "ensembl_gene_id"
        ]
        .map(
            strip_ensembl_version
        )
    )

    df[
        "chromosome_name"
    ] = (
        df[
            "chromosome_name"
        ]
        .map(
            normalize_chromosome
        )
    )

    df[
        "start_position"
    ] = pd.to_numeric(
        df[
            "start_position"
        ],
        errors="raise",
    )

    df[
        "end_position"
    ] = pd.to_numeric(
        df[
            "end_position"
        ],
        errors="raise",
    )

    df[
        "hgnc_symbol"
    ] = (
        df[
            "hgnc_symbol"
        ]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    df = (
        df.loc[
            df[
                "chromosome_name"
            ]
            == chromosome
        ]
        .copy()
    )

    if df.empty:

        fail(
            f"No genes found for chromosome {chromosome} "
            f"in {path}"
        )

    # --------------------------------------------------------------------------
    # Ensembl gene IDs must be unique.
    # --------------------------------------------------------------------------

    duplicated = (
        df[
            "ensembl_gene_id"
        ]
        .duplicated(
            keep=False
        )
    )

    if duplicated.any():

        duplicated_ids = (
            df.loc[
                duplicated,
                "ensembl_gene_id",
            ]
            .dropna()
            .drop_duplicates()
            .head(20)
            .tolist()
        )

        fail(
            "Duplicated Ensembl gene IDs remain in the gene-coordinate "
            "reference. Examples: "
            + ", ".join(
                duplicated_ids
            )
        )

    df = (
        df.sort_values(
            [
                "start_position",
                "ensembl_gene_id",
            ],
            kind="mergesort",
        )
        .reset_index(
            drop=True
        )
    )

    log(
        f"Gene-coordinate reference: "
        f"{len(df):,} genes on chr{chromosome}"
    )

    return df

# ==============================================================================
# Step 7. BRaVa annotation hierarchy
# ==============================================================================

def get_brava_annotation(
    lof,
    revel,
    cadd,
    csq,
    max_ds,
) -> str:
    """
    Apply the BRaVa annotation hierarchy.

    Priority:
        1. LOFTEE HC
        2. missense/protein-altering + REVEL/CADD threshold
        3. SpliceAI DSmax >= 0.20
        4. LOFTEE LC
        5. other missense/protein-altering
        6. non-coding
        7. synonymous
        8. unclassified
    """

    consequences = (
        consequence_set(
            csq
        )
    )

    missense_variant = bool(
        consequences
        & MISSENSE_CSQS
    )

    synonymous_variant = bool(
        consequences
        & SYNONYMOUS_CSQS
    )

    other_variant = bool(
        consequences
        & OTHER_CSQS
    )

    inframe_variant = bool(
        consequences
        & INFRAME_CSQS
    )

    # --------------------------------------------------------------------------
    # 1. High-confidence pLoF
    # --------------------------------------------------------------------------

    if lof == "HC":

        return "pLoF"

    # --------------------------------------------------------------------------
    # 2. Damaging missense/protein-altering by REVEL or CADD
    # --------------------------------------------------------------------------

    elif (
        missense_variant

        and

        (
            (
                pd.notna(
                    revel
                )

                and

                revel
                >= REVEL_CUTOFF
            )

            or

            (
                pd.notna(
                    cadd
                )

                and

                cadd
                >= CADD_CUTOFF
            )
        )
    ):

        return (
            "damaging_missense_or_protein_altering"
        )

    # --------------------------------------------------------------------------
    # 3. Damaging by gene-specific SpliceAI
    # --------------------------------------------------------------------------

    elif (
        pd.notna(
            max_ds
        )

        and

        max_ds
        >= SPLICEAI_CUTOFF
    ):

        return (
            "damaging_missense_or_protein_altering"
        )

    # --------------------------------------------------------------------------
    # 4. LOFTEE LC
    # --------------------------------------------------------------------------

    elif lof == "LC":

        return (
            "damaging_missense_or_protein_altering"
        )

    # --------------------------------------------------------------------------
    # 5. Other missense/protein-altering
    # --------------------------------------------------------------------------

    elif (
        missense_variant
        or
        inframe_variant
    ):

        return (
            "other_missense_or_protein_altering"
        )

    # --------------------------------------------------------------------------
    # 6. Non-coding
    # --------------------------------------------------------------------------

    elif other_variant:

        return "non_coding"

    # --------------------------------------------------------------------------
    # 7. Synonymous
    # --------------------------------------------------------------------------

    elif synonymous_variant:

        return "synonymous"

    else:

        return ""

# ==============================================================================
# Step 8. Assign BRaVa annotations
# ==============================================================================

def annotate_vep_rows(
    df: pd.DataFrame,
) -> pd.DataFrame:

    log(
        "Assigning BRaVa annotations ..."
    )

    df = df.copy()

    df[
        "BRAVA_ANNOTATION"
    ] = [
        get_brava_annotation(
            lof=lof,
            revel=revel,
            cadd=cadd,
            csq=csq,
            max_ds=max_ds,
        )
        for (
            lof,
            revel,
            cadd,
            csq,
            max_ds,
        )
        in zip(
            df[
                "LOF"
            ],
            df[
                "REVEL_SCORE"
            ],
            df[
                "CADD_PHRED"
            ],
            df[
                "CSQ"
            ],
            df[
                "max_DS"
            ],
        )
    ]

    log(
        "BRaVa annotation counts:\n"
        + df[
            "BRAVA_ANNOTATION"
        ]
        .value_counts(
            dropna=False
        )
        .to_string()
    )

    return df

# ==============================================================================
# Step 9. Collapse to one variant x Ensembl gene pair
# ==============================================================================

def collapse_variant_gene(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Collapse residual transcript-level duplicates to one row per:

        variant x Ensembl gene

    BRaVa annotation:
        retain the highest-priority classification.

    Additional consequence flags:
        preserve whether ANY retained annotation contains a BRaVa
        missense or synonymous consequence.
    """

    log(
        "Collapsing to one row per variant x Ensembl gene ..."
    )

    df = df.copy()

    # --------------------------------------------------------------------------
    # Preserve consequence membership.
    # --------------------------------------------------------------------------

    df[
        "IS_BRAVA_MISSENSE_CSQ"
    ] = (
        df[
            "CSQ"
        ]
        .map(
            lambda x:
                bool(
                    consequence_set(x)
                    & MISSENSE_CSQS
                )
        )
    )

    df[
        "IS_BRAVA_SYNONYMOUS_CSQ"
    ] = (
        df[
            "CSQ"
        ]
        .map(
            lambda x:
                bool(
                    consequence_set(x)
                    & SYNONYMOUS_CSQS
                )
        )
    )

    # --------------------------------------------------------------------------
    # Assign annotation priority.
    # --------------------------------------------------------------------------

    df[
        "BRAVA_PRIORITY"
    ] = (
        df[
            "BRAVA_ANNOTATION"
        ]
        .map(
            ANNOTATION_PRIORITY
        )
    )

    if df[
        "BRAVA_PRIORITY"
    ].isna().any():

        fail(
            "Unknown BRaVa annotation encountered during collapse."
        )

    # --------------------------------------------------------------------------
    # Highest-priority representative row.
    # --------------------------------------------------------------------------

    df = df.sort_values(
        [
            "VARIANT_ID_CHR",
            "ENSEMBL_GENE_ID",
            "BRAVA_PRIORITY",
        ],
        ascending=[
            True,
            True,
            False,
        ],
        kind="mergesort",
    )

    representative = (
        df.drop_duplicates(
            subset=[
                "VARIANT_ID_CHR",
                "ENSEMBL_GENE_ID",
            ],
            keep="first",
        )
        .copy()
    )

    # --------------------------------------------------------------------------
    # Preserve ANY consequence flags across residual rows.
    # --------------------------------------------------------------------------

    flags = (
        df.groupby(
            [
                "VARIANT_ID_CHR",
                "ENSEMBL_GENE_ID",
            ],
            as_index=False,
            sort=False,
        )
        .agg(
            ANY_BRAVA_MISSENSE_CSQ=(
                "IS_BRAVA_MISSENSE_CSQ",
                "max",
            ),
            ANY_BRAVA_SYNONYMOUS_CSQ=(
                "IS_BRAVA_SYNONYMOUS_CSQ",
                "max",
            ),
        )
    )

    representative = (
        representative.drop(
            columns=[
                "IS_BRAVA_MISSENSE_CSQ",
                "IS_BRAVA_SYNONYMOUS_CSQ",
            ]
        )
        .merge(
            flags,
            on=[
                "VARIANT_ID_CHR",
                "ENSEMBL_GENE_ID",
            ],
            how="left",
            validate="one_to_one",
        )
    )

    if representative.duplicated(
        subset=[
            "VARIANT_ID_CHR",
            "ENSEMBL_GENE_ID",
        ]
    ).any():

        fail(
            "Duplicate variant-gene pairs remain after collapse."
        )

    log(
        f"Variant-gene pairs after collapse: "
        f"{len(representative):,}"
    )

    return representative

# ==============================================================================
# Step 10. Add Ensembl 114 gene coordinates
# ==============================================================================

def add_gene_coordinates(
    df: pd.DataFrame,
    coordinates: pd.DataFrame,
    outdir: Path,
    chromosome: str,
) -> pd.DataFrame:

    # --------------------------------------------------------------------------
    # Chr16-only exclusion agreed after manual review (2026-09-08).
    # ENSG00000310590 is absent from the shared UKBB/GEL coordinate reference.
    # Keep the original strict coordinate checks for every other gene.
    # --------------------------------------------------------------------------

    if chromosome == "16":

        excluded = df["ENSEMBL_GENE_ID"].eq("ENSG00000310590")

        if excluded.any():

            warning(
                f"Chr16 reference-based exclusion: ENSG00000310590; "
                f"{excluded.sum():,} variant-gene pairs removed."
            )

            log(
                "Excluded annotation classes:\n"
                + df.loc[excluded, "BRAVA_ANNOTATION"]
                .value_counts(dropna=False)
                .to_string()
            )

            df = df.loc[~excluded].copy()

            log(
                f"Variant-gene pairs retained after chr16 exclusion: {len(df):,}"
            )

    coordinate_table = (
        coordinates[
            [
                "ensembl_gene_id",
                "chromosome_name",
                "start_position",
            ]
        ]
        .rename(
            columns={
                "ensembl_gene_id":
                    "ENSEMBL_GENE_ID",

                "chromosome_name":
                    "GENE_CHR",

                "start_position":
                    "GENE_POS",
            }
        )
    )

    df = df.merge(
        coordinate_table,
        on="ENSEMBL_GENE_ID",
        how="left",
        validate="many_to_one",
    )

    missing = (
        df[
            "GENE_POS"
        ]
        .isna()
    )

    if missing.any():

        missing_file = (
            outdir
            / (
                f"chr{chromosome}."
                f"missing_gene_coordinates.tsv.gz"
            )
        )

        df.loc[
            missing
        ].to_csv(
            missing_file,
            sep="\t",
            index=False,
            compression="gzip",
        )

        fail(
            f"{missing.sum():,} variant-gene rows could not be mapped "
            f"to the Ensembl 114 coordinate reference. "
            f"See: {missing_file}"
        )

    df[
        "GENE_POS"
    ] = (
        df[
            "GENE_POS"
        ]
        .astype(int)
    )

    unexpected_chr = (
        df[
            "GENE_CHR"
        ]
        .astype(str)
        != chromosome
    )

    if unexpected_chr.any():

        fail(
            "Gene-coordinate chromosome mismatch detected."
        )

    return df

# ==============================================================================
# Step 11. Define the four rare-coding masks
# ==============================================================================

def define_masks(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """
    pLoF_only:
        BRaVa pLoF

    pLoF_Dmis:
        BRaVa pLoF
        +
        BRaVa damaging_missense_or_protein_altering

    Missense:
        Any selected VEP consequence belonging to BRaVa MISSENSE_CSQS.

        No REVEL/CADD/SpliceAI threshold is required for this broad mask.

    Synonymous:
        Final hierarchical BRaVa annotation == synonymous.

        Therefore synonymous variants promoted to a damaging class,
        for example by SpliceAI >= 0.20, are excluded.
    """

    df = df.copy()

    df[
        "MASK_pLoF_only"
    ] = (
        df[
            "BRAVA_ANNOTATION"
        ]
        == "pLoF"
    )

    df[
        "MASK_pLoF_Dmis"
    ] = (
        df[
            "BRAVA_ANNOTATION"
        ]
        .isin(
            {
                "pLoF",
                "damaging_missense_or_protein_altering",
            }
        )
    )

    df[
        "MASK_Missense"
    ] = (
        df[
            "ANY_BRAVA_MISSENSE_CSQ"
        ]
        .astype(bool)
    )

    df[
        "MASK_Synonymous"
    ] = (
        df[
            "BRAVA_ANNOTATION"
        ]
        == "synonymous"
    )

    return df

# ==============================================================================
# Step 12. Define GEL paths and input-column mapping
# ==============================================================================

VEP_DIRECTORY = Path('/home/vscode/session_data/filesystems/chromosome_level_VEP')
GENE_COORDINATES = Path('/home/vscode/session_data/mounted-data-readonly/Ensembl114_GRCh38_gene_coordinates.tsv')
OUTPUT_DIRECTORY = Path('/home/vscode/session_data/Rare-Coding-Masks')
SPLICEAI_COLUMNS = [f'SpliceAI_pred_DS_{event}' for event in ('AG', 'AL', 'DG', 'DL')]
COLUMN_MAP = {
    'ID': 'SNP_ID', 'Gene': 'GENE', 'LoF': 'LOF', 'REVEL': 'REVEL_SCORE',
    'Consequence': 'CSQ', 'Feature': 'TRANSCRIPT',
}
REQUIRED_COLUMNS = [
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'Gene', 'LoF', 'REVEL',
    'CADD_PHRED', 'Consequence', 'Feature', 'MANE_SELECT', 'CANONICAL',
    'BIOTYPE', *SPLICEAI_COLUMNS,
]

# ==============================================================================
# Step 13. Read and validate GEL VEP annotations in chunks
# ==============================================================================

def read_gel_vep(path: Path, chromosome: str, chunksize: int) -> pd.DataFrame:
    """Read GEL TSVs with explicit missing-value handling and row-specific scores."""

    # --------------------------------------------------------------------------
    # 1. Verify required GEL columns and configure TSV parsing.
    # --------------------------------------------------------------------------

    header = pd.read_csv(path, sep='\t', nrows=0, encoding='utf-8-sig')

    missing = set(REQUIRED_COLUMNS) - set(header.columns)

    if missing:
        fail('GEL VEP input is missing columns: ' + ', '.join(sorted(missing)))

    log(f'Reading GEL VEP: {path}')

    selected, total, retained, missing_genes = [], 0, 0, 0

    reader = pd.read_csv(
        path, sep='\t', usecols=REQUIRED_COLUMNS, dtype=str,
        keep_default_na=False, na_values=['', '.', 'NA', 'NaN', 'nan', 'NAN', '<NA>'],
        encoding='utf-8-sig', chunksize=chunksize,
    )

    for number, chunk in enumerate(reader, 1):
        total += len(chunk)
        # --------------------------------------------------------------------------
        # Select protein-coding MANE/canonical rows, as in UKBB.
        # --------------------------------------------------------------------------

        keep = chunk['BIOTYPE'].eq('protein_coding') & (
            chunk['MANE_SELECT'].notna()
            | (chunk['MANE_SELECT'].isna() & chunk['CANONICAL'].eq('YES'))
        )
        chunk = chunk.loc[keep].copy().rename(columns=COLUMN_MAP)
        retained += len(chunk)
        if chunk.empty:
            log(f'Chunk {number}: {total:,} rows read; {retained:,} selected')
            continue
        # --------------------------------------------------------------------------
        # Normalize Ensembl IDs and exclude rows without a gene.
        # --------------------------------------------------------------------------

        chunk['ENSEMBL_GENE_ID'] = chunk['GENE'].map(strip_ensembl_version)
        no_gene = chunk['ENSEMBL_GENE_ID'].isna()
        missing_genes += int(no_gene.sum())
        chunk = chunk.loc[~no_gene].copy()
        if chunk.empty:
            continue
        # Keep SNP_ID untouched. Internal keys are constructed separately.
        bad_id = chunk['SNP_ID'].isna() | chunk['SNP_ID'].str.contains(r'[\s,]', na=False)
        if bad_id.any():
            fail(f'Chunk {number}: missing or whitespace/comma-containing variant IDs.')
        bad_gene = ~chunk['ENSEMBL_GENE_ID'].str.fullmatch(r'ENSG[0-9]+', na=False)
        if bad_gene.any():
            fail(f'Chunk {number}: invalid Ensembl gene IDs: '
                 + str(chunk.loc[bad_gene, 'ENSEMBL_GENE_ID'].head().tolist()))
        chunk['CHROM'] = chunk['CHROM'].map(normalize_chromosome)
        if not chunk['CHROM'].eq(chromosome).all():
            fail(f'Chunk {number}: selected VEP rows do not all belong to chr{chromosome}.')
        if not chunk['POS'].str.fullmatch(r'[0-9]+', na=False).all():
            fail(f'Chunk {number}: missing or noninteger POS values.')
        chunk['POS'] = pd.to_numeric(chunk['POS'], errors='raise').astype('int64')
        if not chunk['POS'].gt(0).all():
            fail(f'Chunk {number}: nonpositive positions.')
        for allele in ('REF', 'ALT'):
            if (chunk[allele].isna() | chunk[allele].str.contains(r'[\s,:]', na=False)).any():
                fail(f'Chunk {number}: missing/malformed {allele} or unsplit multiallelic rows.')
        chunk['VARIANT_ID_CHR'] = (
            'chr' + chunk['CHROM'] + ':' + chunk['POS'].astype(str)
            + ':' + chunk['REF'] + ':' + chunk['ALT']
        )
        # --------------------------------------------------------------------------
        # Parse embedded REVEL, CADD, and gene-specific SpliceAI scores.
        # --------------------------------------------------------------------------

        chunk['REVEL_SCORE'] = chunk['REVEL_SCORE'].map(parse_revel)
        chunk['CADD_PHRED'] = chunk['CADD_PHRED'].map(parse_cadd)
        for column in SPLICEAI_COLUMNS:
            raw = chunk[column]
            numeric = pd.to_numeric(raw, errors='coerce')
            invalid = raw.notna() & (numeric.isna() | ~np.isfinite(numeric) | ~numeric.between(0, 1))
            if invalid.any():
                fail(f'Chunk {number}: {column} requires scalar scores in [0,1] or missing; '
                     f'examples: {raw.loc[invalid].head().tolist()}')
            chunk[column] = numeric
        chunk['SPLICEAI_N_SCORES'] = chunk[SPLICEAI_COLUMNS].notna().sum(axis=1)
        chunk['max_DS'] = chunk[SPLICEAI_COLUMNS].max(axis=1, skipna=True)
        selected.append(chunk)
        log(f'Chunk {number}: {total:,} rows read; {retained:,} selected')

    if missing_genes:
        warning(f'Excluded {missing_genes:,} selected rows without an Ensembl gene ID.')

    if not selected:
        fail('No usable VEP rows remain after protein-coding/MANE-canonical selection.')

    df = pd.concat(selected, ignore_index=True)
    # Prevent losing or conflating genotype IDs when collapsing annotation rows.

    # --------------------------------------------------------------------------
    # Check the one-to-one relationship between genotype IDs and allele keys.
    # --------------------------------------------------------------------------

    id_keys = df[['SNP_ID', 'VARIANT_ID_CHR']].drop_duplicates()

    if id_keys['SNP_ID'].duplicated().any():
        fail('One genotype ID maps to multiple coordinate/allele keys in selected rows.')

    if id_keys['VARIANT_ID_CHR'].duplicated().any():
        fail('Multiple genotype IDs map to one coordinate/allele key in selected rows.')

    # --------------------------------------------------------------------------
    # Report retained rows, score availability, and threshold counts.
    # --------------------------------------------------------------------------

    log(f'Usable selected rows: {len(df):,} / {total:,}')

    for column in ('REVEL_SCORE', 'CADD_PHRED', 'max_DS'):
        log(f'{column} nonmissing: {df[column].notna().sum():,} / {len(df):,}')

    for column, cutoff in (('REVEL_SCORE', REVEL_CUTOFF),
                           ('CADD_PHRED', CADD_CUTOFF), ('max_DS', SPLICEAI_CUTOFF)):
        log(f'{column} >= {cutoff}: {df[column].ge(cutoff).sum():,} selected annotation rows')

    partial = df['SPLICEAI_N_SCORES'].between(1, 3).sum()

    log(f'SpliceAI rows with 1-3 available scores: {partial:,}; max uses available scores.')

    if df['max_DS'].isna().all():
        warning('All selected SpliceAI scores are missing; no SpliceAI-based promotions are possible.')

    return df

# ==============================================================================
# Step 14. Validate written REGENIE files
# ==============================================================================

def validate_written_mask(annotation_path: Path, setlist_path: Path,
                          maskdef_path: Path, expected: set, mask_name: str,
                          chromosome: str, gene_positions: dict) -> None:
    """Validate the actual files, including empty masks and duplicate records."""

    annotation_pairs, setlist_pairs, seen_genes = set(), set(), set()

    # --------------------------------------------------------------------------
    # Read back annotation rows and check format and uniqueness.
    # --------------------------------------------------------------------------

    with annotation_path.open() as handle:
        for line in handle:
            parts = line.rstrip('\n').split(' ')
            if len(parts) != 3 or parts[2] != mask_name:
                fail(f'Malformed annotation line in {annotation_path}')
            pair = (parts[0], parts[1])
            if pair in annotation_pairs:
                fail(f'Duplicate written annotation pair: {pair}')
            annotation_pairs.add(pair)

    # --------------------------------------------------------------------------
    # Read back set-list rows and check coordinates and variant membership.
    # --------------------------------------------------------------------------

    with setlist_path.open() as handle:
        for line in handle:
            parts = line.rstrip('\n').split(' ')
            if len(parts) != 4:
                fail(f'Malformed set-list line in {setlist_path}')
            gene, chrom, position, variants = parts
            if gene in seen_genes or chrom != chromosome or int(position) != gene_positions[gene]:
                fail(f'Invalid gene/chromosome/position in set-list for {gene}')
            seen_genes.add(gene)
            for variant in variants.split(','):
                pair = (variant, gene)
                if not variant or pair in setlist_pairs:
                    fail(f'Empty/duplicate written set-list variant for {gene}')
                setlist_pairs.add(pair)

    # --------------------------------------------------------------------------
    # Require exact agreement with the expected variant-gene pairs.
    # --------------------------------------------------------------------------

    if annotation_pairs != expected or setlist_pairs != expected:
        fail(f'Written annotation/set-list pairs differ from expected pairs for {mask_name}')

    if maskdef_path.read_text() != f'{mask_name} {mask_name}\n':
        fail(f'Invalid mask definition for {mask_name}')

# ==============================================================================
# Step 15. Write GEL REGENIE files and mask summaries
# ==============================================================================

def write_gel_mask(master: pd.DataFrame, mask_name: str, chromosome: str,
                   staging: Path, outdir: Path) -> dict:
    mask = master.loc[master[f'MASK_{mask_name}']].copy()

    if mask.duplicated(['SNP_ID', 'ENSEMBL_GENE_ID']).any():
        fail(f'Duplicate variant-gene pairs in {mask_name}')

    mask['VARIANT_POS'] = mask['VARIANT_ID_CHR'].map(get_variant_position)

    mask = mask.sort_values(['GENE_POS', 'ENSEMBL_GENE_ID', 'VARIANT_POS', 'SNP_ID'], kind='mergesort')

    prefix = f'chr{chromosome}.{mask_name}'

    annotation = staging / f'{prefix}.annotation.txt'

    setlist = staging / f'{prefix}.setlist.txt'

    maskdef = staging / f'{prefix}.maskdef.txt'

    # --------------------------------------------------------------------------
    # Write annotation file using the original GEL genotype IDs.
    # --------------------------------------------------------------------------

    with annotation.open('w') as handle:
        for variant, gene in mask[['SNP_ID', 'ENSEMBL_GENE_ID']].itertuples(index=False, name=None):
            handle.write(f'{variant} {gene} {mask_name}\n')

    # --------------------------------------------------------------------------
    # Write one set-list record per Ensembl gene.
    # --------------------------------------------------------------------------

    gene_positions = {}

    with setlist.open('w') as handle:
        for gene, group in mask.groupby('ENSEMBL_GENE_ID', sort=False):
            if group['GENE_CHR'].nunique() != 1 or group['GENE_POS'].nunique() != 1:
                fail(f'Conflicting gene coordinates for {gene}')
            position = int(group['GENE_POS'].iloc[0])
            gene_positions[gene] = position
            variants = ','.join(group['SNP_ID'])
            handle.write(f'{gene} {chromosome} {position} {variants}\n')

    # --------------------------------------------------------------------------
    # Write the mask definition and validate all three output files.
    # --------------------------------------------------------------------------

    maskdef.write_text(f'{mask_name} {mask_name}\n')

    expected = set(mask[['SNP_ID', 'ENSEMBL_GENE_ID']].itertuples(index=False, name=None))

    validate_written_mask(annotation, setlist, maskdef, expected, mask_name, chromosome, gene_positions)

    log(f'{mask_name}: {len(mask):,} pairs | {mask.SNP_ID.nunique():,} variants | '
        f'{mask.ENSEMBL_GENE_ID.nunique():,} genes; written-file QC passed')

    return {
        'mask': mask_name, 'variant_gene_pairs': len(mask),
        'unique_variants': mask.SNP_ID.nunique(), 'genes': mask.ENSEMBL_GENE_ID.nunique(),
        'annotation_file': str(outdir / annotation.name),
        'setlist_file': str(outdir / setlist.name),
        'maskdef_file': str(outdir / maskdef.name),
    }

# ==============================================================================
# Step 16. Main workflow
# ==============================================================================

def main() -> None:
    # --------------------------------------------------------------------------
    # 1. Parse arguments and resolve GEL input/output paths.
    # --------------------------------------------------------------------------

    started = time.monotonic()

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--chrom', required=True, help='Autosome 1-22; chr prefix accepted.')

    parser.add_argument('--vep', type=Path, help='Override the default GEL chromosome VEP TSV path.')

    parser.add_argument('--gene-coordinates', type=Path, default=GENE_COORDINATES)

    parser.add_argument('--outdir', type=Path, default=OUTPUT_DIRECTORY)

    parser.add_argument('--chunksize', type=int, default=1_000_000)

    args = parser.parse_args()

    chromosome = normalize_chromosome(args.chrom)

    if chromosome not in {str(i) for i in range(1, 23)}:
        parser.error('--chrom must be an autosome from 1 to 22.')

    if args.chunksize < 1:
        parser.error('--chunksize must be positive.')

    vep_path = args.vep or VEP_DIRECTORY / f'GEL.VEP.chr{chromosome}.tsv'

    # --------------------------------------------------------------------------
    # 2. Pre-flight checks: input files and Parquet dependency.
    # --------------------------------------------------------------------------

    for path in (vep_path, args.gene_coordinates):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f'Required input missing or empty: {path}')

    try:
        import pyarrow  # noqa: F401 -- fail early before chromosome-scale work
    except ImportError:
        fail('pyarrow is required to write the audit Parquet. Install it in the active environment.')

    args.outdir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------------------------
    # 3. Report runtime environment and inputs.
    # --------------------------------------------------------------------------

    log(f'GEL rare-coding mask construction: chr{chromosome}')

    log(f'Start: {datetime.now().astimezone().isoformat(timespec="seconds")}')

    log(f'Host: {platform.node()} | Python: {platform.python_version()} | {sys.executable}')

    log(f'numpy: {np.__version__} | pandas: {pd.__version__} | pyarrow: {pyarrow.__version__}')

    log(f'VEP input: {vep_path}')

    log(f'Gene coordinates: {args.gene_coordinates}')

    log(f'Output directory: {args.outdir} | Chunk size: {args.chunksize:,}')

    log('Original genotype IDs retained; embedded SpliceAI/CADD; no AF or FILTER restrictions.')

    # --------------------------------------------------------------------------
    # 4. Load and validate the Ensembl gene-coordinate reference.
    # --------------------------------------------------------------------------

    coordinates = read_gene_coordinates(args.gene_coordinates, chromosome)

    if coordinates['ensembl_gene_id'].isna().any():
        fail('Missing Ensembl IDs in the chromosome coordinate reference.')

    for field in ('start_position', 'end_position'):
        values = coordinates[field]
        if not (np.isfinite(values) & values.gt(0) & values.eq(np.floor(values))).all():
            fail(f'Gene reference {field} must contain positive integer coordinates.')

    if coordinates['end_position'].lt(coordinates['start_position']).any():
        fail('Gene reference contains end_position < start_position.')

    # --------------------------------------------------------------------------
    # 5. Read GEL annotations, classify, and collapse variant-gene pairs.
    # --------------------------------------------------------------------------

    vep = read_gel_vep(vep_path, chromosome, args.chunksize)

    master = collapse_variant_gene(annotate_vep_rows(vep))

    del vep

    # --------------------------------------------------------------------------
    # 6. Add coordinates and define the four masks.
    # --------------------------------------------------------------------------

    master = add_gene_coordinates(master, coordinates, args.outdir, chromosome)

    master = define_masks(master)

    master['VARIANT_ID_REGENIE'] = master['SNP_ID']

    # --------------------------------------------------------------------------
    # 7. Check the pLoF mask subset relationship.
    # --------------------------------------------------------------------------

    if (master['MASK_pLoF_only'] & ~master['MASK_pLoF_Dmis']).any():
        fail('QC failure: pLoF_only is not a subset of pLoF_Dmis.')

    log('QC passed: pLoF_only is a subset of pLoF_Dmis.')

    # --------------------------------------------------------------------------
    # 8. Stage the audit table, mask files, and summary; validate before publishing.
    # --------------------------------------------------------------------------

    with tempfile.TemporaryDirectory(prefix=f'.chr{chromosome}.masks.', dir=args.outdir) as temp:
        staging = Path(temp)
        master.to_parquet(staging / f'chr{chromosome}.annotation.audit.parquet', index=False,
                          engine='pyarrow', compression='snappy')
        summaries = [write_gel_mask(master, name, chromosome, staging, args.outdir) for name in MASK_NAMES]
        summary = pd.DataFrame(summaries)
        summary_name = f'chr{chromosome}.mask_summary.tsv'
        summary.to_csv(staging / summary_name, sep='\t', index=False)
        # Per-file atomic replacement; the summary is published last.
        for path in sorted(staging.iterdir()):
            if path.name != summary_name:
                path.replace(args.outdir / path.name)
        (staging / summary_name).replace(args.outdir / summary_name)

    # --------------------------------------------------------------------------
    # 9. Print the final mask summary, output sizes, and runtime.
    # --------------------------------------------------------------------------

    log('\n' + summary[['mask', 'variant_gene_pairs', 'unique_variants', 'genes']].to_string(index=False))

    log('All implemented QC checks passed.')

    for path in sorted(args.outdir.glob(f'chr{chromosome}.*')):
        if path.name.endswith(('.annotation.txt', '.setlist.txt', '.maskdef.txt',
                               '.annotation.audit.parquet', '.mask_summary.tsv')):
            log(f'Generated: {path.name} | {path.stat().st_size:,} bytes')

    log(f'End: {datetime.now().astimezone().isoformat(timespec="seconds")}')

    log(f'Elapsed: {time.monotonic() - started:.1f} seconds')

    log(f'Completed: 14 files in {args.outdir}. No AF filtering has been applied.')

# ==============================================================================
# Step 17. Run
# ==============================================================================

if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, ValueError, OSError) as error:
        print(f'[ERROR] {error}', file=sys.stderr, flush=True)
        sys.exit(1)
