nextflow.enable.dsl = 2

/*
===============================================================================
Genomics England CloudOS GWAS / REGENIE

Task 1: Shared genotype QC

Purpose
-------
Generate chromosome-specific high-quality common SNP genotype datasets
from Genomics England AggV3 WGS data for downstream REGENIE Step 1.

Shared QC
---------
1. Restrict to the analysis sample set
2. Restrict to AggV3 variants passing basic site QC
   (FILTER = PASS or LowMLSQ; supplied as chromosome-specific ID lists)
3. Retain A/C/G/T SNPs only
4. MAF >= 0.01
5. Variant missingness <= 0.01
6. Hardy-Weinberg equilibrium P >= 1e-15

LD pruning is intentionally NOT performed here. The resulting chromosome-
specific HQ common genotype datasets are intended to be reused by a separate
LD-pruning pipeline where pruning parameters can be tested iteratively.

Author: Shelley
===============================================================================
*/


/*
-------------------------------------------------------------------------------
Pipeline parameters
-------------------------------------------------------------------------------
*/

params.pgen_root        = null
params.variant_list_dir = null
params.keep_file        = null
params.outdir           = 'results'

params.maf              = 0.01
params.geno             = 0.01
params.hwe              = 1e-15

/*
The container image is defined in nextflow.config:

params.plink2_image =
    'ghcr.io/shiyuzhang0522/gel-gwas-regenie-plink2:alpha7-20260808'
*/


/*
-------------------------------------------------------------------------------
Process: chromosome-specific shared QC
-------------------------------------------------------------------------------
*/

process SHARED_QC {

    tag "chr${chr}"

    container params.plink2_image

    cpus 4
    memory '16 GB'

    publishDir "${params.outdir}/chr${chr}", mode: 'copy'


    /*
    Each task receives:
      - one chromosome number
      - one AggV3 PGEN/PVAR/PSAM fileset
      - one chromosome-specific PASS/LowMLSQ variant-ID list
      - the common sample keep file

    Files are staged with fixed local names so PLINK2 can use:
        --pfile dragen
    */

    input:

    tuple val(chr),
          path(pgen,         name: 'dragen.pgen'),
          path(pvar,         name: 'dragen.pvar'),
          path(psam,         name: 'dragen.psam'),
          path(variant_list, name: 'basic_QC.variant_ids.txt')

    path keep_file, name: 'GEL_CM_REGENIE.keep.txt'


    /*
    Output:
      - filtered chromosome-specific PGEN fileset
      - retained variant ID list
      - per-chromosome QC summary
      - PLINK2 log
    */

    output:

    tuple val(chr),
          path("chr${chr}.HQ_common.pgen"),
          path("chr${chr}.HQ_common.pvar"),
          path("chr${chr}.HQ_common.psam"),
          path("chr${chr}.HQ_common.snplist"),
          path("chr${chr}.HQ_common.summary.tsv"),
          emit: hq_common

    path "chr${chr}.HQ_common.log",
         emit: plink_logs


    script:

    """
    set -euo pipefail

    echo "============================================================"
    echo "GEL REGENIE shared genotype QC"
    echo "Chromosome: chr${chr}"
    echo "============================================================"
    echo
    echo "Parameters:"
    echo "  MAF threshold          : ${params.maf}"
    echo "  Missingness threshold  : ${params.geno}"
    echo "  HWE threshold          : ${params.hwe}"
    echo "  Threads                : ${task.cpus}"
    echo


    # -------------------------------------------------------------------------
    # Pre-QC sanity checks
    # -------------------------------------------------------------------------

    [[ -s dragen.pgen ]] || {
        echo "ERROR: missing or empty PGEN" >&2
        exit 1
    }

    [[ -s dragen.pvar ]] || {
        echo "ERROR: missing or empty PVAR" >&2
        exit 1
    }

    [[ -s dragen.psam ]] || {
        echo "ERROR: missing or empty PSAM" >&2
        exit 1
    }

    [[ -s basic_QC.variant_ids.txt ]] || {
        echo "ERROR: missing or empty basic-QC variant list" >&2
        exit 1
    }

    [[ -s GEL_CM_REGENIE.keep.txt ]] || {
        echo "ERROR: missing or empty sample keep file" >&2
        exit 1
    }


    input_variants=\$(wc -l < basic_QC.variant_ids.txt)
    input_samples=\$(wc -l < GEL_CM_REGENIE.keep.txt)

    duplicate_variant_ids=\$(
        sort basic_QC.variant_ids.txt |
        uniq -d |
        wc -l
    )

    duplicate_sample_iids=\$(
        awk '{print \$2}' GEL_CM_REGENIE.keep.txt |
        sort |
        uniq -d |
        wc -l
    )

    echo "Input basic-QC variants : \${input_variants}"
    echo "Input samples           : \${input_samples}"
    echo "Duplicate variant IDs   : \${duplicate_variant_ids}"
    echo "Duplicate sample IIDs   : \${duplicate_sample_iids}"
    echo


    if [[ "\${duplicate_variant_ids}" -ne 0 ]]; then
        echo "ERROR: duplicated IDs detected in variant list" >&2
        exit 1
    fi

    if [[ "\${duplicate_sample_iids}" -ne 0 ]]; then
        echo "ERROR: duplicated IIDs detected in sample keep file" >&2
        exit 1
    fi


    # -------------------------------------------------------------------------
    # Shared genotype QC
    # -------------------------------------------------------------------------

    plink2 \\
        --pfile dragen \\
        --keep GEL_CM_REGENIE.keep.txt \\
        --extract basic_QC.variant_ids.txt \\
        --snps-only just-acgt \\
        --maf ${params.maf} \\
        --geno ${params.geno} \\
        --hwe ${params.hwe} 0 \\
        --make-pgen \\
        --write-snplist \\
        --threads ${task.cpus} \\
        --out chr${chr}.HQ_common


    # -------------------------------------------------------------------------
    # Post-QC sanity checks
    # -------------------------------------------------------------------------

    for file in \\
        chr${chr}.HQ_common.pgen \\
        chr${chr}.HQ_common.pvar \\
        chr${chr}.HQ_common.psam \\
        chr${chr}.HQ_common.snplist \\
        chr${chr}.HQ_common.log
    do
        [[ -s "\${file}" ]] || {
            echo "ERROR: expected output missing or empty: \${file}" >&2
            exit 1
        }
    done


    output_variants=\$(wc -l < chr${chr}.HQ_common.snplist)

    output_samples=\$(
        awk '
            !/^#/ { n++ }
            END { print n+0 }
        ' chr${chr}.HQ_common.psam
    )


    if [[ "\${output_samples}" -ne "\${input_samples}" ]]; then
        echo "ERROR: output sample count differs from requested keep-list count" >&2
        echo "Input samples : \${input_samples}" >&2
        echo "Output samples: \${output_samples}" >&2
        exit 1
    fi


    if [[ "\${output_variants}" -eq 0 ]]; then
        echo "ERROR: no variants remained after QC" >&2
        exit 1
    fi


    proportion_retained=\$(
        awk \
            -v kept="\${output_variants}" \
            -v total="\${input_variants}" \
            'BEGIN {
                if (total == 0) {
                    print "NA"
                } else {
                    printf "%.6f", kept / total
                }
            }'
    )


    # -------------------------------------------------------------------------
    # Per-chromosome QC summary
    # -------------------------------------------------------------------------

    printf "chromosome\\tinput_basic_QC_variants\\tHQ_common_variants\\tproportion_retained\\tinput_samples\\toutput_samples\\n" \\
        > chr${chr}.HQ_common.summary.tsv

    printf "chr${chr}\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "\${input_variants}" \\
        "\${output_variants}" \\
        "\${proportion_retained}" \\
        "\${input_samples}" \\
        "\${output_samples}" \\
        >> chr${chr}.HQ_common.summary.tsv


    echo
    echo "============================================================"
    echo "chr${chr} shared QC completed"
    echo "============================================================"
    echo "Input basic-QC variants : \${input_variants}"
    echo "HQ common variants      : \${output_variants}"
    echo "Proportion retained     : \${proportion_retained}"
    echo "Samples                 : \${output_samples}"
    echo "============================================================"
    """
}


/*
-------------------------------------------------------------------------------
Workflow
-------------------------------------------------------------------------------
*/

workflow {

    /*
    Required runtime inputs
    */

    if (!params.pgen_root) {
        error """
        Missing required parameter: --pgen_root

        Expected structure:
          <pgen_root>/chrom-1/postproc-pgen/dragen.pgen
          <pgen_root>/chrom-1/postproc-pgen/dragen.pvar
          <pgen_root>/chrom-1/postproc-pgen/dragen.psam
          ...
          <pgen_root>/chrom-22/postproc-pgen/dragen.*
        """.stripIndent()
    }


    if (!params.variant_list_dir) {
        error """
        Missing required parameter: --variant_list_dir

        Expected structure:
          <variant_list_dir>/chr1.PASS_or_LowMLSQ.variant_ids.txt
          ...
          <variant_list_dir>/chr22.PASS_or_LowMLSQ.variant_ids.txt
        """.stripIndent()
    }


    if (!params.keep_file) {
        error """
        Missing required parameter: --keep_file

        Expected PLINK2 two-column FID/IID keep file.
        """.stripIndent()
    }


    /*
    Construct chromosome-specific input tuples.

    Each emitted item has:

      [
        chromosome,
        PGEN,
        PVAR,
        PSAM,
        basic-QC variant list
      ]
    */

    genotype_ch = Channel
        .fromList((1..22).toList())
        .map { chr ->

            def pgen = file(
                "${params.pgen_root}/chrom-${chr}/postproc-pgen/dragen.pgen",
                checkIfExists: true
            )

            def pvar = file(
                "${params.pgen_root}/chrom-${chr}/postproc-pgen/dragen.pvar",
                checkIfExists: true
            )

            def psam = file(
                "${params.pgen_root}/chrom-${chr}/postproc-pgen/dragen.psam",
                checkIfExists: true
            )

            def variant_list = file(
                "${params.variant_list_dir}/chr${chr}.PASS_or_LowMLSQ.variant_ids.txt",
                checkIfExists: true
            )

            tuple(
                chr,
                pgen,
                pvar,
                psam,
                variant_list
            )
        }


    /*
    The same sample keep file is reused by every chromosome task.
    */

    keep_ch = Channel.value(
        file(
            params.keep_file,
            checkIfExists: true
        )
    )


    /*
    Launch chr1-22 independently.
    */

    SHARED_QC(
        genotype_ch,
        keep_ch
    )
}