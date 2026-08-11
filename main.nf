#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
===============================================================================

REGENIE Step 1: Whole-genome regression model

Phenotype:
    Cutaneous melanoma (CM)
    Binary trait:
        0 = control
        1 = case

Dataset:
    Genomics England (GEL)

Input genotype:
    High-quality LD-pruned WGS variants
    PLINK2 pgen format

Input files:
    <prefix>.pgen
    <prefix>.pvar
    <prefix>.psam

Phenotype:
    GEL_CM_REGENIE.phenotype.tsv

Covariates:
    GEL_CM_REGENIE.covariates.tsv

Covariates included:
    - genetic_sex
    - study_source
    - year_of_birth
    - PC1-PC20

Categorical covariates:
    genetic_sex
    study_source

Model:
    Binary trait (--bt)

Software:
    REGENIE v4.1.2

Container:
    ghcr.io/shiyuzhang0522/regenie:4.1.2

===============================================================================
*/


// -----------------------------------------------------------------------------
// Parameters
// -----------------------------------------------------------------------------

params.genotype_prefix = null

params.phenotype = null

params.covariates = null

params.outdir = "REGENIE.Step1.results"

params.out_prefix = "GEL_CM_REGENIE_step1"


params.bsize = 1000


// -----------------------------------------------------------------------------
// Process
// -----------------------------------------------------------------------------

process REGENIE_STEP1 {


    tag "${params.out_prefix}"


    container:
    "ghcr.io/shiyuzhang0522/regenie:4.1.2"


    cpus 16

    memory "100 GB"

    time "7d"



    publishDir "${params.outdir}",
        mode: "copy"



    input:

    tuple val(prefix),
          path(pgen),
          path(pvar),
          path(psam)

    path phenotype

    path covariates



    output:

    path "${params.out_prefix}*"



    script:


    """

    set -euo pipefail


    echo "============================================================"
    echo "REGENIE Step 1"
    echo "Start time: \$(date)"
    echo "Hostname:   \$(hostname)"
    echo "Threads:    ${task.cpus}"
    echo "============================================================"



    echo "Input genotype:"
    echo "${prefix}"


    echo "Samples:"
    awk 'NR>1{count++} END{print count}' ${psam}


    echo "Variants:"
    grep -v '^##' ${pvar} | grep -v '^#CHROM' | wc -l



    mkdir -p tmp



    regenie \\

        --step 1 \\

        --pgen ${prefix} \\

        --phenoFile ${phenotype} \\

        --phenoCol CM \\

        --covarFile ${covariates} \\

        --catCovarList genetic_sex,study_source \\

        --maxCatLevels 30 \\

        --bt \\

        --bsize ${params.bsize} \\

        --lowmem \\

        --lowmem-prefix tmp/regenie_step1_tmp_preds \\

        --threads ${task.cpus} \\

        --out ${params.out_prefix}



    echo "============================================================"
    echo "REGENIE Step 1 completed"
    echo "End time: \$(date)"
    echo "Output prefix:"
    echo "${params.out_prefix}"
    echo "============================================================"

    """

}



// -----------------------------------------------------------------------------
// Workflow
// -----------------------------------------------------------------------------

workflow {


    /*
    ---------------------------------------------------------------------------
    Validate parameters
    ---------------------------------------------------------------------------
    */


    if( !params.genotype_prefix ) {
        error "Missing parameter: --genotype_prefix"
    }


    if( !params.phenotype ) {
        error "Missing parameter: --phenotype"
    }


    if( !params.covariates ) {
        error "Missing parameter: --covariates"
    }



    /*
    ---------------------------------------------------------------------------
    Input channels
    ---------------------------------------------------------------------------
    */


    genotype_prefix = file(params.genotype_prefix)


    genotype_base =
        params.genotype_prefix



    pgen =
        file("${genotype_base}.pgen",
            checkIfExists:true)


    pvar =
        file("${genotype_base}.pvar",
            checkIfExists:true)


    psam =
        file("${genotype_base}.psam",
            checkIfExists:true)



    phenotype =
        file(params.phenotype,
            checkIfExists:true)


    covariates =
        file(params.covariates,
            checkIfExists:true)



    genotype_ch =
        Channel.of(
            tuple(
                genotype_base,
                pgen,
                pvar,
                psam
            )
        )



    /*
    ---------------------------------------------------------------------------
    Run REGENIE Step 1
    ---------------------------------------------------------------------------
    */


    REGENIE_STEP1(
        genotype_ch,
        phenotype,
        covariates
    )

}