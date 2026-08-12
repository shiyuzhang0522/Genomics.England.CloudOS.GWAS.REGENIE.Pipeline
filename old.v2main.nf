nextflow.enable.dsl = 2


/*
===============================================================================

Genomics England CloudOS GWAS / REGENIE Pipeline

Task:
REGENIE Step 2 Genome-wide Association Testing


Phenotype:
-----------
Cutaneous melanoma (CM)

Binary trait:
0 = control
1 = case


Dataset:
--------
Genomics England (GEL) AggV3 WGS


Purpose:
--------
Perform chromosome-specific GWAS using REGENIE Step 2
with LOCO predictions generated from REGENIE Step 1.


Input genotype:
---------------
Chromosome-specific AggV3 PGEN files:

    dragen.pgen
    dragen.pvar
    dragen.psam


Variant QC:
-----------
Only variants passing site-level QC are tested:

    FILTER = PASS or LowMLSQ

Provided through:

    chrN.PASS_or_LowMLSQ.variant_ids.txt


Association model:
------------------
REGENIE Step 2:

    - binary logistic regression
    - additive model
    - LOCO prediction adjustment
    - approximate Firth correction


Covariates:
-----------
Included:

    - genetic_sex
    - study_source
    - year_of_birth
    - PC1-PC20


Categorical covariates:

    - genetic_sex
    - study_source


Parallelization:
----------------
chr1-chr22 independently


Software:
---------
REGENIE v4.1.2


Container:
----------
Defined externally in nextflow.config


Author:
-------
Shelley

===============================================================================
*/


// ============================================================================
// Parameters
// ============================================================================


params.pgen_root = null

params.variant_list_dir = null

params.pheno_file = null

params.covar_file = null


/*
REGENIE Step 1 LOCO prediction file

A local prediction list is generated inside each task
because the original Nextflow-generated pred.list contains
temporary /tmp/nxf.* paths.
*/

params.step1_loco = null


params.outdir = "regenie.step2.GWAS"


params.bsize = 400



// ============================================================================
// Process: REGENIE Step 2
// ============================================================================


process REGENIE_STEP2 {


    tag "chr${chr}"


    /*
    Retry only transient infrastructure failures.

    137:
        SIGKILL, commonly OOM

    143:
        SIGTERM, commonly AWS EC2 termination
    */

    errorStrategy {

        task.exitStatus in [137,143] ?
            'retry' :
            'terminate'

    }


    maxRetries = 3



    cpus 16


    memory "40 GB"


    time "7d"



    publishDir "${params.outdir}/chr${chr}",
        mode: "copy"



    input:


    tuple val(chr),
          path(pgen, name: "dragen.pgen"),
          path(pvar, name: "dragen.pvar"),
          path(psam, name: "dragen.psam"),
          path(variant_list, name: "variant_ids.txt")


    path phenotype,
         name: "GEL_CM_REGENIE.phenotype.tsv"


    path covariates,
         name: "GEL_CM_REGENIE.covariates.tsv"


    path loco,
         name: "GEL_CM_REGENIE_step1_1.loco"



    output:


    tuple val(chr),
          path("chr${chr}.GEL_CM_REGENIE_step2_CM.regenie.gz"),
          emit: regenie_results


    tuple val(chr),
          path("chr${chr}.GEL_CM_REGENIE_step2.log"),
          emit: regenie_logs



    script:


    """

    set -euo pipefail


    echo "============================================================"
    echo "REGENIE Step 2 GWAS"
    echo "============================================================"

    echo "Chromosome: ${chr}"
    echo "Attempt:    ${task.attempt}"
    echo "Task hash:  ${task.hash}"
    echo "Hostname:   \\$(hostname)"
    echo "Date:       \\$(date)"
    echo "Threads:    ${task.cpus}"

    echo
    echo "Memory before REGENIE:"
    free -h || true


    echo
    echo "Working directory:"
    pwd


    echo
    echo "Input files"
    echo "------------------------------------------------------------"

    ls -lh dragen.pgen
    ls -lh dragen.pvar
    ls -lh dragen.psam

    echo
    wc -l variant_ids.txt

    ls -lh GEL_CM_REGENIE_step1_1.loco



    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------


    for file in \\
        dragen.pgen \\
        dragen.pvar \\
        dragen.psam \\
        variant_ids.txt \\
        GEL_CM_REGENIE.phenotype.tsv \\
        GEL_CM_REGENIE.covariates.tsv \\
        GEL_CM_REGENIE_step1_1.loco

    do

        if [[ ! -s "\${file}" ]]; then

            echo "ERROR: missing or empty input file:"
            echo "\${file}"

            exit 1

        fi

    done



    # -------------------------------------------------------------------------
    # Generate local REGENIE prediction list
    # -------------------------------------------------------------------------


    cat > GEL_CM_REGENIE_step1_pred.list << EOF
CM GEL_CM_REGENIE_step1_1.loco
EOF



    echo
    echo "Prediction list:"
    cat GEL_CM_REGENIE_step1_pred.list



    # -------------------------------------------------------------------------
    # REGENIE Step 2
    # -------------------------------------------------------------------------


    regenie \\
        --step 2 \\
        --pgen dragen \\
        --extract variant_ids.txt \\
        --phenoFile GEL_CM_REGENIE.phenotype.tsv \\
        --phenoCol CM \\
        --covarFile GEL_CM_REGENIE.covariates.tsv \\
        --catCovarList genetic_sex,study_source \\
        --maxCatLevels 30 \\
        --pred GEL_CM_REGENIE_step1_pred.list \\
        --bt \\
        --minMAC 20 \\
        --firth \\
        --approx \\
        --firth-se \\
        --pThresh 0.01 \\
        --bsize ${params.bsize} \\
        --threads ${task.cpus} \\
        --gz \\
        --out chr${chr}.GEL_CM_REGENIE_step2



    echo
    echo "============================================================"
    echo "POST-REGENIE OUTPUT CHECK"
    echo "============================================================"


    GWAS_FILE="chr${chr}.GEL_CM_REGENIE_step2_CM.regenie.gz"

    LOG_FILE="chr${chr}.GEL_CM_REGENIE_step2.log"



    echo
    ls -lah



    if [[ ! -s "\${GWAS_FILE}" ]]; then

        echo "ERROR:"
        echo "Missing REGENIE GWAS output:"
        echo "\${GWAS_FILE}"

        exit 1

    fi



    if [[ ! -s "\${LOG_FILE}" ]]; then

        echo "ERROR:"
        echo "Missing REGENIE log file:"
        echo "\${LOG_FILE}"

        exit 1

    fi



    echo
    echo "REGENIE output successfully generated"

    ls -lh "\${GWAS_FILE}"
    ls -lh "\${LOG_FILE}"



    echo
    echo "Memory after REGENIE:"
    free -h || true


    echo
    echo "chr${chr} completed successfully"

    """

}



// ============================================================================
// Workflow
// ============================================================================


workflow {


    // -------------------------------------------------------------------------
    // Validate parameters
    // -------------------------------------------------------------------------


    if (!params.pgen_root)
        error "Missing --pgen_root"


    if (!params.variant_list_dir)
        error "Missing --variant_list_dir"


    if (!params.pheno_file)
        error "Missing --pheno_file"


    if (!params.covar_file)
        error "Missing --covar_file"


    if (!params.step1_loco)
        error "Missing --step1_loco"



    // -------------------------------------------------------------------------
    // Chromosome-specific genotype channel
    // -------------------------------------------------------------------------


    genotype_ch =
        Channel
            .fromList((1..22).toList())
            .map { chr ->


                def pgen =
                    file(
                        "${params.pgen_root}/chrom-${chr}/postproc-pgen/dragen.pgen",
                        checkIfExists: true
                    )


                def pvar =
                    file(
                        "${params.pgen_root}/chrom-${chr}/postproc-pgen/dragen.pvar",
                        checkIfExists: true
                    )


                def psam =
                    file(
                        "${params.pgen_root}/chrom-${chr}/postproc-pgen/dragen.psam",
                        checkIfExists: true
                    )


                def variant_list =
                    file(
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



    phenotype_ch =
        Channel.value(
            file(params.pheno_file, checkIfExists:true)
        )


    covariates_ch =
        Channel.value(
            file(params.covar_file, checkIfExists:true)
        )


    loco_ch =
        Channel.value(
            file(params.step1_loco, checkIfExists:true)
        )



    REGENIE_STEP2(
        genotype_ch,
        phenotype_ch,
        covariates_ch,
        loco_ch
    )


}