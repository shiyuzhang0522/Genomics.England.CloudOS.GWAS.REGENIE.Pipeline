#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
===============================================================================
GEL CloudOS | REGENIE Step 2 rare melanocyte cCRE association analysis
Author: Shelley

Purpose
-------
Reuse the working GEL GWAS PGEN/phenotype/covariate/LOCO framework while
preserving the supplied UKBB rare-cCRE mask and association settings.

Default design
--------------
22 autosomes x 5 masks x 3 AAF thresholds = 330 independent association tasks.
Each task analyzes a complete chromosome, one mask, and one AAF threshold.
The same mask files and fixed gnomAD exclusion list serve all three AAF runs.

Scientific settings
-------------------
--aaf-bins AAF --vc-maxAAF AAF --vc-tests skato,acato-full
--bt --firth --firth-se --approx --pThresh 0.05 --bsize 400
--write-mask-snplist --check-burden-files
As in the tested GEL workflow, --vc-MACthr 10 is explicit.
--minMAC and --build-mask are left at the REGENIE defaults, as in UKBB.
Other REGENIE defaults are retained. In particular, the
standard singleton-mask behavior is not disabled.

Inputs
------
All seven input locations are supplied via CloudOS parameters (see config).
IDs in annotations, set lists, QC lists and exclusions must match PVAR IDs.
Precomputed cCRE masks are used directly; this workflow does not rebuild
annotations or alter the cCRE definitions.

Outputs
-------
UKBB-style <outdir>/<mask>/AAF_<AAF>/chrN.Melanocyte_cCRE.<mask>.AAF_<AAF>* filenames.
Association results and mask SNP lists are uncompressed. REGENIE logs,
burden-file reports, console logs, resolved commands, local prediction lists,
and per-task QC summaries are retained under the corresponding logs/ folder.
Nextflow trace/report/timeline are retained under <outdir>/pipeline_info/.

Validation scope
----------------
Validate file availability, small-file schemas, input row counts, and output
structure. REGENIE performs sample matching and burden-input checking.
The QC TSV records result rows and represented cCRE regions; it does not certify
convergence or a nonmissing P value for every test. Inspect REGENIE logs and
burden reports. Pipeline logs can contain controlled-data identifiers.

Documentation
-------------
https://rgcgithub.github.io/regenie/options/
https://www.nextflow.io/docs/latest/process.html

Examples (supply all seven required input locations through CloudOS)
------------------------------------------------------------------
Full run: defaults in nextflow.config.
Pilot: --chromosomes 21 --masks CADD --aaf_thresholds 0.01
Use Nextflow -resume when restarting with the same work directory.
===============================================================================
*/


// ============================================================================
// Step 1. Validate chromosome and mask/AAF selections
// ============================================================================

def csvValues(value) {
    def items = value instanceof Collection ? value : value.toString().split(',').toList()
    items.collect { it.toString().trim() }
}


def chromosomeValues(value) {
    def result = []
    csvValues(value).each { token ->
        if (!(token ==~ /[0-9]+(-[0-9]+)?/)) {
            error "Invalid chromosome selection '${token}'; use 1-22 or 21,22."
        }
        def bounds = token.split('-').collect { it.toInteger() }
        def lo = bounds[0]
        def hi = bounds.size() == 2 ? bounds[1] : lo
        if (lo < 1 || hi > 22 || lo > hi) {
            error "Chromosome selection must be within 1-22: '${token}'."
        }
        result.addAll((lo..hi).toList())
    }
    if (result.size() != result.unique(false).size()) {
        error 'Chromosome selection contains duplicates or overlapping ranges.'
    }
    result.sort()
}


// ============================================================================
// Step 2. One REGENIE task per chromosome x mask x AAF
// ============================================================================

process REGENIE_RARE_CCRE {

    tag "chr${chr}:${mask}:AAF_${aaf}"

    publishDir "${params.outdir}/${mask}/AAF_${aaf}",
        mode: 'copy',
        overwrite: true,
        failOnError: true,
        saveAs: { name ->
            (name.endsWith('.regenie') || name.endsWith('_masks.snplist')) ?
                name : "logs/${name}"
        }

    input:

    tuple val(chr),
          path(pgen, name: 'dragen.pgen'),
          path(pvar, name: 'dragen.pvar'),
          path(psam, name: 'dragen.psam'),
          path(qc_list, name: 'variant_ids.txt'),
          path(exclude_list, name: 'gnomad_exclude.txt'),
          val(mask),
          val(aaf),
          path(annotation, name: 'ccre.annotation.txt'),
          path(setlist, name: 'ccre.setlist.txt'),
          path(maskdef, name: 'ccre.maskdef.txt')

    path phenotype, name: 'GEL_CM_REGENIE.phenotype.tsv'
    path covariates, name: 'GEL_CM_REGENIE.covariates.tsv'
    path loco, name: 'GEL_CM_REGENIE_step1_1.loco'

    output:

    tuple val(chr), val(mask), val(aaf),
          path("chr${chr}.Melanocyte_cCRE.${mask}.AAF_${aaf}_CM.regenie"), emit: results

    tuple val(chr), val(mask), val(aaf),
          path("chr${chr}.Melanocyte_cCRE.${mask}.AAF_${aaf}_masks.snplist"), emit: mask_snplists

    tuple val(chr), val(mask), val(aaf),
          path("chr${chr}.Melanocyte_cCRE.${mask}.AAF_${aaf}.qc.tsv"), emit: qc_summaries

    tuple val(chr), val(mask), val(aaf),
          path("chr${chr}.Melanocyte_cCRE.${mask}.AAF_${aaf}.log"),
          path("chr${chr}.Melanocyte_cCRE.${mask}.AAF_${aaf}.console.log"),
          path("chr${chr}.Melanocyte_cCRE.${mask}.AAF_${aaf}.command.sh"),
          path("chr${chr}.Melanocyte_cCRE.${mask}.AAF_${aaf}.pred.list"), emit: logs

    // The report may be absent when REGENIE finds nothing to report.
    tuple val(chr), val(mask), val(aaf),
          path("chr${chr}.Melanocyte_cCRE.${mask}.AAF_${aaf}_masks_report.txt"),
          optional: true, emit: burden_reports

    script:

    def prefix = "chr${chr}.Melanocyte_cCRE.${mask}.AAF_${aaf}"

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Run the full task body under a checked Bash script. This captures all
    # diagnostics in a declared output and propagates failures through tee.
    cat > run_task.sh <<'TASK'
    #!/usr/bin/env bash
    set -euo pipefail

    echo "============================================================"
    echo "GEL REGENIE rare-cCRE association: ${prefix}"
    echo "Attempt: ${task.attempt} | CPUs: ${task.cpus}"
    echo "Container: ${params.regenie_image}"
    date -u
    hostname
    command -v regenie
    regenie --version
    free -h || true

    # ------------------------------------------------------------------------
    # 1. Input availability and lightweight format checks
    # ------------------------------------------------------------------------

    for input_file in dragen.pgen dragen.pvar dragen.psam variant_ids.txt \\
        GEL_CM_REGENIE.phenotype.tsv GEL_CM_REGENIE.covariates.tsv \\
        GEL_CM_REGENIE_step1_1.loco ccre.annotation.txt \\
        ccre.setlist.txt ccre.maskdef.txt
    do
        if [[ ! -s "\${input_file}" ]]; then
            echo "[ERROR] Missing or empty input: \${input_file}" >&2
            exit 1
        fi
    done

    # A valid exclusion list can be empty. It must still be supplied/staged.
    [[ -f gnomad_exclude.txt ]] || { echo '[ERROR] Missing exclusion list' >&2; exit 1; }

    awk 'NR==1 {sub(/\\r\$/, ""); if (NF!=3 || \$1!="FID" || \$2!="IID" || \$3!="CM") exit 1; next}
         NF!=3 || (\$3!="0" && \$3!="1" && \$3!="NA") {exit 1}
         END {if(NR<2) exit 1}' GEL_CM_REGENIE.phenotype.tsv || {
        echo '[ERROR] Expected phenotype header FID IID CM and 0/1/NA coding.' >&2; exit 1;
    }

    awk 'NR==1 {sub(/\\r\$/, ""); if (NF!=25 || \$1!="FID" || \$2!="IID") exit 1;
         for(i=3;i<=NF;i++) h[\$i]++;
         if(h["genetic_sex"]!=1 || h["study_source"]!=1 || h["year_of_birth"]!=1) exit 1;
         for(i=1;i<=20;i++) if(h["PC" i]!=1) exit 1; exit 0}' \\
        GEL_CM_REGENIE.covariates.tsv || {
        echo '[ERROR] Expected the agreed 23 covariates plus FID and IID.' >&2; exit 1;
    }

    awk -v mask='${mask}' 'NF!=2 || \$1!=mask || \$2!=mask {bad=1}
         END {exit (bad || NR!=1)}' ccre.maskdef.txt || {
        echo '[ERROR] Mask definition does not match the selected mask.' >&2; exit 1;
    }

    awk -v mask='${mask}' 'NF!=3 || \$3!=mask {exit 1}' ccre.annotation.txt || {
        echo '[ERROR] Malformed annotation rows or incorrect mask label.' >&2; exit 1;
    }

    awk -v chr='${chr}' 'NF!=4 || \$2!=chr || \$3!~/^[0-9]+\$/ || \$3<1 {exit 1}' \\
        ccre.setlist.txt || {
        echo '[ERROR] Malformed set-list row or chromosome mismatch.' >&2; exit 1;
    }

    echo '[INFO] Input row counts (before REGENIE filtering):'
    wc -l variant_ids.txt gnomad_exclude.txt ccre.annotation.txt ccre.setlist.txt
    cat ccre.maskdef.txt

    # ------------------------------------------------------------------------
    # 2. Generate a prediction list using the task-local LOCO filename
    # ------------------------------------------------------------------------

    cat > GEL_CM_REGENIE_step1_pred.list <<'PRED'
    CM GEL_CM_REGENIE_step1_1.loco
    PRED

    cp GEL_CM_REGENIE_step1_pred.list '${prefix}.pred.list'
    cat GEL_CM_REGENIE_step1_pred.list

    # ------------------------------------------------------------------------
    # 3. Preserve the UKBB rare-cCRE analysis settings
    # ------------------------------------------------------------------------

    cat > '${prefix}.command.sh' <<'COMMAND'
    #!/usr/bin/env bash
    set -euo pipefail

    regenie \\
        --step 2 \\
        --pgen dragen \\
        --extract variant_ids.txt \\
        --exclude gnomad_exclude.txt \\
        --phenoFile GEL_CM_REGENIE.phenotype.tsv \\
        --phenoCol CM \\
        --covarFile GEL_CM_REGENIE.covariates.tsv \\
        --catCovarList genetic_sex,study_source \\
        --maxCatLevels 30 \\
        --pred GEL_CM_REGENIE_step1_pred.list \\
        --anno-file ccre.annotation.txt \\
        --set-list ccre.setlist.txt \\
        --mask-def ccre.maskdef.txt \\
        --aaf-bins ${aaf} \\
        --vc-maxAAF ${aaf} \\
        --vc-tests skato,acato-full \\
        --write-mask-snplist \\
        --check-burden-files \\
        --vc-MACthr 10 \\
        --bt \\
        --firth \\
        --firth-se \\
        --approx \\
        --pThresh 0.05 \\
        --bsize 400 \\
        --threads ${task.cpus} \\
        --out '${prefix}'
    COMMAND

    bash '${prefix}.command.sh'

    # ------------------------------------------------------------------------
    # 4. Require REGENIE outputs and summarize represented result cCRE regions
    # ------------------------------------------------------------------------

    for output_file in '${prefix}_CM.regenie' '${prefix}_masks.snplist' '${prefix}.log'
    do
        [[ -s "\${output_file}" ]] || {
            echo "[ERROR] Missing or empty output: \${output_file}" >&2; exit 1;
        }
    done

    # REGENIE may prepend ## mask metadata. Locate the actual header, then
    # count records, distinct cCRE IDs, and missing LOG10P values by TEST.
    awk -v chr='${chr}' -v mask='${mask}' -v aaf='${aaf}' '
        BEGIN {OFS="\\t"; print "chromosome","mask","aaf","test","result_rows","regions","missing_log10p"}
        /^##/ {next}
        !header {
            for(i=1;i<=NF;i++) {name=\$i; sub(/^#/,"",name); col[name]=i}
            if(!col["ID"] || !col["TEST"] || !col["LOG10P"]) {
                print "[ERROR] Missing ID/TEST/LOG10P result columns" > "/dev/stderr"; exit 1
            }
            width=NF; header=1; next
        }
        NF {
            if(NF!=width) {print "[ERROR] Inconsistent result column count" > "/dev/stderr"; bad=1; exit 1}
            test=\$(col["TEST"]); id=\$(col["ID"]); split(id,parts,"[.]"); region=parts[1]
            rows[test]++; total++
            if(!seen[test SUBSEP region]++) regions[test]++
            p=\$(col["LOG10P"])
            if(p=="NA" || p=="NaN" || p=="nan" || p=="-nan" || p==".") missing[test]++
        }
        END {
            if(bad || !header || !total) exit 1
            for(test in rows) print chr,mask,aaf,test,rows[test],regions[test],missing[test]+0
        }' '${prefix}_CM.regenie' > '${prefix}.qc.tsv'

    cat '${prefix}.qc.tsv'
    echo '[PASS] REGENIE finished; required outputs and result structure checked.'
    echo '[INFO] Review missing results, convergence messages, and the burden-file report.'
    date -u
    TASK

    bash run_task.sh 2>&1 | tee '${prefix}.console.log'
    """
}


// ============================================================================
// Step 3. Resolve inputs and construct the task grid
// ============================================================================

workflow {

    if (params.help) {
        log.info '''
GEL REGENIE rare-cCRE pipeline
Required: --pgen_root --variant_list_dir --gnomad_exclude_dir --pheno_file
          --covar_file --step1_loco --mask_dir
Default:  --chromosomes 1-22
          --masks CADD,GERP,JARVIS,FUNC_ALL,ALL
          --aaf_thresholds 0.01,0.001,0.0001
Output:   --outdir Rare-cCRE/Results
Pilot:    --chromosomes 21 --masks CADD --aaf_thresholds 0.01
Resources: --task_cpus 16 --task_memory '40 GB' --task_time 7d
CloudOS supplies the execution backend and work directory.
'''
        return
    }

    ['pgen_root', 'variant_list_dir', 'gnomad_exclude_dir', 'pheno_file',
     'covar_file', 'step1_loco', 'mask_dir'].each { name ->
        if (!params[name]?.toString()?.trim()) error "Missing required parameter --${name}"
    }

    if (!params.outdir?.toString()?.trim()) error '--outdir must not be empty.'
    if (!params.chromosomes || !params.masks || !params.aaf_thresholds) {
        error 'Chromosome, mask, and AAF selections must not be empty.'
    }
    def chromosomes = chromosomeValues(params.chromosomes)
    def masks = csvValues(params.masks)
    // CloudOS/Nextflow may parse a single AAF as a number (e.g. 1.0E-4).
    // Canonical decimal strings preserve the UKBB directory naming.
    def aafs = csvValues(params.aaf_thresholds).collect { token ->
        try {
            new BigDecimal(token).stripTrailingZeros().toPlainString()
        }
        catch (NumberFormatException ignored) {
            error "Invalid AAF threshold '${token}'."
        }
    }

    if (masks.any { !(it in ['CADD', 'GERP', 'JARVIS', 'FUNC_ALL', 'ALL']) } ||
        masks.unique(false).size() != masks.size()) {
        error 'Masks must be unique selections from CADD,GERP,JARVIS,FUNC_ALL,ALL.'
    }
    if (aafs.any { !(it in ['0.01', '0.001', '0.0001']) } ||
        aafs.unique(false).size() != aafs.size()) {
        error 'AAF thresholds must be unique selections from 0.01,0.001,0.0001.'
    }
    ['task_cpus'].each { name ->
        if (!(params[name].toString() ==~ /[1-9][0-9]*/)) error "--${name} must be a positive integer."
    }

    log.info "Planned REGENIE tasks: ${chromosomes.size()} x ${masks.size()} x ${aafs.size()} = ${chromosomes.size() * masks.size() * aafs.size()}"
    log.info "Chromosomes: ${chromosomes.join(',')} | Masks: ${masks.join(',')} | AAF: ${aafs.join(',')}"
    log.info "Container: ${params.regenie_image} | Results: ${params.outdir}"

    // Value channels broadcast the same phenotype/covariates/LOCO to all tasks.
    phenotype_ch = Channel.value(file(params.pheno_file, checkIfExists: true))
    covariates_ch = Channel.value(file(params.covar_file, checkIfExists: true))
    loco_ch = Channel.value(file(params.step1_loco, checkIfExists: true))

    // Resolve each chromosome's source inputs once, then reuse across masks/AAF.
    genotype_ch = Channel.fromList(chromosomes).map { chr ->
        def root = "${params.pgen_root}/chrom-${chr}/postproc-pgen"
        tuple(chr,
            file("${root}/dragen.pgen", checkIfExists: true),
            file("${root}/dragen.pvar", checkIfExists: true),
            file("${root}/dragen.psam", checkIfExists: true),
            file("${params.variant_list_dir}/chr${chr}.PASS_or_LowMLSQ.variant_ids.txt", checkIfExists: true),
            file("${params.gnomad_exclude_dir}/chr${chr}.variant.list", checkIfExists: true))
    }

    task_ch = genotype_ch.flatMap { chr, pgen, pvar, psam, qc, excluded ->
        masks.collectMany { mask ->
            def base = "${params.mask_dir}/chr${chr}/${mask}/chr${chr}.${mask}"
            def annotation = file("${base}.annotation.txt", checkIfExists: true)
            def setlist = file("${base}.setlist.txt", checkIfExists: true)
            def maskdef = file("${params.mask_dir}/chr${chr}/${mask}/${mask}.mask.def", checkIfExists: true)
            aafs.collect { aaf ->
                tuple(chr, pgen, pvar, psam, qc, excluded, mask, aaf, annotation, setlist, maskdef)
            }
        }
    }

    REGENIE_RARE_CCRE(task_ch, phenotype_ch, covariates_ch, loco_ch)
}
