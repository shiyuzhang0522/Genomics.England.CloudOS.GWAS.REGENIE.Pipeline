#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
============================================================
GEL AggV3 VEP Annotation Extractor

Input:
    functional_annotation_shards.csv

Output:
    results/VEP_annotation/*.tsv
    results/VEP_annotation/VEP_annotation_summary.tsv

Design:
    - one subshard = one task
    - one row = one transcript consequence
    - GEL S3 VCFs are staged automatically by Nextflow
    - no manual S3 download or mounting
    - one final summary table after all subshards complete

Author:
    Shiyu Zhang

Version:
    v1.0.3
============================================================
*/


/*
============================================================
Parameters
============================================================
*/

params.manifest = null
params.outdir = "results"

params.extract_script =
    "${projectDir}/WGS-RV-pipeline/RV-protein-coding/Extract.VEP.Annotations/scripts/extract_v8.py"


if( !params.manifest ) {
    error """
Missing required parameter:

--manifest functional_annotation_shards.csv
"""
}


/*
============================================================
Process: extract VEP annotations
============================================================
*/

process EXTRACT_VEP_ANNOTATION {

    tag {
        "${chr}:${region}"
    }

    cpus 4
    memory "8 GB"
    time "48h"

    publishDir "${params.outdir}/VEP_annotation", mode: "copy"


    input:

    tuple val(chr),
          val(start),
          val(end),
          val(region),
          val(shard),
          val(subshard),
          path(vcf),
          path(vcf_index)

    path extract_script


    output:

    tuple val(chr),
          val(start),
          val(end),
          val(region),
          val(shard),
          val(subshard),
          path("*.tsv")


    script:

    def outfile =
        "GEL.VEP.shard${shard}.subshard${subshard}.${region}.tsv"

    """
    python ${extract_script} \
        --vcf ${vcf} \
        --out ${outfile}
    """
}


/*
============================================================
Process: summary collector
============================================================
*/

process SUMMARY {

    tag {
        "VEP_annotation_summary"
    }

    publishDir "${params.outdir}/VEP_annotation", mode: "copy"


    input:

    val(rows)


    output:

    path("VEP_annotation_summary.tsv")


    script:

    """
    cat > VEP_annotation_summary.tsv << 'EOF'
chr\tstart\tend\tregion\tshard\tsubshard\toutput_tsv\tstatus
EOF

    cat >> VEP_annotation_summary.tsv << 'EOF'
${rows.join('\n')}
EOF
    """
}


/*
============================================================
Workflow
============================================================
*/

workflow {

    /*
    --------------------------------------------------------
    Stage extractor script from the GitHub repository
    --------------------------------------------------------
    */

    extract_script_ch = Channel.value(
        file(params.extract_script)
    )


    /*
    --------------------------------------------------------
    Read GEL functional annotation manifest
    --------------------------------------------------------
    */

    manifest_ch = Channel
        .fromPath(params.manifest)
        .splitCsv(header: true)
        .map { row ->

            tuple(
                row.chr,
                row.start,
                row.end,
                row.region,
                row.shard,
                row.subshard,
                file(row.func_anno_vcf),
                file(row.func_anno_vcf_index)
            )

        }


    /*
    --------------------------------------------------------
    Run one extraction task per subshard
    --------------------------------------------------------
    */

    extracted_ch = EXTRACT_VEP_ANNOTATION(
        manifest_ch,
        extract_script_ch
    )


    /*
    --------------------------------------------------------
    Convert each completed task into one summary-table row
    and collect all rows into a single value channel
    --------------------------------------------------------
    */

    summary_rows_ch = extracted_ch
        .map { chr, start, end, region, shard, subshard, tsv ->

            "${chr}\t${start}\t${end}\t${region}\t${shard}\t${subshard}\t${tsv.getName()}\tPASS"

        }
        .collect()


    /*
    --------------------------------------------------------
    Generate one master summary table
    --------------------------------------------------------
    */

    SUMMARY(
        summary_rows_ch
    )
}
