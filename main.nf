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
    - VCF accessed directly from GEL S3
    - no intermediate download management


Author:
    Shiyu Zhang

Version:
    v1.0.2

============================================================
*/


params.manifest = null

params.outdir = "results"


params.extract_script =
"${projectDir}/WGS-RV-pipeline/RV-protein-coding/Extract.VEP.Annotations/scripts/extract_v8.py"



if( !params.manifest ){

    error """
Missing required parameter:

--manifest functional_annotation_shards.csv
"""
}



workflow {


    /*
    ========================================================
    Read GEL functional annotation manifest
    ========================================================
    */


    manifest_ch = Channel
        .fromPath(params.manifest)
        .splitCsv(header:true)
        .map { row ->

            tuple(
                row.chr,
                row.start,
                row.end,
                row.region,
                row.shard,
                row.subshard,
                row.func_anno_vcf,
                row.func_anno_vcf_index
            )
        }


    extracted_ch = EXTRACT_VEP_ANNOTATION(
        manifest_ch
    )


    extracted_ch
        .collect()
        .set { all_results }


    SUMMARY(
        all_results
    )

}




/*
============================================================
Extract VEP annotations
============================================================
*/


process EXTRACT_VEP_ANNOTATION {


    tag {

        "${chr}:${region}"

    }


    cpus 4

    memory "8 GB"

    time "48h"



    publishDir:

        "${params.outdir}/VEP_annotation",
        mode: "copy"



    input:


    tuple val(chr),
          val(start),
          val(end),
          val(region),
          val(shard),
          val(subshard),
          path(vcf),
          path(vcf_index)



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

    python ${params.extract_script} \
        --vcf ${vcf} \
        --out ${outfile}

    """

}



/*
============================================================
Summary collector
============================================================
*/


process SUMMARY {


    tag {

        "VEP_annotation_summary"

    }



    publishDir:

        "${params.outdir}/VEP_annotation",
        mode: "copy"



    input:


    val(records)



    output:


    path("VEP_annotation_summary.tsv")



    script:


    def rows = records.collect { r ->

        "${r[0]}\t${r[1]}\t${r[2]}\t${r[3]}\t${r[4]}\t${r[5]}\t${r[6].getName()}\tPASS"

    }



    """

cat > VEP_annotation_summary.tsv << EOF
chr\tstart\tend\tregion\tshard\tsubshard\toutput_tsv\tstatus
EOF


cat >> VEP_annotation_summary.tsv << EOF
${rows.join('\n')}
EOF

"""

}