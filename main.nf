#!/usr/bin/env nextflow


/*
============================================================
Extract VEP annotations from GEL AggV3 functional annotation
subshard VCFs

Input:
    functional_annotation_shards.csv

Output:
    results/VEP_annotation/*.tsv
    results/VEP_annotation/VEP_annotation_summary.tsv


Design:
    - one subshard = one task
    - one row = one VEP transcript consequence
    - VCFs are accessed through S3 paths
    - summary generated after all tasks complete


Author:
    Shiyu Zhang

Version:
    v1.0.0
============================================================
*/


nextflow.enable.dsl=2



params.manifest = null

params.outdir = "results"

params.extract_script =
    "${projectDir}/scripts/extract_v8.py"



if( !params.manifest ){

    error """

    Missing required parameter:

    --manifest functional_annotation_shards.csv

    """

}



workflow {


    /*
    ==========================================
    Read GEL shard manifest
    ==========================================
    */


    manifest_ch = Channel
        .fromPath(params.manifest)
        .splitCsv(
            header:true
        )
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


    /*
    ==========================================
    Collect all completed tasks
    ==========================================
    */


    extracted_ch
        .collect()
        .set { extracted_all }



    SUMMARY(
        extracted_all
    )

}




process EXTRACT_VEP_ANNOTATION {


    tag {

        "chr${chr}:shard${shard}:subshard${subshard}"

    }


    publishDir:

        "${params.outdir}/VEP_annotation",
        mode: "copy"



    input:


    tuple

    val(chr),
    val(start),
    val(end),
    val(region),
    val(shard),
    val(subshard),
    path(vcf),
    path(vcf_index)



    output:


    tuple

    val(chr),
    val(start),
    val(end),
    val(region),
    val(shard),
    val(subshard),
    val(vcf.simpleName),
    path("*.tsv")



    script:


    def output_name =
        "GEL.VEP.shard${shard}.subshard${subshard}.${region}.tsv"


    """

    python ${params.extract_script} \
        --vcf ${vcf} \
        --out ${output_name}

    """

}




process SUMMARY {


    tag:

    "Generate VEP annotation summary"



    publishDir:

        "${params.outdir}/VEP_annotation",
        mode: "copy"



    input:


    val(records)



    output:


    path(
        "VEP_annotation_summary.tsv"
    )



    script:


    def lines = records.collect {

        row ->


        def chr =
            row[0]

        def start =
            row[1]

        def end =
            row[2]

        def region =
            row[3]

        def shard =
            row[4]

        def subshard =
            row[5]

        def input_vcf =
            row[6]

        def output_tsv =
            row[7].getName()



        "${chr}\t${start}\t${end}\t${region}\t${shard}\t${subshard}\t${input_vcf}\t${output_tsv}\tPASS"


    }



    """

    cat > VEP_annotation_summary.tsv << EOF

chr	start	end	region	shard	subshard	input_vcf	output_tsv	status

EOF


    cat >> VEP_annotation_summary.tsv << EOF

${lines.join('\n')}

EOF


    """

}

