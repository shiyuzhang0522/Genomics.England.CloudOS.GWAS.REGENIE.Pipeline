nextflow.enable.dsl = 2

process PLINK2_VERSION {
    container params.plink2_image

    output:
    path 'plink2.version.txt'

    script:
    """
    plink2 --version | tee plink2.version.txt
    """
}

workflow {
    PLINK2_VERSION()
}
