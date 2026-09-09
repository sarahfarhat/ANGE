#!/usr/bin/env nextflow

process gffreadconvert {
    container 'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4'
    input:
    path gff3

    output:
    path "exonerate_out.gtf"

    script:
    """
    gffread ${gff3} -T -o exonerate_out.gtf
    """
}
