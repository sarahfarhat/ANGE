#!/usr/bin/env nextflow

process minimap2align {
    container 'community.wave.seqera.io/library/minimap2:2.30--dde6b0c5fbc82ebd'
    input:
    path genomemap
    path reads1
    path reads2

    output:
    path "transcripts.sam"

    script:
    """
    minimap2 -ax sr -t ${task.cpus} ${genomemap} ${reads1} ${reads2} -o transcripts.sam
    """
}
