#!/usr/bin/env nextflow

process minimap2index {
    container 'community.wave.seqera.io/library/minimap2:2.30--dde6b0c5fbc82ebd'
    input:
    path genome

    output:
    path "genomemap.mmi"

    script:
    """
    minimap2 -d genomemap.mmi ${genome}
    """
}
