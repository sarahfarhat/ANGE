#!/usr/bin/env nextflow

process fx2tab {
    container 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'
    input:
    path fasta

    output:
    path "${fasta.baseName}.tab"

    script:
    """
    seqkit fx2tab --only-id ${fasta} -o ${fasta.baseName}.tab
    """
}
