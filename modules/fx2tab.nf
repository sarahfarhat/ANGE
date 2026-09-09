#!/usr/bin/env nextflow

// Tuple-based (id, fasta) -> (id, tab) so per-source dataset lookups (fx2tab_dataset /
// fx2tab_transcripts, one call per protein-database / transcriptome-assembly file) can be
// joined back to their matching filterblast output by id in main.nf. The genome-wide call
// (fx2tab_genome) doesn't have a "source" of its own - main.nf passes a constant 'genome' id
// and immediately discards it after this step.
process fx2tab {
    container 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'
    input:
    tuple val(id), path(fasta)

    output:
    tuple val(id), path("${fasta.baseName}.tab")

    script:
    """
    seqkit fx2tab --only-id ${fasta} -o ${fasta.baseName}.tab
    """
}
