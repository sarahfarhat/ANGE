#!/usr/bin/env nextflow

// id carried through untouched (see cleanproteindb.nf) - runs once per source file, so
// redundancy is only removed *within* one source, never across sources.
process clustertranscript {
    container 'community.wave.seqera.io/library/cd-hit:4.8.1--8078bb96af799bf2'
    input:
    tuple val(id), path(transcripts_file)
    val idper

    output:
    tuple val(id), path("cleaned_transcripts.fasta")

    script:
    """
    cd-hit-est -i ${transcripts_file} -o cleaned_transcripts.fasta -c ${idper} -T ${task.cpus} -M 30000
    """
}
