#!/usr/bin/env nextflow

process clustertranscript {
    container 'community.wave.seqera.io/library/cd-hit:4.8.1--8078bb96af799bf2'
    input:
    path transcripts_file
    val idper

    output:
    path "cleaned_transcripts.fasta"

    script:
    """
    cd-hit-est -i ${transcripts_file} -o cleaned_transcripts.fasta -c ${idper} -T ${task.cpus} -M 30000
    """
}
