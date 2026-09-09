#!/usr/bin/env nextflow

process mergetranscripts {
    input:
    path transcript_files

    output:
    path "transcripts_merged.fasta"

    script:
    """
    cat ${transcript_files} > transcripts_merged.fasta
    """
}
