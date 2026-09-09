#!/usr/bin/env nextflow

process samtoolssort {
    container 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'
    input:
    path sam

    output:
    path "transcripts_sorted.bam", emit: bam
    path "transcripts_sorted.bam.bai", emit: bai

    script:
    """
    samtools view -Sb ${sam} -o transcripts.bam
    samtools sort transcripts.bam -o transcripts_sorted.bam -@ ${task.cpus}
    samtools index transcripts_sorted.bam
    rm -f transcripts.bam
    """
}
