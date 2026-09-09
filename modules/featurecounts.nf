#!/usr/bin/env nextflow

process featurecounts {
    container 'quay.io/biocontainers/subread:2.1.1--h577a1d6_0'
    input:
    path gtf
    path bam
    path bai

    output:
    path "exonerate_counts.txt"

    script:
    // NB: the doc's draft reused the .gtf filename for -o, which would overwrite the input;
    // this writes counts to their own file instead.
    """
    featureCounts -a ${gtf} -o exonerate_counts.txt -t transcript -p --donotsort -M -T ${task.cpus} ${bam}
    """
}
