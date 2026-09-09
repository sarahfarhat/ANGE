#!/usr/bin/env nextflow

process buscosummary {
    publishDir { "${params.outdir}/busco" }, mode: 'copy'

    input:
    path summaries
    path 'busco_summary.py'

    output:
    path "busco_summary.tsv"

    script:
    """
    python3 busco_summary.py ${summaries}
    """
}
