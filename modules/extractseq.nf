#!/usr/bin/env nextflow

process extractseq {
    input:
    path genome_tab
    path dataset_tab
    path gff
    path 'extractseq_v4.py'

    output:
    path "matches_query", emit: query_dir
    path "matches_genome", emit: genome_dir

    script:
    """
    cp ${gff} matches.gff
    python3 extractseq_v4.py ${genome_tab} ${dataset_tab} matches.gff
    """
}
