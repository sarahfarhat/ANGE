#!/usr/bin/env nextflow

// genome_tab stays a plain (non-tuple) input: it's the same genome-wide table for every
// source, broadcast to every parallel task by Nextflow. dataset_tab and gff are per-source
// (id, path) results already joined by id in main.nf (fx2tab_dataset.out.join(filterblast.out)),
// so they arrive here as one 3-element tuple sharing that id.
process extractseq {
    input:
    path genome_tab
    tuple val(id), path(dataset_tab), path(gff)
    path 'extractseq_v4.py'

    output:
    tuple val(id), path("matches_query"), path("matches_genome")

    script:
    """
    cp ${gff} matches.gff
    python3 extractseq_v4.py ${genome_tab} ${dataset_tab} matches.gff
    """
}
