#!/usr/bin/env nextflow

// id carried through untouched (see cleanproteindb.nf) - runs once per source file, so
// redundancy is only removed *within* one source, never across sources (that's the whole
// point of the per-source design: no cross-file dedup, full traceability instead).
process clusterprot {
    container 'community.wave.seqera.io/library/cd-hit:4.8.1--8078bb96af799bf2'
    input:
    tuple val(id), path(proteins_file)
    val idper

    output:
    tuple val(id), path("cleaned_proteins.fasta")

    script:
    """
    cd-hit -i ${proteins_file} -o cleaned_proteins.fasta -c ${idper} -M 30000
    """
}
