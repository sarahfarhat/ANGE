#!/usr/bin/env nextflow

// id carried through untouched - runs once per (source, chunk) pair, since main.nf's
// splitFasta(elem:1) preserves the source id on every chunk it produces from a per-source
// fasta. genome_msk stays a plain (non-tuple) input: it's the same file for every source and
// Nextflow automatically broadcasts a single-value input to every parallel task.
process blastprot {
    container 'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd'
    input:
    path genome_msk
    tuple val(id), path(proteins)

    output:
    // NB: named after the input chunk (not a fixed name) - blastprot runs once per chunk
    // from main.nf's splitFasta, and filterblast collects every chunk's output into one
    // task; a fixed output name collides across chunks there ("input file name collision").
    tuple val(id), path("${proteins.baseName}.tblastn.out")

    script:
    """
    makeblastdb -in ${genome_msk} -dbtype nucl -out genome_db
    tblastn -db genome_db -query ${proteins} -out ${proteins.baseName}.tblastn.out -evalue 0.001 -outfmt 6 -num_threads ${task.cpus}
    """
}
