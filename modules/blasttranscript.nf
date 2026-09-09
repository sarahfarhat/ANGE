#!/usr/bin/env nextflow

// id carried through untouched - see blastprot.nf for the same pattern on the protein branch.
process blasttranscript {
    container 'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd'
    input:
    path genome_msk
    tuple val(id), path(transcripts)

    output:
    // NB: named after the input chunk (not a fixed name) - blasttranscript runs once per
    // chunk from main.nf's splitFasta, and filterblast collects every chunk's output into
    // one task; a fixed output name collides across chunks there ("input file name collision").
    tuple val(id), path("${transcripts.baseName}.blastn.out")

    script:
    """
    makeblastdb -in ${genome_msk} -dbtype nucl -out genome_db
    blastn -db genome_db -query ${transcripts} -out ${transcripts.baseName}.blastn.out -evalue 0.001 -outfmt 6 -num_threads ${task.cpus}
    """
}
