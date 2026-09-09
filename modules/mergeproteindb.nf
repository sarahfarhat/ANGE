#!/usr/bin/env nextflow

// Mirrors mergetranscripts.nf: params.proteindatabase can be a glob or a comma-separated
// list of protein fasta files (see the Channel.fromPath(...) call in main.nf), so this
// concatenates them into one file before cleanproteindb runs. Unlike mergetranscripts,
// individual input files may be gzipped (cleanproteindb previously handled that itself for
// the single-file case) - each file is decompressed on the fly here so cleanproteindb
// always receives a single plain-text fasta regardless of how many inputs it was built from.
process mergeproteindb {
    input:
    path protein_files

    output:
    path "proteindb_merged.fasta"

    script:
    """
    for f in ${protein_files}; do
        if [[ "\$f" == *.gz ]]; then
            zcat "\$f" >> proteindb_merged.fasta
        else
            cat "\$f" >> proteindb_merged.fasta
        fi
    done
    """
}
