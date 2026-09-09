#!/usr/bin/env nextflow

// Runs once per input source file (protein database or transcript assembly - see the
// sourceId()-tagged channels in main.nf) rather than once on a pre-merged file, so the id
// travels alongside the sequences through every downstream step and each source's evidence
// stays traceable all the way to the final gff (see the evidence_source attribute added in
// mergegff.nf).
process cleanproteindb {
    input:
    tuple val(id), path(proteins_file)

    output:
    tuple val(id), path("proteins.fasta")

    script:
    """
    if [[ "${proteins_file}" == *.gz ]]; then
        zcat "${proteins_file}" | awk '{if(\$0 ~ /^>/) gsub(/[\\/|,\\-+ \\(\\)\\[\\]]/, "_"); print}' > proteins.fasta
    else
        awk '{if(\$0 ~ /^>/) gsub(/[\\/|,\\-+ \\(\\)\\[\\]]/, "_"); print}' "${proteins_file}" > proteins.fasta
    fi
    """
}
