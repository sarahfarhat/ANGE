#!/usr/bin/env nextflow

process cleanproteindb {
    input:
    path proteins_file

    output:
    path "proteins.fasta"

    script:
    """
    if [[ "${proteins_file}" == *.gz ]]; then
        zcat "${proteins_file}" | awk '{if(\$0 ~ /^>/) gsub(/[\\/|,\\-+ \\(\\)\\[\\]]/, "_"); print}' > proteins.fasta
    else
        awk '{if(\$0 ~ /^>/) gsub(/[\\/|,\\-+ \\(\\)\\[\\]]/, "_"); print}' "${proteins_file}" > proteins.fasta
    fi
    """
}