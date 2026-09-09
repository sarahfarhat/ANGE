#!/usr/bin/env nextflow

process clusterprot {
    container 'community.wave.seqera.io/library/cd-hit:4.8.1--8078bb96af799bf2'
    input:
    path proteins_file
    val idper

    output:
    path "cleaned_proteins.fasta" 

    script:
    """
    cd-hit -i ${proteins_file} -o cleaned_proteins.fasta -c ${idper} -M 30000
    """
}