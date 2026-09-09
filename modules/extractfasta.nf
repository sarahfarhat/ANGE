#!/usr/bin/env nextflow

process extractfasta {
    container 'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4'
    // Closure (lazy) form: evaluated per task at run time, not once when this module is
    // parsed - see the note in nextflow.config on why that ordering matters here.
    publishDir { "${params.outdir}/gmove" }, mode: 'copy'

    input:
    path genome
    path gff3

    output:
    path "gmove_final.cds.fasta",      emit: cds
    path "gmove_final.proteins.fasta", emit: proteins

    script:
    // NB: only run on gmove's own output, not the exonerate-derived intermediate
    // gene_models.*.gff3 files upstream. gmove is a real gene predictor, so its gff3 should
    // carry proper CDS features with phase - which gffread needs to translate correctly.
    // The intermediate files only have gene/exon features (no CDS/phase - see
    // exonerate2gff.nf), so gffread can't reliably extract a coding sequence from those; the
    // per-branch busco steps use extractgenefasta.nf's cruder raw-span extraction instead.
    // I couldn't run gmove here to confirm its gff3 actually has well-formed CDS/phase
    // fields - if gffread comes back with empty/near-empty fasta files, send me the head of
    // gmove_out.gff3 and I'll adjust.
    """
    gffread -g ${genome} -x gmove_final.cds.fasta -y gmove_final.proteins.fasta ${gff3}
    """
}
