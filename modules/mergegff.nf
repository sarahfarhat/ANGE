#!/usr/bin/env nextflow

// The `val branchLabel` ('protein'/'transcript') is gone - output is now named after the
// source id itself, so results/gene_models/ gets one file per input database/assembly file
// (e.g. gene_models.uniprotkb_D1.gff3, gene_models.spades_mimit05_transcripts_renamed.gff3)
// instead of one per branch. On top of that, every feature gets an evidence_source=<id>
// GFF3 attribute appended, so origin stays recoverable even after main.nf's combineevidence
// step concatenates every source's file into the single --prot / --rna evidence gmove needs.
process mergegff {
    // Closure (lazy) form: evaluated per task at run time, not once when this module is
    // parsed - see the note in nextflow.config on why that ordering matters here.
    publishDir { "${params.outdir}/gene_models" }, mode: 'copy'

    input:
    tuple val(id), path(exonerate_gff3)
    path to_remove_list
    path 'gff_parser.py'

    output:
    tuple val(id), path("gene_models.${id}.gff3")

    script:
    """
    python3 gff_parser.py ${exonerate_gff3} ${to_remove_list}
    awk 'BEGIN{OFS="\t"}{\$9=\$9";evidence_source=${id}"; print}' \$(ls *_NR.gff3 | head -1) > gene_models.${id}.gff3
    """
}
