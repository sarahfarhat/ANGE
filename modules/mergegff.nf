#!/usr/bin/env nextflow

process mergegff {
    // Closure (lazy) form: evaluated per task at run time, not once when this module is
    // parsed - see the note in nextflow.config on why that ordering matters here.
    publishDir { "${params.outdir}/gene_models" }, mode: 'copy'

    input:
    path exonerate_gff3
    path to_remove_list
    path 'gff_parser.py'
    val branchLabel

    output:
    path "gene_models.${branchLabel}.gff3"

    script:
    // gff_parser.py names its own output after its input file's basename (e.g.
    // "<input>_NR.gff3"), which is identical for both branches upstream (both exonerate2gff
    // calls write "exonerate_out.gff3") - so it's copied here to a branch-labeled name
    // rather than relying on that derived name, which would collide once both branches'
    // outputs meet downstream (e.g. at gmove, or here in the shared publishDir).
    """
    python3 gff_parser.py ${exonerate_gff3} ${to_remove_list}
    cp \$(ls *_NR.gff3 | head -1) gene_models.${branchLabel}.gff3
    """
}
