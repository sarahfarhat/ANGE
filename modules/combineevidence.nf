#!/usr/bin/env nextflow

// Concatenates several per-source gff3 files (e.g. one gene_models.<id>.gff3 per protein
// database or transcript assembly file - see mergegff.nf) into one evidence file. Each
// feature still carries the evidence_source=<id> attribute mergegff.nf added, so origin
// remains recoverable even after this merge - used to build the single --prot / --rna
// evidence file gmove needs (and, only when params.run_low_count_filter is on, the single
// combined gff that step's RNA-seq read-count filter needs before mergegff runs).
process combineevidence {
    // Closure (lazy) form: evaluated per task at run time, not once when this module is
    // parsed - see the note in nextflow.config on why that ordering matters here.
    publishDir { "${params.outdir}/evidence_combined" }, mode: 'copy'

    input:
    path gff_files
    val label

    output:
    path "${label}.gff3"

    script:
    """
    cat ${gff_files} > ${label}.gff3
    """
}
