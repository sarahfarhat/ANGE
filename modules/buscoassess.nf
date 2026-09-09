#!/usr/bin/env nextflow

process buscoassess {
    container 'community.wave.seqera.io/library/busco:6.1.0--6d1f7006d91892b3'
    // Closure (lazy) form: evaluated per task at run time, not once when this module is
    // parsed - see the note in nextflow.config on why that ordering matters here.
    publishDir { "${params.outdir}/busco" }, mode: 'copy'

    input:
    tuple val(label), path(fasta) // label names this run's output dir and its
                                   // aggregation-friendly summary copy, e.g. 'genome',
                                   // 'final_gmove', 'protein_branch_uniprotkb_D1' - now
                                   // carried alongside its fasta in one tuple so per-source
                                   // busco runs (one per protein-database/transcriptome file,
                                   // see main.nf) stay correctly paired with their own label.
    val mode         // busco -m: 'genome', 'transcriptome', or 'proteins'
    path lineage_dir // pre-downloaded lineage dataset directory - see buscodownload.nf. A
                      // local path (rather than a bare lineage name) tells busco to use it
                      // directly instead of trying to download its own copy.

    output:
    path "${label}_busco",             emit: dir
    path "${label}.short_summary.txt", emit: summary

    script:
    // The renamed copy of busco's own short_summary file (busco embeds the lineage name and
    // dataset name in its own filename, which varies run to run) gives buscosummary.nf a
    // predictable "<label>.short_summary.txt" name to key off of when building the combined
    // table across steps, instead of having to parse busco's own naming scheme.
    //
    // NB: "./${lineage_dir}" not "${lineage_dir}" - Nextflow stages the buscodownload
    // directory here using just its basename (e.g. "mollusca_odb12"), which is textually
    // indistinguishable from a bare lineage name. Busco's local-vs-download decision for -l
    // appears to key off whether the argument looks like a path (contains a "/"), so a bare
    // basename with no "/" was silently falling through to its own auto-download - the whole
    // point of buscodownload.nf, gone. The "./" prefix forces it to be recognized as the
    // staged local directory instead.
    """
    busco -i ${fasta} -o ${label}_busco -l ./${lineage_dir} -m ${mode} --cpu ${task.cpus}
    cp ${label}_busco/short_summary.specific.*.txt ${label}.short_summary.txt
    """
}
