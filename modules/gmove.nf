#!/usr/bin/env nextflow

process gmove {
    conda 'bioconda::gmove=1.3'
    // Closure (lazy) form: evaluated per task at run time, not once when this module is
    // parsed - see the note in nextflow.config on why that ordering matters here.
    publishDir { "${params.outdir}/gmove" }, mode: 'copy'

    input:
    path genome
    path protein_gff, stageAs: 'protein_evidence.gff3'
    path transcript_gff, stageAs: 'transcript_evidence.gff3'
    val cdsMin

    output:
    path "gmove_out.gff3"

    script:
    // NB: forcing distinct stageAs names above because the protein and transcript branches'
    // mergegff outputs both end up named identically ("exonerate_out_NR.gff3") - without
    // this, staging both into one gmove task collides ("input file name collision").
    // transcript_evidence.gff3 is a placeholder empty file (assets/NO_TRANSCRIPT_EVIDENCE)
    // when no transcript evidence was built for this run; the bash `-s` test below (real
    // file AND non-empty) decides whether to pass --rna, rather than checking the filename,
    // since stageAs renames it the same either way.
    // NB: I couldn't run gmove here to confirm --prot and --rna can be combined in one call;
    // the doc's examples only ever show them used separately. If gmove rejects using both at
    // once, split this into two `gmove` calls (one per evidence type) and `cat` their
    // gmove_out.gff3 outputs together instead.
    """
    RNA_FLAG=""
    [ -s transcript_evidence.gff3 ] && RNA_FLAG="--rna transcript_evidence.gff3"
    gmove -f ${genome} --prot protein_evidence.gff3 \$RNA_FLAG --cds ${cdsMin} -o gmove_out
    #NB: using cat rather than the doc's `cp .../filter/* dest`, which only works when
    #exactly one file lands in filter/ - cat is safe either way.
    cat gmove_out/filter/* > gmove_out.gff3
    """
}
