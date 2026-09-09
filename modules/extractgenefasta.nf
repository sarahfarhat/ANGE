#!/usr/bin/env nextflow

process extractgenefasta {
    input:
    path genome_tab
    tuple val(id), path(gff)
    path 'extractseq_genomev4.py'

    output:
    tuple val(id), path("${id}.fasta")

    script:
    // Extracts each gene/mRNA/HSP/exon feature's raw genomic span as nucleotide sequence -
    // NOT spliced, NOT translated. This is a rough proxy, not a real coding sequence.
    // NB: fed exonerate2gff's *raw*, per-source output directly (before mergegff's monoexonic
    // merge, and before combineevidence's cross-source concat - see main.nf's BUSCO section),
    // so - unlike mergegff's gene_models.*.gff3, which (per the note this replaced) only ever
    // has "exon" rows - this input has BOTH "gene" and "exon" rows for every feature, and
    // extractseq_genomev4.py's type filter extracts both. That means a feature's full span
    // AND its individual exon spans both end up as separate sequences in the output fasta -
    // worth knowing when reading BUSCO's n_total for this step, since it isn't "one sequence
    // per candidate gene". Say if you'd rather this only count exon rows.
    // extractseq_genomev4.py always names its output after the input gff's basename
    // ("exonerate_out.fasta", since exonerate2gff.nf's output file is always literally named
    // "exonerate_out.gff3") - renamed to "${id}.fasta" below so parallel per-source runs
    // don't collide and BUSCO's output stays traceable to its source.
    """
    python3 extractseq_genomev4.py ${genome_tab} ${gff}
    mv ${gff.baseName}.fasta ${id}.fasta
    """
}
