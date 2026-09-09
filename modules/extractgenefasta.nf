#!/usr/bin/env nextflow

process extractgenefasta {
    input:
    path genome_tab
    path gff
    path 'extractseq_genomev4.py'

    output:
    path "${gff.baseName}.fasta"

    script:
    // Extracts each gene/mRNA/HSP/exon feature's raw genomic span as nucleotide sequence -
    // NOT spliced, NOT translated. This is a rough proxy, not a real coding sequence.
    // NB: gene_models.*.gff3 (this process's usual input, from mergegff/gff_parser.py) turns
    // out to only ever contain "exon"-type rows in practice - gff_parser.py groups rows by
    // the id it parses out of each attributes column, and a "gene" row's id (from its
    // `ID=EXONAME_x`) never matches its own exon rows' id (from `Parent=EXOTRANS1_x` - a
    // different prefix), so every "gene" row is a singleton and gets silently dropped; only
    // exon rows (which share their Parent id across a transcript's exons) survive. Harmless
    // for gmove, which wants spliced alignment evidence rather than pre-formed genes, but it
    // means extractseq_genomev4.py's type filter needs "exon" in it or this produces an
    // empty fasta (which is what made BUSCO fail with "does not contain nucleotide
    // sequences" on a real run - not a busco problem, an empty-input problem).
    // extractseq_genomev4.py names its output after the gff's basename (swapping the
    // extension for .fasta), hence the output path below matching that exactly.
    """
    python3 extractseq_genomev4.py ${genome_tab} ${gff}
    """
}
