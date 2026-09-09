#!/usr/bin/env nextflow

// Translates the protein branch's per-source, pre-mergegff extracted genomic sequence (see
// extractgenefasta.nf) to protein, in frame 1 with the standard genetic code, so that
// source's BUSCO check can run in -m proteins mode instead of -m transcriptome.
//
// NB: the input here is still an unspliced, raw genomic span (see extractgenefasta.nf's own
// note) - it includes intron sequence for anything but a single-exon (intron-free) hit, and
// frame 1 is only a guess (the real reading frame isn't known at this point in the pipeline,
// unlike gmove's own output which carries proper CDS/phase - see extractfasta.nf). So a
// multi-exon gene's translation will typically run into frameshifts and premature stop
// codons partway through. This is a rough proxy, same spirit as the rest of this per-branch
// BUSCO signal, not a substitute for gmove's real, spliced, correctly-framed translation.
// -x (--allow-unknown-codon) translates any codon seqkit can't resolve (e.g. containing an N)
// to 'X' instead of erroring - verified locally: `seqkit translate -f 1 -T 1 -x` on a small
// test fasta, including one sequence with an ambiguous base, produced the expected
// frame-1/standard-code translation without error.
process translateprot {
    container 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'
    input:
    tuple val(id), path(fasta)

    output:
    tuple val(id), path("${id}.protein.fasta")

    script:
    """
    seqkit translate -f 1 -T 1 -x ${fasta} -o ${id}.protein.fasta
    """
}
