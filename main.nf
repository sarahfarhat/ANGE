#!/usr/bin/env nextflow


include { cleanproteindb }                     from "${baseDir}/modules/cleanproteindb.nf"
include { cleanproteindb as cleantranscripts }  from "${baseDir}/modules/cleanproteindb.nf"
include { clusterprot }                        from "${baseDir}/modules/clusterprot.nf"
include { repbasefilter }                      from "${baseDir}/modules/repbasefilter.nf"
include { blastprot }                          from "${baseDir}/modules/blastprot.nf"
include { filterblast as filterblastprot }     from "${baseDir}/modules/filterblast.nf"
include { fx2tab as fx2tab_genome }            from "${baseDir}/modules/fx2tab.nf"
include { fx2tab as fx2tab_dataset }           from "${baseDir}/modules/fx2tab.nf"
include { fx2tab as fx2tab_transcripts }       from "${baseDir}/modules/fx2tab.nf"
include { extractseq as extractseqprot }       from "${baseDir}/modules/extractseq.nf"
include { exonerate as exonerateprot }         from "${baseDir}/modules/exonerate.nf"
include { exonerate2gff as exonerate2gffprot } from "${baseDir}/modules/exonerate2gff.nf"
include { mergegff as mergegffprot }           from "${baseDir}/modules/mergegff.nf"
include { minimap2index }                      from "${baseDir}/modules/minimap2index.nf"
include { minimap2align }                      from "${baseDir}/modules/minimap2align.nf"
include { samtoolssort }                       from "${baseDir}/modules/samtoolssort.nf"
include { gffreadconvert }                     from "${baseDir}/modules/gffreadconvert.nf"
include { featurecounts }                      from "${baseDir}/modules/featurecounts.nf"
include { lowcountthreshold }                  from "${baseDir}/modules/lowcountthreshold.nf"
include { gmove }                              from "${baseDir}/modules/gmove.nf"
include { extractfasta }                       from "${baseDir}/modules/extractfasta.nf"
include { extractgenefasta as extractgenefastaprot }  from "${baseDir}/modules/extractgenefasta.nf"
include { extractgenefasta as extractgenefastatrans } from "${baseDir}/modules/extractgenefasta.nf"
include { buscodownload }                     from "${baseDir}/modules/buscodownload.nf"
include { buscoassess as buscoassessgenome }   from "${baseDir}/modules/buscoassess.nf"
include { buscoassess as buscoassessprot }     from "${baseDir}/modules/buscoassess.nf"
include { buscoassess as buscoassesstrans }    from "${baseDir}/modules/buscoassess.nf"
include { buscoassess as buscoassessfinal }    from "${baseDir}/modules/buscoassess.nf"
// The two "original file" busco checks below run on params.proteindatabase/transcriptomes
// exactly as given - a sanity check on the input evidence itself, independent of the genome
// or this pipeline's own gene-finding - see the BUSCO section of the workflow below.
include { buscoassess as buscoassessprotorig }  from "${baseDir}/modules/buscoassess.nf"
include { buscoassess as buscoassesstransorig } from "${baseDir}/modules/buscoassess.nf"
include { buscosummary }                       from "${baseDir}/modules/buscosummary.nf"
include { translateprot }                      from "${baseDir}/modules/translateprot.nf"

// Concatenates per-source gff3 files into one labeled evidence file, right before gmove/BUSCO
// need a single file - see combineevidence.nf. Aliased three ways: the final protein-branch
// evidence, the final transcript-branch evidence, and (only when run_low_count_filter is on)
// the pre-mergegff protein-branch exonerate output that step's read-count filter needs.
include { combineevidence as combineProteinEvidence }    from "${baseDir}/modules/combineevidence.nf"
include { combineevidence as combineTranscriptEvidence } from "${baseDir}/modules/combineevidence.nf"
include { combineevidence as combineProteinExonerate }   from "${baseDir}/modules/combineevidence.nf"

// Transcript-evidence branch (mirrors the protein branch above with blastn/cd-hit-est/est2genome)
include { clustertranscript }                  from "${baseDir}/modules/clustertranscript.nf"
include { blasttranscript }                    from "${baseDir}/modules/blasttranscript.nf"
include { filterblast as filterblasttrans }    from "${baseDir}/modules/filterblast.nf"
include { extractseq as extractseqtrans }      from "${baseDir}/modules/extractseq.nf"
include { exonerate as exoneratetrans }        from "${baseDir}/modules/exonerate.nf"
include { exonerate2gff as exonerate2gfftrans } from "${baseDir}/modules/exonerate2gff.nf"
include { mergegff as mergegfftrans }          from "${baseDir}/modules/mergegff.nf"

// NB: params.outdir's default lives in nextflow.config, not here - see the comment there.
// Setting it this late (after the `include` statements above, which is when each module's
// directives - including publishDir - get parsed) was the cause of the "publishDir path
// contains a variable with a null value" bug / the stray "null/" output directory.
// NB: chunkSize's default lives in nextflow.config (params{ chunkSize=... }), not here -
// same reasoning as outdir above. A second `params.chunkSize = ...` line used to sit here
// too and silently overrode the config value on every run (99000 instead of the 10000 in
// nextflow.config) - harmless while -profile test was always overriding it anyway, but
// exactly the kind of thing to catch before a real production run.
params.nb_cpus        = 28
params.idprot          = 0.95
params.idpercentprot   = 50
params.topHits         = 3

// Protein-evidence alignment (doc step 7, peptide branch)
params.exoneratePercent = 50
params.gmoveCdsMin      = 30

// Transcript-evidence branch: set params.transcriptomes to a glob (e.g.
// "test_data/*.fasta") or a comma-separated list of assembled-transcript fasta files to
// turn this on. Mirrors the protein branch (blastn instead of tblastn, cd-hit-est instead
// of cd-hit, exonerate est2genome instead of protein2genome) and its output is fed into
// gmove as --rna evidence alongside the protein branch's --prot evidence.
params.transcriptomes        = ""
params.idtranscript          = 0.99
params.idpercenttranscript   = 90
params.exoneratePercentTranscript = 80

// Optional step (doc step 3, "TBD" in the doc): skipped unless params.repbase_db is set.
params.repbase_db     = null
params.repbaseEvalue  = '1e-10'

// Optional step (doc step 9): needs raw paired-end RNA-seq reads (fastq), NOT the
// assembled transcripts used by the transcript-evidence branch above - this is a separate,
// independent way of using RNA-seq data (flagging poorly-supported genes for removal
// rather than feeding gmove directly). Skipped unless explicitly turned on.
params.run_low_count_filter = false
params.reads1                = null
params.reads2                = null
params.lowCountThreshold     = 10

// Optional step (doc step 12): skipped unless explicitly turned on. Set --busco_lineage to
// any other BUSCO dataset name (e.g. 'metazoa_odb10') to override the default.
params.run_busco     = false
params.busco_lineage = 'mollusca_odb12'

// Derives a filesystem-safe, human-readable "source id" from an input fasta path - used to
// tag which of several protein-database / transcriptome-assembly files each piece of
// evidence came from, all the way through to the final gff (see the NB in mergegff.nf).
// Strips a trailing .gz first, then the remaining last extension - "uniprotkb_D1.fasta" ->
// "uniprotkb_D1", "spades_mimit01_transcripts_renamed.fa" -> "spades_mimit01_transcripts_renamed".
def sourceId(f) {
    def n = f.getFileName().toString()
    if (n.endsWith('.gz')) { n = n[0..-4] }
    def dot = n.lastIndexOf('.')
    return dot > 0 ? n[0..<dot] : n
}

log.info """\
        ANGE - ANnotation of GEnes
        ===========================
        Genome (unmasked) : ${params.genome}
        Genome (masked)   : ${params.genome_msk}
        species           : ${params.species}
        Protein database  : ${params.proteindatabase}
        protein clustering id rate : ${params.idprot}
        Transcriptomes    : ${params.transcriptomes ?: '(none - protein evidence only)'}
        Number of sequences per chunk : ${params.chunkSize}
        CPUs number       : ${params.nb_cpus}
        REPbase filter    : ${params.repbase_db ? 'on (' + params.repbase_db + ')' : 'off'}
        Low-count filter  : ${params.run_low_count_filter}
        Busco assessment  : ${params.run_busco}
        Results directory : ${params.outdir}
        """
        .stripIndent()

workflow {
        /* Treat proteins: params.proteindatabase can be a glob (e.g. "db/*.fasta") or a
           comma-separated list of protein fasta files. Each file is carried through its own
           independent branch (clean -> cluster -> [repbase] -> blast -> extract -> exonerate
           -> mergegff), tagged with its own source id (see sourceId() above) - so every
           gene_models.<id>.gff3 published under results/gene_models/ is traceable back to
           exactly one input database file, with NO cross-file redundancy removal (that's the
           deliberate trade-off of this design: full traceability over the compute savings
           cd-hit would otherwise get by deduping near-identical sequences across sources).
           All per-source protein evidence is only combined into one file right before gmove
           (see combineProteinEvidence below), and even then every feature still carries an
           evidence_source=<id> attribute (added in mergegff.nf) so origin is recoverable
           from the merged file too. */
        Channel.fromPath(params.proteindatabase.split(',').collect { it.trim() }, checkIfExists: true)
               .map { f -> tuple(sourceId(f), f) }
               .set { protein_sources_ch }

        cleanproteindb(protein_sources_ch)
        clusterprot(cleanproteindb.out, params.idprot)

        /* Optional: drop sequences that match a repeat/REPbase protein database */
        def proteins_final_ch
        if (params.repbase_db) {
            repbasefilter(clusterprot.out, file(params.repbase_db), params.repbaseEvalue)
            proteins_final_ch = repbasefilter.out
        } else {
            proteins_final_ch = clusterprot.out
        }

        /* index genome and blast each source's proteins onto it, chunked for parallelism
           within a source - elem:1 splits only the fasta element of the (id, fasta) tuple,
           tagging every resulting chunk with its source's id. */
        proteins_final_ch
                   .splitFasta(by: params.chunkSize, file: true, elem: 1)
                   .set { prot_ch }

        blastprot(params.genome_msk, prot_ch)
        blastprot.out
                   .groupTuple()
                   .set { protein_blast_grouped_ch }
        filterblastprot(params.genome_msk, protein_blast_grouped_ch, params.idpercentprot, params.topHits)

        /* Extract genome regions and their match sequences, per source. genome_tab is the
           same genome-wide table for every source (id discarded right after this call - it
           has no "source" of its own); fx2tab_dataset's and filterblastprot's per-source
           outputs are explicitly joined by id (not just paired by arrival order, since each
           source finishes its own clean/cluster/blast/filter chain at a different pace). */
        fx2tab_genome(tuple('genome', file(params.genome)))
        fx2tab_genome.out
                   .map { id, tab -> tab }
                   .set { genome_tab_ch }

        fx2tab_dataset(proteins_final_ch)
        fx2tab_dataset.out
                   .join(filterblastprot.out)
                   .set { extractseqprot_in_ch }
        extractseqprot(genome_tab_ch, extractseqprot_in_ch, file("${baseDir}/modules/extractseq_v4.py"))

        /* Align blast matches to the genome, per source */
        exonerateprot(extractseqprot.out, 'protein2genome', params.exoneratePercent)
        exonerate2gffprot(exonerateprot.out.exo_dir)

        /* Optional: flag poorly RNA-seq-supported genes for removal (needs raw fastq reads).
           NB: computed once from ALL protein-branch sources combined (not per source) - RNA-
           seq read support for a candidate gene doesn't depend on which protein database
           first suggested it, so there is no reason to run this once per source. */
        def to_remove
        if (params.run_low_count_filter) {
            if (!params.reads1 || !params.reads2) {
                error "params.run_low_count_filter is true but params.reads1/params.reads2 are not set"
            }
            minimap2index(file(params.genome))
            minimap2align(minimap2index.out, file(params.reads1), file(params.reads2))
            samtoolssort(minimap2align.out)

            exonerate2gffprot.out
                       .map { id, gff -> gff }
                       .collect()
                       .set { protein_exonerate_gff_ch }
            combineProteinExonerate(protein_exonerate_gff_ch, 'protein_exonerate_combined')
            gffreadconvert(combineProteinExonerate.out)
            featurecounts(gffreadconvert.out, samtoolssort.out.bam, samtoolssort.out.bai)
            lowcountthreshold(featurecounts.out, params.lowCountThreshold)
            to_remove = lowcountthreshold.out
        } else {
            to_remove = file("${baseDir}/assets/empty_remove_list.txt")
        }

        /* Filter result per source: merge monoexonic genes and (optionally) drop low-count
           IDs - each output is gene_models.<id>.gff3, tagged evidence_source=<id> internally
           (see mergegff.nf). to_remove is a single broadcast value reused for every source. */
        mergegffprot(exonerate2gffprot.out, to_remove, file("${baseDir}/modules/gff_parser.py"))

        /* Combine every protein-branch source into one file for gmove/BUSCO - each feature
           still carries its evidence_source=<id> attribute (see mergegff.nf), so origin
           survives this merge too. */
        mergegffprot.out
                   .map { id, gff -> gff }
                   .collect()
                   .set { protein_gff_files_ch }
        combineProteinEvidence(protein_gff_files_ch, 'protein_evidence')

        /* Transcript-evidence branch: same per-source shape as the protein branch above,
           using assembled transcripts mapped onto the genome as additional evidence for
           gmove. params.transcriptomes can be a glob or a comma-separated list, same as
           params.proteindatabase above - each file gets its own independent branch. */
        def transcript_gff_for_gmove
        if (params.transcriptomes) {
            Channel.fromPath(params.transcriptomes.split(',').collect { it.trim() }, checkIfExists: true)
                   .map { f -> tuple(sourceId(f), f) }
                   .set { transcript_sources_ch }

            cleantranscripts(transcript_sources_ch)
            clustertranscript(cleantranscripts.out, params.idtranscript)

            clustertranscript.out
                   .splitFasta(by: params.chunkSize, file: true, elem: 1)
                   .set { trans_ch }

            blasttranscript(params.genome_msk, trans_ch)
            blasttranscript.out
                       .groupTuple()
                       .set { transcript_blast_grouped_ch }
            filterblasttrans(params.genome_msk, transcript_blast_grouped_ch, params.idpercenttranscript, params.topHits)

            fx2tab_transcripts(clustertranscript.out)
            fx2tab_transcripts.out
                       .join(filterblasttrans.out)
                       .set { extractseqtrans_in_ch }
            extractseqtrans(genome_tab_ch, extractseqtrans_in_ch, file("${baseDir}/modules/extractseq_v4.py"))

            exoneratetrans(extractseqtrans.out, 'est2genome', params.exoneratePercentTranscript)
            exonerate2gfftrans(exoneratetrans.out.exo_dir)
            mergegfftrans(exonerate2gfftrans.out, file("${baseDir}/assets/empty_remove_list.txt"), file("${baseDir}/modules/gff_parser.py"))

            mergegfftrans.out
                       .map { id, gff -> gff }
                       .collect()
                       .set { transcript_gff_files_ch }
            combineTranscriptEvidence(transcript_gff_files_ch, 'transcript_evidence')

            transcript_gff_for_gmove = combineTranscriptEvidence.out
        } else {
            transcript_gff_for_gmove = file("${baseDir}/assets/NO_TRANSCRIPT_EVIDENCE")
        }

        /* Run gmove to get the final gene models, combining protein + transcript evidence */
        gmove(file(params.genome), combineProteinEvidence.out, transcript_gff_for_gmove, params.gmoveCdsMin)

        /* Final gene models as fasta: spliced CDS (nucleotide) and translated protein.
           Always produced (not gated on run_busco) - this is the main sequence-level output
           people actually want out of the final gene set. */
        extractfasta(file(params.genome), gmove.out)

        /* Optional: Busco completeness assessment at several points along the pipeline - not
           just the final gene set, and not just the genome/pipeline's own output, but also
           the raw input evidence itself, so a low score can be traced to "the input database
           was already incomplete" vs. "the genome/pipeline lost something the input had".
           Every step below lands in the SAME results/busco/busco_summary.tsv (one row each). */
        if (params.run_busco) {
            // Downloaded once and reused by every buscoassess call below (each via a local
            // path rather than a bare lineage name) - see buscodownload.nf. Since this call's
            // only input is a plain value (not a channel), its output is automatically a
            // Nextflow "value channel" and safe to reference from multiple downstream
            // processes/parallel task instances, however many there end up being.
            buscodownload(params.busco_lineage)

            /* 1) Raw, unmasked genome assembly - a completeness baseline to compare
               everything else to. mode=genome (busco's own internal gene prediction over the
               whole assembly) - see buscoassess.nf for the exact command line. */
            buscoassessgenome(tuple('genome', file(params.genome)), 'genome', buscodownload.out)
            def busco_summary_ch = buscoassessgenome.out.summary

            /* 2) Each protein-database file exactly as given (params.proteindatabase),
               untouched by cleaning/clustering/blast/exonerate - does this input database
               contain a complete BUSCO set on its own, independent of the genome or this
               pipeline? Already real protein sequences, so mode=proteins directly - no
               extraction or translation needed. One run per file (protein_sources_ch is the
               same id-tagged channel the protein branch itself starts from, up at the top of
               this workflow). */
            protein_sources_ch
                       .map { id, fasta -> tuple("protein_original_${id}", fasta) }
                       .set { protein_original_busco_in_ch }
            buscoassessprotorig(protein_original_busco_in_ch, 'proteins', buscodownload.out)
            busco_summary_ch = busco_summary_ch.mix(buscoassessprotorig.out.summary)

            /* 3) Protein-branch evidence the genome/pipeline actually recovered: per source,
               extracted from exonerate2gffprot's raw gff (BEFORE mergegff's monoexonic merge
               and BEFORE combineevidence's cross-source concat), then translated to protein
               (see translateprot.nf) so this also runs mode=proteins - directly comparable to
               (2) above. extractgenefasta.nf's own note covers why this extraction is only a
               rough proxy (raw genomic span, "gene" and "exon" rows both extracted); on top
               of that, translateprot.nf's note covers why the translation itself is only a
               frame-1 guess, not a real spliced/correctly-framed CDS translation. */
            extractgenefastaprot(genome_tab_ch, exonerate2gffprot.out, file("${baseDir}/modules/extractseq_genomev4.py"))
            translateprot(extractgenefastaprot.out)
            translateprot.out
                       .map { id, fasta -> tuple("protein_branch_${id}", fasta) }
                       .set { protein_busco_in_ch }
            buscoassessprot(protein_busco_in_ch, 'proteins', buscodownload.out)
            busco_summary_ch = busco_summary_ch.mix(buscoassessprot.out.summary)

            if (params.transcriptomes) {
                /* 4) Each transcriptome-assembly file exactly as given
                   (params.transcriptomes), untouched by cleaning/clustering/blast/exonerate -
                   same idea as (2), but these are already real assembled transcript
                   sequences (properly spliced, unlike the raw genomic spans in (5) below), so
                   mode=transcriptome is the textbook-correct BUSCO mode here. Re-globs the
                   same files params.transcriptomes matches at the top of this workflow,
                   rather than reusing that section's channel variable across the two
                   independent `if (params.transcriptomes)` blocks. */
                Channel.fromPath(params.transcriptomes.split(',').collect { it.trim() }, checkIfExists: true)
                           .map { f -> tuple("transcriptome_original_${sourceId(f)}", f) }
                           .set { transcript_original_busco_in_ch }
                buscoassesstransorig(transcript_original_busco_in_ch, 'transcriptome', buscodownload.out)
                busco_summary_ch = busco_summary_ch.mix(buscoassesstransorig.out.summary)

                /* 5) Transcript-branch evidence the genome/pipeline actually recovered: per
                   source, extracted from exonerate2gfftrans's raw gff (BEFORE mergegff) -
                   same treatment as (3) but NOT translated: stays mode=transcriptome, since
                   this extraction is still nucleotide (an unspliced genomic span, not a real
                   transcript - see extractgenefasta.nf). */
                extractgenefastatrans(genome_tab_ch, exonerate2gfftrans.out, file("${baseDir}/modules/extractseq_genomev4.py"))
                extractgenefastatrans.out
                           .map { id, fasta -> tuple("transcript_branch_${id}", fasta) }
                           .set { transcript_busco_in_ch }
                buscoassesstrans(transcript_busco_in_ch, 'transcriptome', buscodownload.out)
                busco_summary_ch = busco_summary_ch.mix(buscoassesstrans.out.summary)
            }

            /* 6) gmove's final gene models, mode=proteins on the translated aa fasta - the
               one that actually matters for judging the final annotation, since only
               gmove's own output has proper CDS/phase for a real translation. */
            buscoassessfinal(tuple('final_gmove', extractfasta.out.proteins), 'proteins', buscodownload.out)
            busco_summary_ch = busco_summary_ch.mix(buscoassessfinal.out.summary)

            buscosummary(busco_summary_ch.collect(), file("${baseDir}/modules/busco_summary.py"))
        }
}
