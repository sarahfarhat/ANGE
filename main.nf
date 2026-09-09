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
include { buscosummary }                       from "${baseDir}/modules/buscosummary.nf"

// Transcript-evidence branch (mirrors the protein branch above with blastn/cd-hit-est/est2genome)
include { mergetranscripts }                   from "${baseDir}/modules/mergetranscripts.nf"
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
        /* Treat proteins : Clean protein database: rename sequences and remove redundancy */
        cleanproteindb(params.proteindatabase)
        clusterprot(cleanproteindb.out, params.idprot)

        /* Optional: drop sequences that match a repeat/REPbase protein database */
        def proteins_final
        if (params.repbase_db) {
            repbasefilter(clusterprot.out, file(params.repbase_db), params.repbaseEvalue)
            proteins_final = repbasefilter.out
        } else {
            proteins_final = clusterprot.out
        }

        proteins_final
                   .splitFasta(by: params.chunkSize, file:true)
                   .set { prot_ch }

        /* index genome and blast proteins onto it */
        blastprot(params.genome_msk, prot_ch)
        filterblastprot(params.genome_msk, blastprot.out.collect(), params.idpercentprot, params.topHits, 'protein')

        /* Extract genome regions and their match sequences */
        fx2tab_genome(file(params.genome))
        fx2tab_dataset(proteins_final)
        extractseqprot(fx2tab_genome.out, fx2tab_dataset.out, filterblastprot.out, file("${baseDir}/modules/extractseq_v4.py"))

        /* Align blast matches to the genome */
        exonerateprot(extractseqprot.out.query_dir, extractseqprot.out.genome_dir, 'protein2genome', params.exoneratePercent)
        exonerate2gffprot(exonerateprot.out.exo_dir)

        /* Optional: flag poorly RNA-seq-supported genes for removal (needs raw fastq reads) */
        def to_remove
        if (params.run_low_count_filter) {
            if (!params.reads1 || !params.reads2) {
                error "params.run_low_count_filter is true but params.reads1/params.reads2 are not set"
            }
            minimap2index(file(params.genome))
            minimap2align(minimap2index.out, file(params.reads1), file(params.reads2))
            samtoolssort(minimap2align.out)
            gffreadconvert(exonerate2gffprot.out)
            featurecounts(gffreadconvert.out, samtoolssort.out.bam, samtoolssort.out.bai)
            lowcountthreshold(featurecounts.out, params.lowCountThreshold)
            to_remove = lowcountthreshold.out
        } else {
            to_remove = file("${baseDir}/assets/empty_remove_list.txt")
        }

        /* Filter result: merge monoexonic genes and (optionally) drop low-count IDs */
        mergegffprot(exonerate2gffprot.out, to_remove, file("${baseDir}/modules/gff_parser.py"), 'protein')

        /* Transcript-evidence branch: same shape as the protein branch above, using
           assembled transcripts mapped onto the genome as additional evidence for gmove. */
        def transcript_gff_for_gmove
        if (params.transcriptomes) {
            Channel.fromPath(params.transcriptomes.split(',').collect { it.trim() }, checkIfExists: true)
                   .collect()
                   .set { transcript_files_ch }
            mergetranscripts(transcript_files_ch)
            cleantranscripts(mergetranscripts.out)
            clustertranscript(cleantranscripts.out, params.idtranscript)

            clustertranscript.out
                   .splitFasta(by: params.chunkSize, file:true)
                   .set { trans_ch }

            blasttranscript(params.genome_msk, trans_ch)
            filterblasttrans(params.genome_msk, blasttranscript.out.collect(), params.idpercenttranscript, params.topHits, 'transcript')

            fx2tab_transcripts(clustertranscript.out)
            extractseqtrans(fx2tab_genome.out, fx2tab_transcripts.out, filterblasttrans.out, file("${baseDir}/modules/extractseq_v4.py"))

            exoneratetrans(extractseqtrans.out.query_dir, extractseqtrans.out.genome_dir, 'est2genome', params.exoneratePercentTranscript)
            exonerate2gfftrans(exoneratetrans.out.exo_dir)
            mergegfftrans(exonerate2gfftrans.out, file("${baseDir}/assets/empty_remove_list.txt"), file("${baseDir}/modules/gff_parser.py"), 'transcript')

            transcript_gff_for_gmove = mergegfftrans.out
        } else {
            transcript_gff_for_gmove = file("${baseDir}/assets/NO_TRANSCRIPT_EVIDENCE")
        }

        /* Run gmove to get the final gene models, combining protein + transcript evidence */
        gmove(file(params.genome), mergegffprot.out, transcript_gff_for_gmove, params.gmoveCdsMin)

        /* Final gene models as fasta: spliced CDS (nucleotide) and translated protein.
           Always produced (not gated on run_busco) - this is the main sequence-level output
           people actually want out of the final gene set. */
        extractfasta(file(params.genome), gmove.out)

        /* Optional: Busco completeness assessment at several points along the pipeline, not
           just the final gene set - so you can see whether/where completeness improves.
             - genome: raw assembly, mode=genome - a baseline to compare everything else to.
             - final_gmove: gmove's final gene models, mode=proteins on the translated aa
               fasta above - the one that actually matters for judging the annotation.
             - protein_branch / transcript_branch: the pre-gmove exonerate-derived gene
               models from each evidence branch, mode=transcriptome on a raw (unspliced,
               untranslated) per-feature sequence extraction - see extractgenefasta.nf for
               why these are only an approximate signal, not a real proteome. */
        if (params.run_busco) {
            // Downloaded once and reused by all four calls below (each via a local path
            // rather than a bare lineage name) - see buscodownload.nf. Since this call's only
            // input is a plain value (not a channel), its output is automatically a Nextflow
            // "value channel" and safe to reference from multiple downstream processes.
            buscodownload(params.busco_lineage)

            buscoassessgenome(file(params.genome), 'genome', buscodownload.out, 'genome')
            buscoassessfinal(extractfasta.out.proteins, 'proteins', buscodownload.out, 'final_gmove')

            extractgenefastaprot(fx2tab_genome.out, mergegffprot.out, file("${baseDir}/modules/extractseq_genomev4.py"))
            buscoassessprot(extractgenefastaprot.out, 'transcriptome', buscodownload.out, 'protein_branch_gene_models')

            def busco_summary_ch = buscoassessgenome.out.summary
                    .mix(buscoassessfinal.out.summary, buscoassessprot.out.summary)

            if (params.transcriptomes) {
                extractgenefastatrans(fx2tab_genome.out, mergegfftrans.out, file("${baseDir}/modules/extractseq_genomev4.py"))
                buscoassesstrans(extractgenefastatrans.out, 'transcriptome', buscodownload.out, 'transcript_branch_gene_models')
                busco_summary_ch = busco_summary_ch.mix(buscoassesstrans.out.summary)
            }

            buscosummary(busco_summary_ch.collect(), file("${baseDir}/modules/busco_summary.py"))
        }
}
