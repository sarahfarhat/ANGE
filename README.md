# ANGE - ANnotation of GEnes

A Nextflow (DSL2) pipeline that predicts gene models from a genome assembly plus a
protein dataset and, optionally, assembled transcripts: clean and de-duplicate each
dataset, blast it onto the genome, align matches with exonerate, merge/clean the
resulting gene models, and run [gmove](https://github.com/institut-de-genomique/Gmove)
to produce final gene predictions - combining protein and transcript evidence when both
are available - with optional REPbase filtering, RNA-seq-based low-confidence gene
removal, and Busco completeness assessment on evidences.

## Pipeline steps

**Protein-evidence branch** (always runs):

1. **cleanproteindb** - rename dataset sequences (strip characters that break downstream tools).
2. **clusterprot** - reduce redundancy with cd-hit.
3. **repbasefilter** *(optional)* - drop sequences that blast-hit a repeat protein database.
4. **blastprot** - `tblastn` the protein dataset against the masked genome (chunked for speed).
5. **filterblast** - filter hits by %identity, merge nearby HSPs into candidate gene regions, keep the top N hits per protein, pad each region by 5kb.
6. **fx2tab / extractseq** - pull out the matching genome region and query sequence for every candidate.
7. **exonerate** (`protein2genome`) - align each query back to its genome region.
8. **exonerate2gff** - convert exonerate's target-gff output into a proper gff3 with absolute genome coordinates.
9. **mergegff** - merge overlapping monoexonic genes.

**Transcript-evidence branch** *(optional - on when `params.transcriptomes` is set)*, same
shape as above but for assembled transcripts: **mergetranscripts** (concatenate your
per-sample fasta files) → **cleantranscripts** → **clustertranscript** (`cd-hit-est`) →
**blasttranscript** (`blastn`) → **filterblast** → **fx2tab / extractseq** → **exonerate**
(`est2genome`) → **exonerate2gff** → **mergegff**.

**Then, regardless of the transcript branch:**

- **minimap2 / samtools / gffread / featureCounts** *(optional, independent of the
  transcript-evidence branch above - see the params table)* - flag genes with weak
  RNA-seq *read* support, using raw paired-end fastq rather than assembled transcripts.
- **gmove** - final gene prediction, combining `--prot` evidence (always) with `--rna`
  evidence (when the transcript branch ran).
- **extractfasta** - pull the final gene models' CDS (nucleotide) and translated protein
  sequences out of gmove's gff3, via `gffread`.
- **buscodownload** *(optional)* - pre-downloads the busco lineage dataset once, shared by
  all four busco steps below via a local path rather than each downloading its own copy by
  name (see "Known limitations" - getting this actually working end-to-end took two fixes,
  after redundant downloads caused a real run to fail with "No space left on device").
- **buscoassess** *(optional)* - completeness check at several points: the raw genome, each
  branch's pre-gmove gene models, and (the one that matters most) the final gene set's
  protein sequences - see the Output section below for what "per step" means here.
- **buscosummary** *(optional)* - combines all of the above busco runs into one table.

## Requirements

- Nextflow (`module load nextflow` or see https://www.nextflow.io/)
- A container engine: Singularity/Apptainer on the HPC (`-profile singularity`), or
  Docker locally (`-profile test` turns this on - Docker Desktop, or Colima on older
  macOS that Docker Desktop won't install on)
- Internet access from the execution host the first time each container is needed
  (most tools use pre-built images from `community.wave.seqera.io` / `quay.io`; exonerate
  and gmove have no pre-existing public container, so they're built on the fly from a
  `conda` spec via Nextflow's Wave integration - see "Known limitations" below)

## Running it

```bash
# on your HPC (see run_ANGE for the SLURM wrapper) - uses Singularity
nextflow run main.nf -profile singularity -resume

# a small local test (MacOS) run once you've added sample data under test_data/mini/ (see below) - uses Docker
nextflow run main.nf -profile test -resume
```

No container engine is enabled by default - it's only ever turned on by a profile (see
`nextflow.config`), so each run picks exactly one engine via `-profile`.

## Output

Everything Nextflow actually runs lives under `work/` (named by task hash, not human-
readable) - that's normal Nextflow behavior and not where you should go looking. The
pipeline copies the results worth keeping into `${params.outdir}` (default `results/`,
relative to wherever you launch `nextflow run` from) as it runs:

| Path | Contents |
|---|---|
| `results/candidate_regions/candidate_regions.protein.gff` | Blast-derived candidate regions, protein branch |
| `results/candidate_regions/candidate_regions.transcript.gff` | Same, transcript branch (only if `transcriptomes` is set) |
| `results/gene_models/gene_models.protein.gff3` | Exonerate-derived gene models feeding gmove's `--prot`, protein branch |
| `results/gene_models/gene_models.transcript.gff3` | Same, transcript branch (only if `transcriptomes` is set) |
| `results/gmove/gmove_out.gff3` | **The final gene predictions** - this is the main output |
| `results/gmove/gmove_final.cds.fasta` | Final gene models, spliced coding (nucleotide) sequence |
| `results/gmove/gmove_final.proteins.fasta` | Final gene models, translated protein sequence |
| `results/busco/genome_busco/` | Busco on the raw genome assembly (`-m genome`) - a completeness baseline (only if `run_busco` is true) |
| `results/busco/protein_branch_gene_models_busco/` | Busco on the protein branch's pre-gmove gene models (only if `run_busco` is true) |
| `results/busco/transcript_branch_gene_models_busco/` | Same, transcript branch (only if `run_busco` and `transcriptomes` are both set) |
| `results/busco/final_gmove_busco/` | Busco on the **final** gene models' protein sequences (`-m proteins`) - the number that matters most (only if `run_busco` is true) |
| `results/busco/busco_summary.tsv` | All of the above combined into one table - `step`, `lineage`, `complete_pct`, `single_copy_pct`, `duplicated_pct`, `fragmented_pct`, `missing_pct`, `n_total` (only if `run_busco` is true) |

Set `--outdir some/other/path` to change where these land. (If you already have a stray
`null/` directory from before this fix - `rm -rf null` to clean it up; nothing under it is
needed, it's the exact same content that now correctly lands under `results/`.)

A note on what "per step" busco means here: only gmove's own gff3 output carries proper
`CDS`/phase annotations, so only `final_gmove` is scored on a real, correctly-translated
protein set (`-m proteins`) - which is also why that's the number to trust most. The
`genome`, `protein_branch_gene_models` and `transcript_branch_gene_models` runs use
`-m genome` / `-m transcriptome` on raw, unspliced sequence extractions instead (see
extractgenefasta.nf) - useful as a rough trend across the pipeline, but not directly
comparable to a "real" transcriptome/proteome BUSCO score.

Core parameters (set in `nextflow.config` or with `--paramName value`):

| Param | Default | Meaning |
|---|---|---|
| `outdir` | `results` | Where published output files are copied to |
| `genome` | - | Unmasked genome fasta (used for sequence extraction, gmove, Busco) |
| `genome_msk` | - | Masked genome fasta (used for the blast database) |
| `proteindatabase` | - | Protein (or transcript) dataset fasta, optionally gzipped |
| `species` | - | Label used in log output |
| `nb_cpus` | 28 | CPUs given to every process |
| `idprot` | 0.95 | cd-hit clustering identity |
| `idpercentprot` | 50 | Minimum %identity kept from blast |
| `topHits` | 3 | Best-scoring regions kept per protein |
| `chunkSize` | 10000 | Sequences per blast chunk (both branches) |
| `exoneratePercent` | 50 | exonerate `--percent` threshold, protein branch |
| `gmoveCdsMin` | 30 | gmove `--cds` minimum |

Transcript-evidence branch (off unless `transcriptomes` is set):

| Param | Default | Meaning |
|---|---|---|
| `transcriptomes` | `""` | Glob (e.g. `"test_data/*.fasta"`) or comma-separated list of assembled-transcript fasta files. Setting this turns the whole branch on. |
| `idtranscript` | 0.99 | cd-hit-est clustering identity |
| `idpercenttranscript` | 90 | Minimum %identity kept from the transcript blast |
| `exoneratePercentTranscript` | 80 | exonerate `--percent` threshold, transcript branch |

Other optional steps (all off by default):

| Param | Meaning |
|---|---|
| `repbase_db` | Path to a repeat protein fasta (e.g. your own RepBase copy - it's not publicly redistributable, so you supply it). Setting this turns the step on. |
| `run_low_count_filter` | Set `true` + `reads1`/`reads2` (paired-end RNA-seq **fastq reads**, not the assembled transcripts above) to flag and drop genes with weak read support. Independent of the transcript-evidence branch. |
| `lowCountThreshold` | 10 - genes with fewer than this many reads are dropped. |
| `run_busco` | Set `true` to run Busco at each step (genome, each branch's gene models, final gene set) and produce `busco_summary.tsv`. |
| `busco_lineage` | `mollusca_odb12` by default - set `--busco_lineage some_other_odb12` (or any other BUSCO dataset name) to use a different lineage. Applied to every busco step. |

## Test data

`nextflow.config` includes a `test` profile pointing at `test_data/mini/test_genome.fasta`,
`test_data/mini/test_genome_masked.fasta`, `test_data/mini/test_proteins.fasta` and
`test_data/mini/test_transcripts.fasta` - a small, checked-in subset kept well under
GitHub's 100MB file size limit so `-profile test` works out of the box on a fresh clone.

## Fixes and design calls made while wiring this up

The original design doc (made by Giacomo Chiappa) had draft awk/shell snippets with a handful of small bugs. All
noted inline in the relevant module, summarized here:

- **filterblast**: the genome chromosome-length calculation was missing its `END` block,
  so the last sequence in the genome fasta was silently dropped from the length table
  (its region got clamped to a 0-length feature). Fixed.
- **blastprot**: switched from `blastx` (genome as query) to `tblastn` (protein as query,
  matching the doc) - `filterblast`'s column mapping was updated to match.
- **exonerate2gff**: the coordinate remap (`region_start + local_position`) was off by
  one; it needed `- 1` since local position 1 maps onto the region's start position.
  Also replaced the doc's hardcoded `^Scaffold|^contig` header filter (specific to one
  genome's naming) with a filter on the `chr:start-end` region-id format the pipeline
  itself generates, so it isn't tied to any particular assembly's sequence names.
- **lowcountthreshold**: the doc's draft (`$7>10`) would have flagged well-supported
  genes for removal under a step titled "list low count genes". Flipped to `$7<threshold`
  to match the step's stated intent - genes get removed for having *too few* reads.
- **process.cpus**: `params.nb_cpus` was previously only printed in the log banner and
  never actually applied - every process ran at the default of 1 cpu. Now wired into
  `nextflow.config`.
- **publishDir resolving to `null/`**: `params.outdir`'s default was originally set with a
  plain `params.outdir = "results"` line in `main.nf`. That line runs too late: the
  `include` statements above it are what parse each module's process definitions -
  including any plain (non-closure) `"${params.outdir}/..."` publishDir string - so every
  publishDir had already locked in `outdir == null` before that default was set, and
  everything got copied to a literal `null/` directory. Fixed by (a) moving the default
  into `nextflow.config`, which is fully loaded before any module is parsed, and (b)
  wrapping every `publishDir` in a closure (`publishDir { "${params.outdir}/x" }`) so it's
  evaluated lazily per task instead of once at parse time.
- **REPbase filtering**: the doc left this step as "TBD" with no concrete recipe. It's
  implemented as a blastp of the dataset against a user-supplied repeat protein database,
  dropping anything with a significant hit - swap in a different approach if you had
  something else in mind.
- **per-branch busco fed an empty fasta**: `gene_models.*.gff3` (mergegff's output) turns
  out to only ever contain `exon`-type rows in practice, never `gene`/`mRNA` - see the note
  in `extractgenefasta.nf` for why (`gff_parser.py` drops every `gene` row as a singleton,
  since its id never matches its own exon rows' id). Harmless for gmove, which wants exon-
  level alignment evidence anyway, but `extractseq_genomev4.py`'s type filter didn't
  originally include `exon`, so it silently wrote an empty fasta for the protein/transcript
  branch busco steps, which busco then rejected as "does not contain nucleotide sequences".
  Fixed by adding `exon` to that filter.

## Known limitations / things to double-check on your first real run

- **exonerate2gff's attribute parsing** (splitting exonerate's `--showtargetgff`
  attribute string into `t[4]`, `t[7]`, `t[10]`...) was ported as-is from the doc since I
  couldn't run exonerate here to confirm the exact token layout. Check `exonerate_out.gff3`
  after your first real run - if the attribute fields look shifted or wrong, send me an
  example `.exo` file and I'll fix the indices.
- **busco lineage downloaded 4x (two bugs, not one)**: the first version of this feature had
  all four `buscoassess*` calls independently download+decompress their own copy of the
  lineage dataset in their own isolated task work dir. That's both slow and, worse, means
  all 4 copies stick around permanently under `work/` (Nextflow doesn't clean up completed
  task dirs) - which is what caused a real run to fail with "No space left on device"
  partway through. Fixed with a new `buscodownload` step that fetches the dataset once via
  `busco --download`, with all four `buscoassess` calls pointed at that local directory
  instead of a bare lineage name. That alone wasn't enough, though: Nextflow stages the
  downloaded directory into each consuming task using just its basename (e.g.
  `mollusca_odb12`), which is textually identical to a bare lineage name once interpolated
  into the busco command - busco's local-vs-download logic for `-l` appears to key off
  whether the argument looks like a path (contains a `/`), so it was still silently falling
  through to its own auto-download every time. Fixed by passing `-l ./mollusca_odb12`
  (prefixed with `./`) instead of `-l mollusca_odb12`, forcing busco to recognize it as the
  already-staged local directory.
- **exonerate and gmove containers**: no pre-built public container image could be found
  for either, so both use Nextflow's Wave "build from conda" integration instead of a
  pinned image hash like the other modules. This needs outbound internet from wherever
  Nextflow runs (fine on an HPC login/compute node with normal internet, less fine on a
  fully air-gapped cluster - let me know if that's your setup and I'll switch these to a
  locally-built container instead).
- **gmove with combined `--prot` + `--rna` evidence**: the doc's examples only ever show
  gmove called with one evidence type at a time, so I couldn't confirm it accepts both
  flags in a single call. If a real run shows gmove rejecting that, say so and I'll split
  it into two `gmove` calls (one per evidence type) with their outputs concatenated instead.
- **extractfasta assumes gmove's gff3 has proper CDS/phase**: I couldn't run gmove here to
  confirm this. If `gmove_final.proteins.fasta` / `gmove_final.cds.fasta` come out empty or
  look wrong, send me the head of `gmove_out.gff3` and I'll adjust the extraction.
- **busco_summary.py's parsing of busco's short_summary text**: written against BUSCO's
  standard, long-stable one-line summary format, but I couldn't run busco here to confirm
  the exact wording for v6.1.0. If `busco_summary.tsv` comes out with WARNING lines about
  unparseable summaries (check the pipeline's stdout) or missing rows, send me one of the
  `*_busco/short_summary.specific.*.txt` files and I'll fix the regex.
- **genome-mode busco runtime**: `-m genome` runs its own internal gene prediction (metaeuk
  by default) over the whole assembly, which can be considerably slower than the
  transcriptome/protein modes on a real-size genome - budget more walltime/resources for
  that step than the others once you move past the tiny test data.
- This pipeline was built and reviewed without being able to actually execute Nextflow in
  this environment (network restrictions), so please run it and send me any errors -
  they're usually quick fixes.
