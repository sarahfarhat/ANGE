#!/usr/bin/env nextflow

// Optional step (doc step 3, left as "TBD" in the original doc - the doc has no concrete
// recipe here, so this design is mine: blastp the protein dataset against a user-supplied
// repeat protein database (e.g. RepBase peptide library, which is not publicly redistributable
// so you must provide your own copy via params.repbase_db) and drop any sequence with a
// significant hit before it goes into the genome blast.
process repbasefilter {
    container 'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd'
    input:
    path proteins_file
    path repbase_db
    val evalue

    output:
    path "proteins_norepeat.fasta"

    script:
    """
    makeblastdb -in ${repbase_db} -dbtype prot -out repbase_db
    blastp -query ${proteins_file} -db repbase_db -evalue ${evalue} -outfmt 6 -num_threads ${task.cpus} -out repbase_hits.out
    cut -f1 repbase_hits.out | sort -u > repbase_hits.list
    awk 'BEGIN{while((getline line < "repbase_hits.list")>0) bad[line]=1}
         /^>/{split(\$0,a," "); id=substr(a[1],2); keep = !(id in bad)}
         keep' ${proteins_file} > proteins_norepeat.fasta
    """
}
