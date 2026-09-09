#!/usr/bin/env nextflow

process filterblast {
    // Closure (lazy) form: evaluated per task at run time, not once when this module is
    // parsed - see the note in nextflow.config on why that ordering matters here.
    publishDir { "${params.outdir}/candidate_regions" }, mode: 'copy'

    input:
    path genomelength
    path blastout
    val idpercentprot
    val topHits
    val branchLabel

    output:
    path "candidate_regions.${branchLabel}.gff"

    script:
    """
    #Get genome sequences length
    awk '{OFS="\t"; if(\$0~/^>/){if(NR!=1){print l,id;} id=\$1; l=0}else{l=l+length(\$1)}}END{print l,id}' genomelength | sed 's/>//' > genome.lg
    #concatenate all blast results
    cat ${blastout} > tblastn_results_proteins.out
    #Filter HSPs at idpercentprot % identity and convert to gff
    #tblastn columns: 1=qseqid(protein) 2=sseqid(genome) 3=pident 7=qstart 8=qend 9=sstart 10=send 12=bitscore
    awk -F"\t" -v OFS='\t' '{if(\$3>='${idpercentprot}'){if(\$10<\$9){strand="-"; tmp=\$9; \$9=\$10; \$10=tmp}else{strand="+"} print \$2,"BLAST","HSP",\$9,\$10,\$12,strand,".","HSP="\$1;}}' tblastn_results_proteins.out | sort -k9,9 -k1,1 -k4,4n > filter1.gff
    #Group blast HSPs into "genes" merging exons within 100kb, compute a score for each gene
    awk '{OFS="\t"; if(\$9!=id || \$1!=chr){if(NR!=1){print chr,"BLAST","HSP",deb,end,mscor/cpt,strand,".",id} deb=\$4; end=\$5; chr=\$1; strand=\$7; id=\$9; mscor=\$6; cpt=1}else{if(\$5<deb+100000){end=\$5; mscor=mscor+\$6; cpt=cpt+1}else{print chr,"BLAST","HSP",deb,end,mscor/cpt,strand,".",id; deb=\$4; end=\$5; chr=\$1; strand=\$7; id=\$9; mscor=\$6; cpt=1}}}END{print chr,"BLAST","HSP",deb,end,mscor/cpt,strand,".",id;}' filter1.gff > filter2.gff
    #Keep only the top N best scores matches per feature
    sort -k9,9 -k6,6nr filter2.gff | awk -v topn='${topHits}' '{if(id!=\$9){print \$0; cpt=1; id=\$9}else{if(cpt<topn){print \$0; cpt=cpt+1}}}' > filter3.gff
    #Add 5kb before and after the feature
    #NB: falls back to just padding by 5kb (no upper clamp) when a hit's chromosome isn't
    #found in genome.lg - the original version left `newend` blank in that case (an ID
    #mismatch between the length table and blast's sseqid), which crashed extractseq_v4.py
    #downstream on an unparseable coordinate. This also logs a warning so a real ID mismatch
    #is still visible instead of silently working around it.
    awk -F"\t" 'BEGIN{while(getline < "genome.lg" > 0){len[\$2]=\$1}}{OFS="\t"; chr=\$1; if((\$4-5000) < 0){newstart=1}else{newstart=\$4-5000} if(!(chr in len)){print "WARNING: no length found in genome.lg for chromosome:", chr, "- not clamping this feature end coordinate" > "/dev/stderr"; newend=\$5+5000}else if((\$5+5000)>len[chr]){newend=len[chr]}else{newend=\$5+5000} print \$1,\$2,\$3,newstart,newend,\$6,\$7,\$8,\$9}' filter3.gff > candidate_regions.${branchLabel}.gff
    """
}
