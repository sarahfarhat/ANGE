#!/usr/bin/env nextflow

// id carried through untouched. Published (unlike before) so every source's raw, per-feature
// alignment evidence - BEFORE mergegff's monoexonic merge, and well before combineevidence's
// cross-source concat - is available as its own file: one gff3 per protein-database and
// transcriptome-assembly file, in absolute genome coordinates, ready to load as a track in a
// genome browser (IGV, JBrowse, ...) alongside results/gmove/gmove_out.gff3 to see exactly
// what each individual source of evidence looked like before any merging happened.
process exonerate2gff {
    publishDir { "${params.outdir}/exonerate" }, mode: 'copy'

    input:
    tuple val(id), path(exo_dir)

    output:
    tuple val(id), path("exonerate.${id}.gff3")

    script:
    """
    #Keep only real targetgff feature lines (9 tab fields, target id is our region format chr:start-end),
    #dropping exonerate's banner/header/footer text. NB: unlike the original doc script (which grepped on
    #hardcoded "^Scaffold|^contig" prefixes), this matches on the region-id format we control from
    #extractseq_v4.py, so it isn't tied to any particular genome's sequence naming.
    for f in ${exo_dir}/*.exo; do
        awk -F'\t' 'NF==9 && \$1 ~ /:[0-9]+-[0-9]+\$/ && (\$3=="gene" || \$3=="exon")' "\$f"
    done > exonerate_out.gff2

    #Remap region-local coordinates back to absolute genome coordinates and renumber gene/transcript IDs.
    #NB: coordinates are 1-based; local position 1 of the extracted region corresponds to the region's
    #start position, hence the "-1" (the original doc draft omitted this and was off by one base).
    awk -F"\t" 'BEGIN{OFS="\t"}
    {
        split(\$1,i,":|-");
        split(\$9,t," ");
        if(\$3=="gene"){
            cpt=cpt+1;
            \$9="ID=EXONAME_"sprintf("%08d",cpt)";"t[4]"="t[5]";"t[7]"="t[8]";"t[10]"="t[11];
        } else {
            \$9="Parent=EXOTRANS1_"sprintf("%08d",cpt)";"t[7]"="t[8]";"t[10]"="t[11];
        }
        print i[1],\$2,\$3,i[2]+\$4-1,i[2]+\$5-1,\$6,\$7,\$8,\$9
    }' exonerate_out.gff2 > exonerate.${id}.gff3
    """
}
