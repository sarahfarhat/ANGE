#!/usr/bin/env nextflow

process buscodownload {
    container 'community.wave.seqera.io/library/busco:6.1.0--6d1f7006d91892b3'

    input:
    val lineage

    output:
    path "busco_downloads/lineages/${lineage}"

    script:
    // Downloads+decompresses the lineage dataset exactly once, shared by all four
    // buscoassess* calls below (each passed this directory via -l /path/to/dataset instead
    // of a bare lineage name) - without this, each of the 4 calls independently re-downloads
    // and re-decompresses its own full copy in its own isolated task work dir, which is both
    // slow and was the direct cause of a "No space left on device" failure on a real run.
    """
    busco --download ${lineage}
    """
}
