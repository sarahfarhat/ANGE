#!/usr/bin/env nextflow

process buscodownload {
    // Switched from the Seqera "Community Containers" auto-built busco image to the official
    // ezlabgva/busco image from the BUSCO team's own Docker Hub (verified current tag as of
    // this edit - https://hub.docker.com/r/ezlabgva/busco/tags). The community-built image hit
    // "busco: command not found" (exit 127) even though singularity/apptainer clearly launched
    // (its own INFO lines about SINGULARITYENV_* vars appeared in the log) - the most likely
    // explanation is that busco's bioconda recipe pulls in several dependencies (augustus,
    // metaeuk, hmmer, sepp...) whose PATH/env setup normally happens via a conda-activate hook,
    // and an auto-built mulled container doesn't always bake that into static image ENV the way
    // a purpose-built image does - `singularity exec` doesn't source that hook, so busco wasn't
    // on PATH. The official ezlabgva image is built and maintained by the BUSCO developers
    // specifically for container use (this exact image is what nf-core's busco module uses) and
    // doesn't have that failure mode. Every other container in this pipeline (seqkit, blast,
    // etc.) is a single, dependency-light static binary and isn't affected by this - no need to
    // change those.
    container 'ezlabgva/busco:v6.1.0_cv2'

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
