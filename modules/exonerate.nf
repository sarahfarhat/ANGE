#!/usr/bin/env nextflow

process exonerate {
    conda 'bioconda::exonerate=2.4.0'
    input:
    path query_dir
    path genome_dir
    val mode
    val percentThreshold

    output:
    path "exonerate_out", emit: exo_dir

    script:
    """
    mkdir -p exonerate_out
    ls ${query_dir} | xargs -I{} -P ${task.cpus} bash -c '
        base=\$(basename "{}" .faa)
        exonerate -m ${mode} \\
            -q ${query_dir}/{} \\
            -t ${genome_dir}/\${base}.fa \\
            -Q \$([ "${mode}" = "est2genome" ] && echo dna || echo protein) \\
            -T dna -c 1 --percent ${percentThreshold} \\
            --showvulgar no --showalignment no --showtargetgff yes \\
            > exonerate_out/\${base}.exo
    '
    """
}
