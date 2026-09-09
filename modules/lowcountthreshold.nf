#!/usr/bin/env nextflow

process lowcountthreshold {
    input:
    path counts
    val threshold

    output:
    path "to_remove.list"

    script:
    // NB: the doc's draft awk used `$7>10` under the step title "list low count genes", which
    // would flag well-supported genes for removal instead of poorly-supported ones. This keeps
    // the step's stated intent - genes below the read-count threshold go on the removal list.
    // featureCounts output: line 1 = program comment, line 2 = column header, so skip both.
    """
    awk -F"\t" -v t='${threshold}' 'NR>2 && \$7<t {print \$1}' ${counts} > to_remove.list
    """
}
