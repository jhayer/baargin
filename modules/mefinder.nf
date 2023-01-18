process mefinder {
    label 'mefinder'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        val(deconta)
    output:
        path("${id}_mefinder_${deconta}.csv"), emit: mefinder_tab
    script:
        """
         mefinder index
         mefinder find -c ${contigs} -g --db-path /tmp/mge_finder/database ${id}_mefinder_${deconta}
        """
}
