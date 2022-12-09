process platon {
    label 'platon'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        path(platon_db)
        val(deconta)
    output:
        path("platon_accu_${deconta}")
        path("${id}_platon_${deconta}_results.tsv"), emit: platon_tab
    script:
        """
        platon --db ${platon_db} -o platon_accu_${deconta} ${contigs}

        mv platon_accu_${deconta}/*.tsv ${id}_platon_${deconta}_results.tsv
        """
}
