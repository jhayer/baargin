process plasmidfinder {
    label 'plasmidfinder'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        path(pf_db)
        val(deconta)
    output:
        path("plasmidfinder_${deconta}")
        path("${id}_plasmidfinder_${deconta}_results.tsv"), emit: plasmidfinder_tab
    script:
        """
        mkdir plasmidfinder_${deconta}
        plasmidfinder.py -i ${contigs} -o plasmidfinder_${deconta} -p ${pf_db} -mp blastn -x

        plasmidfinder.py -i ${contigs} -o plasmidfinder_${deconta} -p ${pf_db} -mp kma -x

        mv plasmidfinder_${deconta}/results_tab.tsv ${id}_plasmidfinder_${deconta}_results.tsv
        """
}
