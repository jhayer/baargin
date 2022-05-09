process busco {
    label 'busco'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        val(lineage)
        val(deconta)
    output:
        path("busco_${deconta}")
        path "${id}_${deconta}_summary_enterobacterales_odb10_busco.txt", emit: busco_sum
    script:
        """
        busco -i ${contigs} -o busco_${deconta} --mode genome --lineage_dataset ${lineage}

        mv busco_${deconta}/short_summary.specific.enterobacterales_odb10.busco_${deconta}.txt ${id}_${deconta}_summary_enterobacterales_odb10_busco.txt
        """
}

process busco_auto_prok {
    label 'busco_auto_prok'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        val(deconta)
    output:
        path("busco_${deconta}")
        path "${id}_${deconta}_summary_autoprok_busco.txt", emit: busco_sum
    script:
        """
        busco -i ${contigs} -o busco_${deconta} --mode genome --auto-lineage-prok

        mv busco_${deconta}/short_summary*.busco_${deconta}.txt ${id}_${deconta}_summary_autoprok_busco.txt
        """
}
