process busco {
    label 'busco'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        val(lineage)
        val(deconta)
        path(busco_dl_db)
    output:
        path("busco_${deconta}")
        path "${id}_${deconta}_summary_enterobacterales_odb10_busco.txt", emit: busco_sum
    script:
        """
        if [ -z "${busco_dl_db}" ]
        then
          busco -i ${contigs} -o busco_${deconta} --mode genome -f --lineage_dataset ${lineage}
        else
          busco -i ${contigs} -o busco_${deconta} --mode genome \
            --offline --download_path ${busco_dl_db}  -f --lineage_dataset ${lineage}
        fi
        
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
        busco -i ${contigs} -o busco_${deconta} --mode genome -f --auto-lineage-prok

        mv busco_${deconta}/short_summary*.busco_${deconta}.txt ${id}_${deconta}_summary_autoprok_busco.txt
        """
}
