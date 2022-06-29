process card_rgi {
    label 'card_rgi'
    publishDir "${params.output}/${id}/AMR/card_rgi", mode: 'copy'
    conda '/home/hayer/.conda/envs/rgi'
//    conda "${params.conda_card_rgi}"

    input:
        tuple val(id), path(contigs)
        path(card_json_db)
        val(deconta)
    output:
        path("${id}_${deconta}_RGI_main.json")
        path("${id}_${deconta}_RGI_main.txt")
    script:
        """
        rgi load --card_json ${card_json_db}

        rgi main -i ${contigs} -o ${id}_${deconta}_RGI_main --debug -a BLAST -d wgs -n ${task.cpus}
        rgi tab -i ${id}_${deconta}_RGI_main.json

        rm ${contigs}.temp.*
        rm ${contigs}.temp
        """
}
