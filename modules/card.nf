process card_rgi {
    label 'card_rgi'
    publishDir "${params.output}/${id}/AMR/card_rgi", mode: 'copy'
    conda 'rgi'

    input:
        tuple val(id), path(contigs)
        path(card_json_db)
    output:
        path("${id}_RGI_main.json")
        path("${id}_RGI_main.txt")
    script:
        """
        rgi load --card_json ${card_json_db}

        rgi main -i ${contigs} -o ${id}_RGI_main --debug -a BLAST -d wgs -n ${task.cpus}
        rgi tab -i ${id}_RGI_main.json

        rm ${contigs}.temp.*
        rm ${contigs}.temp
        """
}
