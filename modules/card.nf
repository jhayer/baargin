process card_rgi {
    label 'card'
    publishDir "${params.output}/${id}/AMR/card_rgi", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        path(card_json_db)
        val(deconta)
    output:
        path("${id}_${deconta}_RGI_main.json"), emit: card_json
        path("${id}_${deconta}_RGI_main.txt")
    script:
        """
        rgi load --card_json ${card_json_db} --local

        rgi main -i ${contigs} -o ${id}_${deconta}_RGI_main --debug \
          --local -a BLAST -d wgs -n ${task.cpus}
        rgi tab -i ${id}_${deconta}_RGI_main.json

        rm ${contigs}.temp.*
        rm ${contigs}.temp
        """
}

process compile_card {
    label 'card'
    publishDir "${params.output}/compile_results", mode: 'copy'

    input:
        path(json_files)
        val(deconta)
    output:
        path("hm_drugclass_samples_${deconta}*.csv"), emit: card_hm_dc
        path("hm_genefamily_samples_${deconta}*.csv"), emit: card_hm_gene
        path("hm_*.png")
    script:
        """
        mkdir card_${deconta}
        cp ${json_files} card_${deconta}

        # Heatmap of samples clustered by similarity of resistome and AMR genes organized by drug class
        rgi heatmap -i card_${deconta} -cat drug_class -o hm_drugclass_samples_${deconta} -clus samples

        # Heatmap of samples clustered by similarity of resistome and AMR genes organized by AMR gene families
        rgi heatmap -i card_${deconta} -cat gene_family -o hm_genefamily_samples_${deconta} -clus samples
        """
}
