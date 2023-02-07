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
        path("platon_accu_${deconta}/*.json"), emit: platon_json
        val(id), emit: platon_id
    script:
        """
        platon --db ${platon_db} -o platon_accu_${deconta} ${contigs}

        mv platon_accu_${deconta}/*.tsv ${id}_platon_${deconta}_results.tsv
        """
}

process platon_json2tsv {
    label 'python3'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'

    input:
        path(json_platon)
        val(deconta)
        val(id)
    output:
        path("${id}_platon_inctypes_${deconta}.tsv"), emit: platon_inc
        path("${id}_platon_plasmids_${deconta}.tsv"), emit: platon_plasmid
        path("${id}_platon_amr_${deconta}.tsv"), emit: platon_amr

    script:
        """
        json2tsv_platon.py \
          -j ${json_platon} -i ${id}_platon_inctypes_${deconta}.tsv -p ${id}_platon_plasmids_${deconta}.tsv -a ${id}_platon_amr_${deconta}.tsv
        
        """

}
