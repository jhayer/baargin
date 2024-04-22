process platon {
    label 'platon'
    errorStrategy 'ignore'
    publishDir "${params.output}/${id}/plasmids", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        path(platon_db)
        val(deconta)
    output:
        path("platon_accu_${deconta}")
        tuple val(id), path("platon_accu_${deconta}/*.json"), emit: platon_json optional true
        tuple val(id), path("platon_accu_${deconta}/*.tsv"), emit: tp_platon_id_tsv optional true
    script:
        """
        platon --db ${platon_db} -o platon_accu_${deconta} ${contigs}
        """
}

process platon_json2tsv {
    label 'python3'
    publishDir "${params.output}/${id}/plasmids", mode: 'copy'

    input:
        tuple val(id), path(json_platon)
        val(deconta)
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

process compile_platon {
    label 'python3'
    publishDir "${params.output}/compile_results", mode: 'copy'

    input:
        path(inc_files)
        path(plamids_files)
        path(amr_files)
        val(deconta)
    output:
        path("platon_inc_${deconta}.tsv"), emit: inc_hm
        path("platon_plasmids_${deconta}.tsv"), emit: plasmids_hm
        path("platon_amr_${deconta}.tsv"), emit: amr_hm
    script:
        """
        mkdir platon_inc_${deconta}
        cp ${inc_files} platon_inc_${deconta}

        platon_compile.py \
          -i platon_inc_${deconta} -o platon_inc_${deconta}.tsv -t inc -s _platon_inctypes_${deconta}.tsv


        mkdir platon_plasmids_${deconta}
        cp ${plamids_files} platon_plasmids_${deconta}

        platon_compile.py \
          -i platon_plasmids_${deconta} -o platon_plasmids_${deconta}.tsv -t plasmid -s _platon_plasmids_${deconta}.tsv


        mkdir platon_amr_${deconta}
        cp ${amr_files} platon_amr_${deconta}

        platon_compile.py \
          -i platon_amr_${deconta} -o platon_amr_${deconta}.tsv -t amr -s _platon_amr_${deconta}.tsv
        """
}
