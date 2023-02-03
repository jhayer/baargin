process mlst {
    label 'mlst'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        val(deconta)
    output:
        path("${id}_mlst_${deconta}.tsv"), emit: mlst_tab
    script:
        """
        mlst ${contigs} > ${id}_mlst_${deconta}.tsv
        """
}

process compile_mlst {

  publishDir "${params.output}/compile_results", mode: 'copy'

  input:
        path(mlst_files)
        val(deconta)
    output:
        path("mlst_compile_${deconta}.tsv"), emit: mlst_compile
    script:
        """bash 
        cat ${mlst_files} > mlst_compile_${deconta}.tsv
        """
}

