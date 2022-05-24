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
