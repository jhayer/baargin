process busco {
    label 'busco'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        val(lineage)
    output:
        path("busco")
    script:
        """
        busco -i ${contigs} -o busco --mode genome --lineage_dataset ${lineage}
        """
}

process busco_auto_prok {
    label 'busco_auto_prok'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
    output:
        tuple val(id), path("busco/..")
        path("busco")
    script:
        """
        busco -i ${contigs} -o busco --mode genome --auto-lineage-prok
        """
}
