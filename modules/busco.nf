process busco {
    label 'busco'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs), lineage
    output:
        tuple val(id), path("busco/..")
        path("busco")
    script:
        """
        busco -i ${contigs} -o busco --mode genome --lineage_dataset ${lineage}
        """
}
