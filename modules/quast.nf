process quast {
    label 'quast'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs), path(illumina)
    output:
        path("quast")
    script:
        """
        quast.py -o quast -t ${task.cpus} --conserved-genes-finding \
          --gene-finding --pe1 ${illumina[0]} --pe2 ${illumina[1]} ${contigs}
        """
}
