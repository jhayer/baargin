process spades {
    label 'spades'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(illumina)
    output:
        tuple val(id), path("spades/scaffolds.fasta")
        path("spades")
    script:
        """
        spades.py -1 ${illumina[0]} -2 ${illumina[1]} \
         -o spades -t ${task.cpus}
        """
}
