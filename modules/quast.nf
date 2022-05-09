process quast {
    label 'quast'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        tuple val(id), path(illumina)
        val(deconta)
    output:
        path("quast_${deconta}")
        path("${id}_quast_${deconta}_report.tsv"), emit: quast_report
    script:
        """
        quast.py -o quast_${deconta} -t ${task.cpus} --no-plots --no-icarus \
            --pe1 ${illumina[0]} --pe2 ${illumina[1]} ${contigs}

        mv quast_${deconta}/report.tsv ${id}_quast_${deconta}_report.tsv
        """
}
//quast.py -o quast -t ${task.cpus} --conserved-genes-finding \
  //  --gene-finding --pe1 ${illumina[0]} --pe2 ${illumina[1]} ${contigs}
