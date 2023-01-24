process quast {
    label 'quast'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(illumina), path(contigs)
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
process quast_contigs_only {
    label 'quast_contigs_only'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        val(deconta)
    output:
        path("quast_${deconta}")
        path("${id}_quast_${deconta}_report.tsv"), emit: quast_report
    script:
        """
        quast.py -o quast_${deconta} -t ${task.cpus} --no-plots --no-icarus \
            ${contigs}

        mv quast_${deconta}/report.tsv ${id}_quast_${deconta}_report.tsv
        """
}
process quast_hybrid {
    label 'quast_hybrid'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        tuple val(id), path(illuminaR1), path(illuminaR2), path(ont)
        val(deconta)
    output:
        path("quast_${deconta}")
        path("${id}_quast_${deconta}_report.tsv"), emit: quast_report
    script:
        """
        quast.py -o quast_${deconta} -t ${task.cpus} --no-plots --no-icarus \
            --pe1 ${illuminaR1} --pe2 ${illuminaR2} --nanopore ${ont} ${contigs}

        mv quast_${deconta}/report.tsv ${id}_quast_${deconta}_report.tsv
        """
}
