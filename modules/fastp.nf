process fastp {
    label 'fastp'
    publishDir "${params.output}/${id}/qc", mode: 'copy'
    input:
        tuple val(id), path(illumina)
    output:
        tuple val(id), path("*_R?_clean.fastq")
        path("${id}_fastp_report.html")
    script:
        """
        fastp -i ${illumina[0]} -I ${illumina[1]} \
            -o ${id}_R1_clean.fastq -O ${id}_R2_clean.fastq \
            --detect_adapter_for_pe --html ${id}_fastp_report.html
        """
}
