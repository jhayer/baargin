process quast {
    label 'quast'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(illumina), path(contigs)
        val(deconta)
    output:
        path("${id}_quast_${deconta}_transposed_report.tsv"), emit: quast_transpo
        path("${id}_quast_${deconta}_report.tsv"), emit: quast_report
    script:
        """
        quast.py -o quast_${deconta} -t ${task.cpus} --no-plots --no-icarus \
            --pe1 ${illumina[0]} --pe2 ${illumina[1]} ${contigs}

        mv quast_${deconta}/report.tsv ${id}_quast_${deconta}_report.tsv
        mv quast_${deconta}/transposed_report.tsv ${id}_quast_${deconta}_transposed_report.tsv

        rm -r quast_${deconta}
        rm ${illumina[0]}
        rm ${illumina[1]}
        """
}
process quast_contigs_only {
    label 'quast'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        val(deconta)
    output:
        path("${id}_quast_${deconta}_transposed_report.tsv"), emit: quast_transpo
        path("${id}_quast_${deconta}_report.tsv"), emit: quast_report
    script:
        """
        quast.py -o quast_${deconta} -t ${task.cpus} --no-plots --no-icarus \
            ${contigs}

        mv quast_${deconta}/report.tsv ${id}_quast_${deconta}_report.tsv
        mv quast_${deconta}/transposed_report.tsv ${id}_quast_${deconta}_transposed_report.tsv

        rm -r quast_${deconta}
        """
}
process quast_hybrid {
    label 'quast'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        tuple val(id), path(illuminaR1), path(illuminaR2), path(ont)
        val(deconta)
    output:
        path("${id}_quast_${deconta}_transposed_report.tsv"), emit: quast_transpo
        path("${id}_quast_${deconta}_report.tsv"), emit: quast_report

    script:
        """
        quast.py -o quast_${deconta} -t ${task.cpus} --no-plots --no-icarus \
            --pe1 ${illuminaR1} --pe2 ${illuminaR2} --nanopore ${ont} ${contigs}

        mv quast_${deconta}/report.tsv ${id}_quast_${deconta}_report.tsv
        mv quast_${deconta}/transposed_report.tsv ${id}_quast_${deconta}_transposed_report.tsv

        rm -r quast_${deconta}
        rm ${illuminaR1}
        rm ${illuminaR2}        
        """
}

process compile_quast {
    publishDir "${params.output}/compile_results", mode: 'copy'

    input:
        path(quast_files)
        val(deconta)
    output:
        path("quast_${deconta}.tsv"), emit: quast_compile
    script:
        """bash 
        mkdir quast_out_${deconta}
        cp ${quast_files} quast_out_${deconta}

        echo "Assembly #contigs(>=0bp) #contigs(>=10000bp) #contigs(>=50000bp) #contigs(>=500bp) Largest_contig Total_length GC(%) N50 N75 L50 L75 Avg.coverage_depth" > quast_${deconta}.tsv 
       
        for i in quast_out_${deconta}/*.tsv
        do 
            cat \${i} | awk -F'\t' '{ print \$1,\$2,\$5,\$7,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21,\$27 }' | grep -v 'Assembly' >> quast_${deconta}.tsv 
        done

        """
}
