process kraken2nt_contigs {
    label 'kraken2nt_contigs'
    publishDir "${params.output}/${id}/taxonomic_classif/contigs", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        path(db_k2nt)
    output:
        path("*_kn2_nt-re*.txt")
    script:
        """
        kraken2 --db ${db_k2nt} --memory-mapping \
            --threads ${task.cpus} --output ${id}_kn2_nt-res.txt \
            --report-minimizer-data \
            --report ${id}_kn2_nt-report.txt ${contigs}
        """
}

//add here KrakenTools - retrieval of species of interest
