process plasmidfinder {
    label 'plasmidfinder'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        path(pf_db)
    output:
        path("plasmidfinder")
    script:
        """
        mkdir plasmidfinder
        plasmidfinder.py -i ${contigs} -o plasmidfinder -p ${pf_db} -mp blastn -x

        plasmidfinder.py -i ${contigs} -o plasmidfinder -p ${pf_db} -mp kma -x
        """
}
