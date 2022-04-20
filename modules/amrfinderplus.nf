process amrfinderplus {
    label 'amrfinderplus'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'
    conda "${params.conda_amrfinder}"

    input:
        tuple val(id), path(contigs)
        val(species)
    output:
        path("${id}_AMRfinder.txt")
        path("${id}_AMRfinder_all_mut.txt")
    script:
        """
        amrfinder --nucleotide ${contigs} -o ${id}_AMRfinder.txt --plus \
            --organism ${species} --mutation_all ${id}_AMRfinder_all_mut.txt
        """
}

process amrfinderplus_no_species {
    label 'amrfinderplus_no_species'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'
    conda 'amrfinderplus'

    input:
        tuple val(id), path(contigs)
    output:
        path("${id}_AMRfinder.txt")
    script:
        """
        amrfinder --nucleotide ${contigs} -o ${id}_AMRfinder.txt --plus
        """
}
