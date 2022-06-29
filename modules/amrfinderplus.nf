process amrfinderplus {
    label 'amrfinderplus'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'
  //  scratch './tmpdir/'

    input:
        tuple val(id), path(contigs)
        val(species)
        val(deconta)
    output:
        path("${id}_${deconta}_AMRfinder.txt")
        path("${id}_${deconta}_AMRfinder_all_mut.txt")
    script:
        """
        amrfinder --nucleotide ${contigs} -o ${id}_${deconta}_AMRfinder.txt --plus \
            --organism ${species} --mutation_all ${id}_${deconta}_AMRfinder_all_mut.txt
        """
}

process amrfinderplus_no_species {
    label 'amrfinderplus_no_species'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        val(deconta)
    output:
        path("${id}_${deconta}_AMRfinder.txt")
    script:
        """
        amrfinder --nucleotide ${contigs} -o ${id}_${deconta}_AMRfinder.txt --plus
        """
}
