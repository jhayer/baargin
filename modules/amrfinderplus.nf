process amrfinderplus {
    label 'amrfinderplus'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        val(species)
        path(local_amr_db)
        val(deconta)
    output:
        path("${id}_${deconta}_AMRfinder.txt"), emit: amrfile
        path("${id}_${deconta}_AMRfinder_all_mut.txt"), emit: amrfile_allmut
        tuple val(id), path("${id}_${deconta}_AMRfinder.txt"), path("${id}_${deconta}_AMRfinder_all_mut.txt"), emit: tp_id_amrf
    script:
        """
        amrfinder --database ${local_amr_db} --nucleotide ${contigs} -o ${id}_${deconta}_AMRfinder.txt --plus \
          --organism ${species} --mutation_all ${id}_${deconta}_AMRfinder_all_mut.txt
        """
}

process amrfinderplus_no_db {
    label 'amrfinderplus'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        val(species)
        val(deconta)
    output:
        path("${id}_${deconta}_AMRfinder.txt"), emit: amrfile
        path("${id}_${deconta}_AMRfinder_all_mut.txt"), emit: amrfile_allmut
        tuple val(id), path("${id}_${deconta}_AMRfinder.txt"), path("${id}_${deconta}_AMRfinder_all_mut.txt"), emit: tp_id_amrf
    script:
        """
        amrfinder --nucleotide ${contigs} -o ${id}_${deconta}_AMRfinder.txt --plus \
          --organism ${species} --mutation_all ${id}_${deconta}_AMRfinder_all_mut.txt
        """
}

process amrfinderplus_no_species {
    label 'amrfinderplus'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        path(local_amr_db)
        val(deconta)
    output:
        path("${id}_${deconta}_AMRfinder.txt"), emit: amrfile
        tuple val(id), path("${id}_${deconta}_AMRfinder.txt"), emit: tp_id_amrf
    script:
        """
        amrfinder --database ${local_amr_db} --nucleotide ${contigs} -o ${id}_${deconta}_AMRfinder.txt --plus
        """
}

process amrfinderplus_no_species_no_db {
    label 'amrfinderplus'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        val(deconta)
    output:
        path("${id}_${deconta}_AMRfinder.txt"), emit: amrfile
        tuple val(id), path("${id}_${deconta}_AMRfinder.txt"), emit: tp_id_amrf
    script:
        """
        amrfinder --nucleotide ${contigs} -o ${id}_${deconta}_AMRfinder.txt --plus
        """
}
