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
    script:
        """
        if [ -z "${local_amr_db}" ]
        then
          amrfinder --nucleotide ${contigs} -o ${id}_${deconta}_AMRfinder.txt --plus \
            --organism ${species} --mutation_all ${id}_${deconta}_AMRfinder_all_mut.txt
        else
          amrfinder --database ${local_amr_db} --nucleotide ${contigs} -o ${id}_${deconta}_AMRfinder.txt --plus \
            --organism ${species} --mutation_all ${id}_${deconta}_AMRfinder_all_mut.txt

        fi
        """
}

process amrfinderplus_no_species {
    label 'amrfinderplus_no_species'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        path(local_amr_db)
        val(deconta)
    output:
        path("${id}_${deconta}_AMRfinder.txt")
    script:
        """
        if [ -z "${local_amr_db}" ]
        then
          amrfinder --nucleotide ${contigs} -o ${id}_${deconta}_AMRfinder.txt --plus
        else
          amrfinder --database ${local_amr_db} --nucleotide ${contigs} -o ${id}_${deconta}_AMRfinder.txt --plus
        fi
        """
}
