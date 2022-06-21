process prokka {
    label 'prokka'
    publishDir "${params.output}/${id}/annotation", mode: 'copy'


    input:
        tuple val(id), path(contigs)
        val(genus)
        val(species)
    output:
        path("${id}_prokka")
        path("${id}_prokka.gff"), emit: prokka_gff
    script:
        """
        prokka --force --genus ${genus} --species ${species} \
          --kingdom Bacteria --usegenus \
          --notrna --prefix ${id} --outdir ${id}_prokka ${contigs}

        mv ${id}_prokka/${id}.gff ${id}_prokka.gff
        """
}
