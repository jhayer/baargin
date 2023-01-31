process prokka {
    label 'prokka'
    publishDir "${params.output}/${id}/annotation", mode: 'copy'


    input:
        tuple val(id), path(contigs)
        val(genus)
        val(species)
    output:
        path("${id}_prokka")
        path("${id}_prokka.gff"), emit: annot_gff
        path("${id}_prokka/${id}_scaffolds.faa"), emit: annot_faa
    script:
        """
        prokka --force --genus ${genus} --species ${species} \
          --kingdom Bacteria --usegenus \
          --notrna --prefix ${id} --outdir ${id}_prokka ${contigs}

        mv ${id}_prokka/${id}.gff ${id}_prokka.gff
        """
}
