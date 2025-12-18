process bakta {
    label 'bakta'
    publishDir "${params.output}/${id}/annotation", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        path(bakta_db)
        val(genus)
        val(species)
    output:
        path("${id}_bakta")
        path("${id}_bakta.gff"), emit: annot_gff
        tuple val(id), path("${id}_bakta/${id}.faa"), emit: annot_faa
    script:
        """
        bakta --db ${bakta_db} --keep-contig-headers --genus ${genus} \
          --species ${species} --skip-trna --prefix ${id} --output ${id}_bakta \
          --verbose --skip-rrna ${contigs} --skip-plot --force

        mv ${id}_bakta/${id}.gff3 ${id}_bakta.gff
        """
}

process bakta_genus {
    label 'bakta'
    publishDir "${params.output}/${id}/annotation", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        path(bakta_db)
        val(genus)
    output:
        path("${id}_bakta")
        path("${id}_bakta.gff"), emit: annot_gff
        tuple val(id), path("${id}_bakta/${id}.faa"), emit: annot_faa
    script:
        """
        bakta --db ${bakta_db} --keep-contig-headers --genus ${genus} \
          --skip-trna --prefix ${id} --output ${id}_bakta \
          --verbose --skip-rrna ${contigs} --skip-plot --force

        mv ${id}_bakta/${id}.gff3 ${id}_bakta.gff
        """
}

process bakta_plasmids {
    label 'bakta'
    publishDir "${params.output}/${id}/annotation", mode: 'copy'

    input:
        tuple val(id), path(plasmid_contigs)
        path(bakta_db)
        val(genus)
    output:
        path("${id}_bakta_plasmids")
        path("${id}_bakta_plasmids.gff"), emit: annot_gff
        tuple val(id), path("${id}_bakta_plasmids/${id}.faa"), emit: annot_faa
    script:
        """
        bakta --db ${bakta_db} --keep-contig-headers --plasmid unnamed \
          --skip-trna --prefix ${id} --output ${id}_bakta_plasmids \
          --verbose --skip-rrna --skip-plot --force ${plasmid_contigs}

        mv ${id}_bakta_plasmids/${id}.gff3 ${id}_bakta_plasmids.gff
        """
}
