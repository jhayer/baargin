process bakta {
    label 'bakta'
    publishDir "${params.output}/${id}/annotation", mode: 'copy'
  //  conda '/home/hayer/.conda/envs/bakta'

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
          --verbose --skip-rrna ${contigs} --skip-plot

        mv ${id}_bakta/${id}.gff3 ${id}_bakta.gff
        """
}
