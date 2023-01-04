process roary {
    label 'roary'
    publishDir "${params.output}/pangenome", mode: 'copy'

    input:
        path(prokka_gff)

    output:
        path("roary_out"), emit: roary
    script:
        """
        roary -p ${task.cpus} -o roary_clust -f roary_out -e -n -r -v  ${prokka_gff}
        """
}
