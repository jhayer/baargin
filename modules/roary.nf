process roary {
    label 'roary'
    publishDir "${params.output}/pangenome", mode: 'copy'

    input:
        path(prokka_gff)

    output:
        path("roary_out"), emit: roary
        path("roary_out/core_gene_alignment.aln"), emit: core_aln
        path("roary_out/gene_presence_absence.csv"), emit: gene_matrix
    script:
        """
        roary -p ${task.cpus} -o roary_clust -f roary_out -e -n -r -v  ${prokka_gff}
        """
}

process roary_fasttree {
    label 'fasttree'
    publishDir "${params.output}/pangenome", mode: 'copy'
    input:
        path(core_gene_aln)
    output: 
        path("core_gene_fasttree.newick"), emit: fasttree_nwck
    script:
    """
    FastTree -nt -gtr ${core_gene_aln} > core_gene_fasttree.newick
    """
}

process roary_plots {
    label 'roary_plots'
    publishDir "${params.output}/pangenome", mode: 'copy'
    input:
        path(fasttree_nwck)
        path(roary_gene_matrix)
    output: 
        path("pangenome_*.png")
    script:
    """
    roary_plots.py --labels ${fasttree_nwck} ${roary_gene_matrix} 
    """
}
