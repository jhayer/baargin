process compile_card {
    label 'compile_card'
    publishDir "${params.output}/compile_results", mode: 'copy'

    input:
        path(json_files)
        val(deconta)
    output:
        path("hm_drugclass_samples_${deconta}*.csv"), emit: card_hm_dc
        path("hm_genefamily_samples_${deconta}*.csv"), emit: card_hm_gene
        path("hm_*.png")
    script:
        """
        mkdir card_${deconta}
        cp ${json_files} card_${deconta}

        # Heatmap of samples clustered by similarity of resistome and AMR genes organized by drug class
        rgi heatmap -i card_${deconta} -cat drug_class -o hm_drugclass_samples_${deconta} -clus samples

        # Heatmap of samples clustered by similarity of resistome and AMR genes organized by AMR gene families
        rgi heatmap -i card_${deconta} -cat gene_family -o hm_genefamily_samples_${deconta} -clus samples
        """
}
