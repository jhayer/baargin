process compile_card_split_plasmid {
    label 'python3'
    publishDir "${params.output}/compile_results", mode: 'copy'

    input:
        path(amr_files_plas)
        path(amr_files_chrom)
        val(deconta)
    output:
        path("card_heatmap_${deconta}_plasmid.tsv"), emit: card_hm_p
        path("card_heatmap_${deconta}_chrom.tsv"), emit: card_hm_c
    script:
        """
        mkdir cardfiles_${deconta}_plasmid
        cp ${amr_files_plas} cardfiles_${deconta}_plasmid

        card_compile.py -i cardfiles_${deconta}_plasmid -o card_heatmap_${deconta}_plasmid.tsv -s _RGI_main_plasmid.txt

        
        mkdir cardfiles_${deconta}_chrom
        cp ${amr_files_chrom} cardfiles_${deconta}_chrom

        card_compile.py -i cardfiles_${deconta}_chrom -o card_heatmap_${deconta}_chrom.tsv -s _RGI_main_chrom.txt
        """
}