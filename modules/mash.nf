process mash_screen {
    label 'mash_screen'
    publishDir "${params.output}/${id}/taxonomic_classif/contigs", mode: 'copy'

    input:
        tuple val(id), path(contigs)
        tuple(species)
        path(mash_sketch)
    output:
        path("${id}_${species}_screen_sort_gr.tab")
        path("${id}_3_mash_top_hits")
    script:
        """
        mash screen -p 4 ${mash_sketch} ${contigs} > ${id}_${species}_screen.tab
        sort -gr ${id}_${species}_screen.tab > ${id}_${species}_screen_sort_gr.tab

        # getting the 3 top hits of closest relative that are complete genomes
        grep -i 'complete genome' ${id}_${species}_screen_sort_gr.tab | head -n10 > ${id}_10_mash_top_hits

        # maybe later add retrieval of top AC from NCBI
        # file looks like this (col 5, row 1 to use as AC number):
        # 0.999761	995/1000	1	0	CP020354.1	Staphylococcus aureus strain FORC59 chromosome, complete genome
        """
}
