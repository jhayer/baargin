process compile_amrfinder {
    label 'python3'
    publishDir "${params.output}/compile_results", mode: 'copy'

    input:
        path(amr_files)
        path(amr_files_allmut)
        val(deconta)
    output:
        path("amrfinder_heatmap_${deconta}.tsv"), emit: amrf_hm
        path("amrfinder_allmut_heatmap_${deconta}.tsv"), emit: amrf_allmut_hm
    script:
        """
        mkdir amrfiles_${deconta}
        cp ${amr_files} amrfiles_${deconta}

        mkdir amrfiles_allmut_${deconta}
        cp ${amr_files_allmut} amrfiles_allmut_${deconta}

        amrfinder_compile_heatmap.py \
          -i amrfiles_${deconta} -o amrfinder_heatmap_${deconta}.tsv -s _AMRfinder.txt

        amrfinder_compile_heatmap.py \
          -i amrfiles_allmut_${deconta} -o amrfinder_allmut_heatmap_${deconta}.tsv -s _AMRfinder_all_mut.txt
        """
}

process compile_amrfinder_no_species {
    label 'python3'
    publishDir "${params.output}/compile_results", mode: 'copy'

    input:
        path(amr_files)
        val(deconta)
    output:
        path("amrfinder_heatmap_${deconta}.tsv"), emit: amrf_hm
    script:
        """
        mkdir amrfiles_${deconta}
        cp ${amr_files} amrfiles_${deconta}

        amrfinder_compile_heatmap.py \
          -i amrfiles_${deconta} -o amrfinder_heatmap_${deconta}.tsv -s _AMRfinder.txt
        """
}

process compile_amrfinder_plasmid_split {
    label 'python3'
    publishDir "${params.output}/compile_results", mode: 'copy'

    input:
        path(amr_files_plas)
        path(amr_files_allmut_plas)
        path(amr_files_chrom)
        path(amr_files_allmut_chrom)
        val(deconta)
    output:
        path("amrfinder_heatmap_${deconta}_plasmid.tsv"), emit: amrf_hm_p
        path("amrfinder_allmut_heatmap_${deconta}_plasmid.tsv"), emit: amrf_allmut_hm_p
        path("amrfinder_heatmap_${deconta}_chrom.tsv"), emit: amrf_hm_c
        path("amrfinder_allmut_heatmap_${deconta}_chrom.tsv"), emit: amrf_allmut_hm_c
    script:
        """
        mkdir amrfiles_${deconta}_plasmid
        cp ${amr_files_plas} amrfiles_${deconta}_plasmid

        mkdir amrfiles_allmut_${deconta}_plasmid
        cp ${amr_files_allmut_plas} amrfiles_allmut_${deconta}_plasmid

        amrfinder_compile_heatmap.py \
          -i amrfiles_${deconta}_plasmid -o amrfinder_heatmap_${deconta}_plasmid.tsv -s _AMRfinder_plasmid.txt

        amrfinder_compile_heatmap.py \
          -i amrfiles_allmut_${deconta}_plasmid -o amrfinder_allmut_heatmap_${deconta}_plasmid.tsv -s _AMRfinder_all_mut_plasmid.txt

        
        mkdir amrfiles_${deconta}_chrom
        cp ${amr_files_chrom} amrfiles_${deconta}_chrom

        mkdir amrfiles_allmut_${deconta}_chrom
        cp ${amr_files_allmut_chrom} amrfiles_allmut_${deconta}_chrom

        amrfinder_compile_heatmap.py \
          -i amrfiles_${deconta}_chrom -o amrfinder_heatmap_${deconta}_chrom.tsv -s _AMRfinder_chrom.txt

        amrfinder_compile_heatmap.py \
          -i amrfiles_allmut_${deconta}_chrom -o amrfinder_allmut_heatmap_${deconta}_chrom.tsv -s _AMRfinder_all_mut_chrom.txt
        """
}


process compile_amrfinder_no_sp_plasmid_split {
    label 'python3'
    publishDir "${params.output}/compile_results", mode: 'copy'

    input:
        path(amr_files_plas)
        path(amr_files_chrom)
        val(deconta)
    output:
        path("amrfinder_heatmap_${deconta}_plasmid.tsv"), emit: amrf_hm_p
        path("amrfinder_heatmap_${deconta}_chrom.tsv"), emit: amrf_hm_c
    script:
        """
        mkdir amrfiles_${deconta}_plasmid
        cp ${amr_files_plas} amrfiles_${deconta}_plasmid

        amrfinder_compile_heatmap.py \
          -i amrfiles_${deconta}_plasmid -o amrfinder_heatmap_${deconta}_plasmid.tsv -s _AMRfinder_plasmid.txt

        
        mkdir amrfiles_${deconta}_chrom
        cp ${amr_files_chrom} amrfiles_${deconta}_chrom

        amrfinder_compile_heatmap.py \
          -i amrfiles_${deconta}_chrom -o amrfinder_heatmap_${deconta}_chrom.tsv -s _AMRfinder_chrom.txt
        """
}
