process compile_amrfinder {
    label 'python3.9'
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
    label 'compile_amrfinder_no_species'
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
