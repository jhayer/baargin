process compile_plasmidfinder {
    label 'python3.9'
    publishDir "${params.output}/compile_results", mode: 'copy'

    input:
        path(pf_files)
        val(deconta)
    output:
        path("plasmidfinder_heatmap_${deconta}.tsv"), emit: pf_hm
    script:
        """
        mkdir plasmidfinder_${deconta}
        cp ${pf_files} plasmidfinder_${deconta}

        plasmidfinder_compile.py \
          -i plasmidfinder_${deconta} -o plasmidfinder_heatmap_${deconta}.tsv -s _results.tsv
        """
}
