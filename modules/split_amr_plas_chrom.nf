process split_amr_plas_chrom_amrfinder_sp {
    label 'python3'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'

    input:
        tuple val(id), path(amr_file), path(amr_file_allmut), path(platon_sum_tsv)
        val(deconta)
    output:
        path("${id}_${deconta}_AMRfinder_plasmid.txt"), emit: amrf_plasmid
        path("${id}_${deconta}_AMRfinder_all_mut_plasmid.txt"), emit: amrf_allmut_plasmid
        path("${id}_${deconta}_AMRfinder_chrom.txt"), emit: amrf_chrom
        path("${id}_${deconta}_AMRfinder_all_mut_chrom.txt"), emit: amrf_allmut_chrom
    script:
        """
        split_AMR_out_by_plasmid.py -i ${platon_sum_tsv} -a ${amr_file} -t amrfinder -p ${id}_${deconta}_AMRfinder_plasmid.txt -c ${id}_${deconta}_AMRfinder_chrom.txt

        split_AMR_out_by_plasmid.py -i ${platon_sum_tsv} -a ${amr_file_allmut} -t amrfinder -p ${id}_${deconta}_AMRfinder_all_mut_plasmid.txt -c ${id}_${deconta}_AMRfinder_all_mut_chrom.txt
        """
}

process split_amr_plas_chrom_amrfinder_no_sp {
    label 'python3'
    publishDir "${params.output}/${id}/AMR/amrfinderplus", mode: 'copy'

    input:
        tuple val(id), path(amr_file), path(platon_sum_tsv)
        val(deconta)
    output:
        path("${id}_${deconta}_AMRfinder_plasmid.txt"), emit: amrf_plasmid
        path("${id}_${deconta}_AMRfinder_chrom.txt"), emit: amrf_chrom
    script:
        """
        split_AMR_out_by_plasmid.py -i ${platon_sum_tsv} -a ${amr_file} -t amrfinder -p ${id}_${deconta}_AMRfinder_plasmid.txt -c ${id}_${deconta}_AMRfinder_chrom.txt
        """
}
