process busco {
    label 'busco'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        val(deconta)
        val(command)
    output:
        path("busco_${deconta}")
        path("${id}_busco_${deconta}.txt"), emit: busco_sum
    script:
        """
        busco -i ${contigs} -o busco_${deconta} --mode genome -f ${command}

        mv busco_${deconta}/short_summary.*.busco_${deconta}.txt ${id}_busco_${deconta}.txt
        """
}

process busco_proteins {
    label 'busco'
    publishDir "${params.output}/${id}/annotation", mode: 'copy'
    input:
        tuple val(id), path(proteins)
        val(prot_command)
    output:
        path("busco_proteins")
        path "${id}_busco_proteins.txt", emit: busco_prot
    script:
        """
        busco -i ${proteins} -o busco_proteins --mode proteins -f ${prot_command}

        mv busco_proteins/short_summary.*.busco_proteins.txt ${id}_busco_proteins.txt
        """
}

process compile_busco {

  publishDir "${params.output}/compile_results", mode: 'copy'

  input:
        path(busco_files)
        val(deconta)
    output:
        path("busco_assembly_${deconta}.tsv"), emit: busco_compile
    script:
        """bash 
        mkdir busco_files_${deconta}
        cp ${busco_files} busco_files_${deconta}

        for i in busco_files_${deconta}/*.txt
        do 
          echo \$(basename \${i}) \$(grep 'C:' \${i}) >> busco_assembly_${deconta}_preformat.tsv
        done

        cat busco_assembly_${deconta}_preformat.tsv | awk -F':' '{ print \$1,\$2,\$3,\$4,\$5,\$6,\$7}'| awk -F'[' '{ print \$1,\$2}' | awk -F']' '{ print \$1,\$2}' > busco_assembly_${deconta}.tsv
        """
}

process compile_busco_prot {

  publishDir "${params.output}/compile_results", mode: 'copy'

  input:
        path(busco_files)
    output:
        path("busco_annotation.tsv"), emit: busco_prot_compile
    script:
        """bash 
        mkdir busco_annot
        cp ${busco_files} busco_annot

        for i in busco_annot/*.txt
        do 
          echo \$(basename \${i}) \$(grep 'C:' \${i}) >> busco_annot_preformat.tsv
        done

        cat busco_annot_preformat.tsv | awk -F':' '{ print \$1,\$2,\$3,\$4,\$5,\$6,\$7}'| awk -F'[' '{ print \$1,\$2}' | awk -F']' '{ print \$1,\$2}' > busco_annotation.tsv
        """
}
