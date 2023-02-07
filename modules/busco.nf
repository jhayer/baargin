process busco {
    label 'busco'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        val(lineage)
        val(deconta)
        path(busco_dl_db)
    output:
        path("busco_${deconta}")
        path("${id}_busco_${deconta}.txt"), emit: busco_sum
    script:
        """
        if [ -z "${busco_dl_db}" ]
        then
          busco -i ${contigs} -o busco_${deconta} --mode genome -f --lineage_dataset ${lineage}
        else
          busco -i ${contigs} -o busco_${deconta} --mode genome \
            --offline --download_path ${busco_dl_db}  -f --lineage_dataset ${lineage}
        fi

        mv busco_${deconta}/short_summary.*.busco_${deconta}.txt ${id}_busco_${deconta}.txt
        """
}

process busco_auto_prok {
    label 'busco'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        val(deconta)
        path(busco_dl_db)
    output:
        path("busco_${deconta}")
        path "${id}_busco_${deconta}.txt", emit: busco_sum
    script:
        """
        if [ -z "${busco_dl_db}" ]
        then
          busco -i ${contigs} -o busco_${deconta} --mode genome -f --auto-lineage-prok
        else
          busco -i ${contigs} -o busco_${deconta} --mode genome \
            --offline --download_path ${busco_dl_db}  -f --auto-lineage-prok
        fi

        mv busco_${deconta}/short_summary*.busco_${deconta}.txt ${id}_busco_${deconta}.txt
        """
}

process busco_proteins {
    label 'busco'
    publishDir "${params.output}/${id}/annotation", mode: 'copy'
    input:
        tuple val(id), path(proteins)
        val(lineage)
        path(busco_dl_db)
    output:
        path("busco_proteins")
        path "${id}_busco_proteins.txt", emit: busco_prot
    script:
        """
        if [ -z "${busco_dl_db}" ]
        then
          busco -i ${proteins} -o busco_proteins --mode proteins -f --lineage_dataset ${lineage}
        else
          busco -i ${proteins} -o busco_proteins --mode proteins \
            --offline --download_path ${busco_dl_db}  -f --lineage_dataset ${lineage}
        fi

        mv busco_proteins/short_summary.*.busco_proteins.txt ${id}_busco_proteins.txt
        """
}

process busco_proteins_auto_prok {
  label 'busco'
  publishDir "${params.output}/${id}/annotation", mode: 'copy'

  input:
        tuple val(id), path(proteins)
        path(busco_dl_db)
    output:
        path("busco_proteins")
        path "${id}_busco_proteins.txt", emit: busco_prot
    script:
        """
        if [ -z "${busco_dl_db}" ]
        then
          busco -i ${proteins} -o busco_proteins --mode proteins -f --auto-lineage-prok
        else
          busco -i ${proteins} -o busco_proteins --mode proteins \
            --offline --download_path ${busco_dl_db}  -f --auto-lineage-prok
        fi

        mv busco_proteins/short_summary*.busco_proteins.txt ${id}_busco_proteins.txt
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
