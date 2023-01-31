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
        path "${id}_${deconta}_busco.txt", emit: busco_sum
    script:
        """
        if [ -z "${busco_dl_db}" ]
        then
          busco -i ${contigs} -o busco_${deconta} --mode genome -f --lineage_dataset ${lineage}
        else
          busco -i ${contigs} -o busco_${deconta} --mode genome \
            --offline --download_path ${busco_dl_db}  -f --lineage_dataset ${lineage}
        fi

        mv busco_${deconta}/short_summary.*.busco_${deconta}.txt ${id}_${deconta}_busco.txt
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
        path "${id}_${deconta}_busco.txt", emit: busco_sum
    script:
        """
        if [ -z "${busco_dl_db}" ]
        then
          busco -i ${contigs} -o busco_${deconta} --mode genome -f --auto-lineage-prok
        else
          busco -i ${contigs} -o busco_${deconta} --mode genome \
            --offline --download_path ${busco_dl_db}  -f --auto-lineage-prok
        fi

        mv busco_${deconta}/short_summary*.busco_${deconta}.txt ${id}_${deconta}_busco.txt
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
