process unicycler {
    label 'unicycler'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(illuminaR1), path(illuminaR2), path(ont)
    output:
        tuple val(id), path("${id}_scaffolds.fasta"), emit: assembly
    //    tuple val(id), path(illumina), path("${id}_scaffolds.fasta"), emit: quast
        path "${id}_unicycler.log"  , emit: log
        path "${id}_graph.gfa"      , emit: graph
        path "unicycler.version.txt", emit: version
    script:
        """

        unicycler -1 ${illuminaR1} -2 ${illuminaR2} -l ${ont} -o unicycler \
          -t ${task.cpus} --keep 0
        mv unicycler/assembly.fasta ${id}_scaffolds.fasta
        mv unicycler/assembly.gfa ${id}_graph.gfa
        mv unicycler/unicycler.log ${id}_unicycler.log

        unicycler --version > unicycler.version.txt
        """
}
