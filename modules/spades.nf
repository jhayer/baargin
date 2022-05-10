process spades {
    label 'spades'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(illumina)
    output:
        tuple val(id), path("${id}_scaffolds.fasta"), emit: assembly
        tuple val(id), path(illumina), path("${id}_scaffolds.fasta"), emit: quast
        path "${id}_spades.log"                    , emit: log
        path "${id}_contigs.fasta"                 , emit: contigs
        path "${id}_graph.gfa"                     , emit: graph
        path "spades.version.txt"                  , emit: version
    script:
        """
        spades.py -1 ${illumina[0]} -2 ${illumina[1]} \
         -o spades --isolate -t ${task.cpus}

        mv spades/assembly_graph_with_scaffolds.gfa ${id}_graph.gfa
        mv spades/scaffolds.fasta ${id}_scaffolds.fasta
        mv spades/contigs.fasta ${id}_contigs.fasta
        mv spades/spades.log ${id}_spades.log

        spades.py --version | sed "s/SPAdes v//; s/ \\[.*//" > spades.version.txt
        """
}
