process emblmygff3 {
    label 'emblmygff3'
    publishDir "${params.output}/emblmygff3", mode: 'copy'
    input:
        path annotation
    output:
        file *.embl
    script:
        """
        # define proper emblmygff3 command
				...
        """
}
