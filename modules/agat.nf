process agat_sp_complement_annotations.pl {
    label 'agat'
    publishDir "${params.output}/meta_annotation", mode: 'copy'
    input:
        path(annotations)
    output:
        file("meta_annotation.gff3")
    script:
        """
				# prepare parameters for agat script
				command=""
				foreach annotation $annotations;do
					command="${command} --gff \${annotation}"
				done
				# run the complement agat script
				agat_sp_complement_annotations.pl $command -o meta_annotation.gff3
        """
}
