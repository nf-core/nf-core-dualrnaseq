#!/usr/bin/env nextflow
process EXTRACT_PROCESSED_READS {
    tag "extract_processed_reads_${process}"
    label 'process_high'

    conda "python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nfcore/dualrnaseq:dev' :
        'nfcore/dualrnaseq:dev' }"
   
    input: 
    tuple val(meta), file (json_file)
    val(process)

    output:
    path("${meta.id}.txt"), emit: collect_results

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ ${process} == "salmon" ]; then # for Salmon extract 'num_processed' from meta_info.json file
	    processed=\$(grep "num_processed" ${json_file} | sed 's/num_processed//g'| sed 's/[^a-zA-Z0-9]//g') 
	    echo -e "${prefix}\t\${processed}" > ${prefix}.txt
    elif [ ${process} == "salmon_alignment" ]; then # for Salmon alignment-based mode extract 'num_mapped' from meta_info.json file
	    processed=\$(grep "num_mapped" ${json_file} | sed 's/num_mapped//g'| sed 's/[^a-zA-Z0-9]//g') 
	    echo -e "${prefix}\t\${processed}" > ${prefix}.txt
    elif [ ${process} == "star" ]; then  # for STAR extract "Number of input reads" from *Log.final.out file
	    processed=\$(grep "Number of input reads" ${json_file} | sed 's/Number of input reads//g'| sed 's/[^a-zA-Z0-9]//g')
	    echo -e "${prefix}\t\${processed}" > ${prefix}.txt
    fi
    """
}
