process COLLATE_PROCESSED_READS {
    label 'process_high'
    conda "python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'nfcore/dualrnaseq:dev' : 'nfcore/dualrnaseq:dev' }"

    input:
    file partial_results

    output:
    path 'total_processed_reads.tsv'

    script:
    """
    # Concatenate all partial result files into the master file
     cat ${partial_results.join(' ')} > total_processed_reads.tsv
    """
}
