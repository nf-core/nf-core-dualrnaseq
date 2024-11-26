process COLLATE_PROCESSED_READS {
    label 'process_high'
    conda "python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'nfcore/dualrnaseq:dev' : 'nfcore/dualrnaseq:dev' }"

    input:
    file partial_results
    val(process)

    output:
    path 'total_processed_reads.tsv'

    // Set publishDir to the process using inputs.process
    publishDir "${params.outdir}/mapping_statistics/${process}/", mode: params.publish_dir_mode

    script:
    """
    # Concatenate all partial result files into the master file
     cat ${partial_results.join(' ')} > total_processed_reads.tsv
    """
}
