process SALMON_INDEX {
    tag "$transcript_fasta"
    label "process_medium"

    conda "bioconda::salmon=1.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.1--h7e5ed60_0' :
        'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0' }"

    input:
    path genome_fasta
    path transcript_fasta

    output:
    path "salmon_index"      , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // set additional args to to additional_params 
    def additional_params = params.salmon_sa_index_args ?: ''
    


    def get_decoy_ids = "grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
    def gentrome      = "gentrome.fa"
    if (genome_fasta.endsWith('.gz')) {
        get_decoy_ids = "grep '^>' <(gunzip -c $genome_fasta) | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
        gentrome      = "gentrome.fa.gz"
    }

   

    """
    $get_decoy_ids
    sed -i.bak -e 's/>//g' decoys.txt
    cat $transcript_fasta $genome_fasta > $gentrome

    # Debug: Print the salmon index command to logs
    echo "Running salmon index with:"
    echo "  Threads: $task.cpus"
    echo "  Gentrome: $gentrome"
    echo "  Decoy file: decoys.txt"
    echo "  Additional params: $additional_params"

    salmon \\
        index \\
        --threads $task.cpus \\
        -t $gentrome \\
        -d decoys.txt \\
        $args \\
        $additional_params \\
        -i salmon_index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
