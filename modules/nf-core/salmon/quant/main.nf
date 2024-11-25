process SALMON_QUANT {
    tag "$meta.id"
    label "process_medium"

    conda "bioconda::salmon=1.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.1--h7e5ed60_0' :
        'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf
    path  transcript_fasta
    val   alignment_mode
    val   lib_type

    output:
    tuple val(meta), path("${prefix}") , emit: results
    tuple val(meta), path("*info.json"), emit: json_info, optional: true
    tuple val(meta), path("${prefix}/quant.sf")  , emit: quant
    path  "versions.yml"               , emit: versions
    tuple val(meta), path("${prefix}/aux_info/meta_info.json") , emit: json_results

    when:
    task.ext.when == null || task.ext.when

    script:
    // set additional args to to additional_params 
    def additional_params = params.salmon_sa_args ?: ''


    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    // command for SA quant, whic passes the index and reads
    def reference   = "--index $index"
    def input_reads = meta.single_end ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    
    // for AB, use transcriptome 
    if (alignment_mode) {
        reference   = "-t $transcript_fasta"
        input_reads = "-a $reads"
        additional_params = params.salmon_ab_args ?: ''
    }

    def strandedness_opts = [
        'A', 'U', 'SF', 'SR',
        'IS', 'IU' , 'ISF', 'ISR',
        'OS', 'OU' , 'OSF', 'OSR',
        'MS', 'MU' , 'MSF', 'MSR'
    ]
    def strandedness =  'A'
    if (lib_type) {
        if (strandedness_opts.contains(lib_type)) {
            strandedness = lib_type
        } else {
            log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to auto-detection with '--libType=A'."
        }
    } else {
        strandedness = meta.single_end ? 'U' : 'IU'
        if (meta.strandedness == 'forward') {
            strandedness = meta.single_end ? 'SF' : 'ISF'
        } else if (meta.strandedness == 'reverse') {
            strandedness = meta.single_end ? 'SR' : 'ISR'
        }
    }
    """
    salmon quant \\
        --geneMap $gtf \\
        --threads $task.cpus \\
        --libType=$strandedness \\
        $reference \\
        $input_reads \\
        $args \\
        $additional_params \\
        -o $prefix

    if [ -f $prefix/aux_info/meta_info.json ]; then
        cp $prefix/aux_info/meta_info.json "${prefix}_meta_info.json"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
