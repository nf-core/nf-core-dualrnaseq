process EXTRACT_ANNOTATIONS {
    tag "${organism}_${quantifier}"
    label 'process_high'

    conda "bioconda::conda-forge::python=3.8.3=3.11.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nfcore/dualrnaseq:dev':
        'nfcore/dualrnaseq:dev' }"

    input:
        path(gff)
        val(gene_feature)
        val(gene_attribute)
        val(organism)    // host or pathogen
        val(quantifier)  // htseq or salmon

    output:
        path("*.tsv")      , emit: annotations
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "extracted_annotations_${organism}_${quantifier}"

    """
    python $workflow.projectDir/bin/extract_annotations_from_gff.py  \\
        --gff $gff  \\
        --gene_feature $gene_feature  \\
        --gene_attribute $gene_attribute  \\
        --organism $organism  \\
        --quantifier $quantifier  \\
        --output $prefix  \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{ print \$2 }')
    END_VERSIONS
    """
}
