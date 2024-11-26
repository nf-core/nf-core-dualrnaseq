include { CREATE_TRANSCRIPTOME_FASTA } from '../../modules/local/create_transcriptome_fasta'
include { COMBINE_FILES } from '../../modules/local/combine_files'




workflow PREPARE_PATHOGEN_TRANSCRIPTOME {
    take:
        uncompressed_fasta_genome
        uncompressed_gff
    main:

        parameters = Channel.value(
                [
                    params.gene_feature_gff_to_create_transcriptome_pathogen,
                    params.gene_attribute_gff_to_create_transcriptome_pathogen
                ]
        )
        CREATE_TRANSCRIPTOME_FASTA(
                uncompressed_fasta_genome,
                uncompressed_gff,
                parameters
       )

       ch_transcriptome = CREATE_TRANSCRIPTOME_FASTA.out
    

    emit:
        transcriptome = ch_transcriptome

}