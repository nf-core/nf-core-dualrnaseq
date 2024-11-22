include { CREATE_TRANSCRIPTOME_FASTA_GFFREAD } from '../../modules/local/create_transcriptome_fasta_gffread'
// include { CREATE_TRANSCRIPTOME_FASTA as COMBINE_HOST_GFF_FILES } from '../../modules/local/create_transcriptome_fasta'
// include { CREATE_TRANSCRIPTOME_FASTA as CREATE_TRANSCRIPTOME_TRNA_FASTA } from '../../modules/local/create_transcriptome_fasta'
// include { COMBINE_FILES as COMBINE_TRANSCRIPT_FASTA_HOST }  from '../../modules/local/combine_files'
// include { REPLACE_GENE_FEATURE_GFF_SALMON } from '../../modules/local/replace_gene_feature'



// Prepare host transcriptome
workflow PREPARE_HOST_TRANSCRIPTOME {
  take:
    uncompressed_fasta_genome
    uncompressed_gff_host
    // uncompressed_gff_host_tRNA

  main:
    
        // Create host transcriptome using GFFREAD 
        CREATE_TRANSCRIPTOME_FASTA_GFFREAD(
            uncompressed_fasta_genome,
            uncompressed_gff_host
        )  

        // Store out file
        ch_transcriptome = CREATE_TRANSCRIPTOME_FASTA_GFFREAD.out

        // NOTE:
        // Code to be removed - related to including host tRNA file - excluded now
        // if(params.gff_host_tRNA){
        //  //   UNCOMPRESS_GFF_TRNA(ch_gff_host)
        // //    UNCOMPRESS_GFF_TRNA_FILE(ch_gff_host_tRNA)
        //     parameters = Channel.value(
        //         [
        //             params.gene_feature_gff_to_create_transcriptome_host,
        //             params.gene_attribute_gff_to_create_transcriptome_host
        //         ]
        //     )
        //     CREATE_TRANSCRIPTOME_TRNA_FASTA(
        //         uncompressed_fasta_genome,
        //         uncompressed_gff_host_tRNA,
        //         parameters
        //     )

        //     // COMBINE_TRANSCRIPT_FASTA_HOST(CREATE_TRANSCRIPTOME_FASTA_GFFREAD.out, CREATE_TRANSCRIPTOME_TRNA_FASTA.out, 'host_transcriptome.fasta' )
        //     COMBINE_TRANSCRIPT_FASTA_HOST(CREATE_TRANSCRIPTOME_FASTA_GFFREAD.out, 'host_transcriptome.fasta' )


        //     ch_transcriptome = COMBINE_TRANSCRIPT_FASTA_HOST.out
        //  //   ch_uncompress_gff_trna_file = UNCOMPRESS_GFF_TRNA_FILE.out
        //  //   ch_uncompress_gff_trna = UNCOMPRESS_GFF_TRNA.out
        //  //   ch_uncompress_gff = Channel.empty()
        // } else {
        //  //   ch_uncompress_gff_trna_file = Channel.empty()
        //  //   ch_uncompress_gff_trna = Channel.empty()
        //  //   ch_uncompress_gff = UNCOMPRESS_GFF.out
        //     ch_transcriptome = CREATE_TRANSCRIPTOME_FASTA_GFFREAD.out
        // }
        
 emit:
      transcriptome = ch_transcriptome
 //   uncompress_gff_trna_file = ch_uncompress_gff_trna_file
 //   uncompress_gff_trna = ch_uncompress_gff_trna
  //  uncompress_gff = ch_uncompress_gff

}