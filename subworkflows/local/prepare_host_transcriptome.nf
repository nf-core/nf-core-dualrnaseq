include { CREATE_TRANSCRIPTOME_FASTA_GFFREAD } from '../../modules/local/create_transcriptome_fasta_gffread'

// Prepare host transcriptome
workflow PREPARE_HOST_TRANSCRIPTOME {
  take:
    uncompressed_fasta_genome
    uncompressed_gff_host


  main:
    
        // Create host transcriptome using GFFREAD 
        CREATE_TRANSCRIPTOME_FASTA_GFFREAD(
            uncompressed_fasta_genome,
            uncompressed_gff_host
        )  

        // Store out file
        ch_transcriptome = CREATE_TRANSCRIPTOME_FASTA_GFFREAD.out

        
 emit:
      transcriptome = ch_transcriptome


}