include {
    UNZIPFILES as UNCOMPRESS_HOST_FASTA_GENOME;
    UNZIPFILES as UNCOMPRESS_PATHOGEN_FASTA_GENOME;
    UNZIPFILES as UNCOMPRESS_HOST_GFF;
    UNZIPFILES as UNCOMPRESS_PATHOGEN_GFF;
    UNZIPFILES as UNCOMPRESS_HOST_TRANSCRIPTOME;
    UNZIPFILES as UNCOMPRESS_PATHOGEN_TRANSCRIPTOME
} from '../../modules/nf-core/unzipfiles/main'

include {
    REPLACE_ATTRIBUTE_GFF_STAR_SALMON as REPLACE_ATTRIBUTE_GFF_STAR_SALMON_HOST;
    REPLACE_ATTRIBUTE_GFF_STAR_SALMON as REPLACE_ATTRIBUTE_GFF_STAR_SALMON_PATHOGEN;
} from '../../modules/local/replace_attribute'

include {
    REPLACE_GENE_FEATURE_GFF_SALMON as REPLACE_GENE_FEATURE_GFF_PATHOGEN_SALMON;
    REPLACE_GENE_FEATURE_GFF_SALMON as REPLACE_GENE_FEATURE_GFF_HOST_SALMON
 } from '../../modules/local/replace_gene_feature'

include {
    // COMBINE_FILES as COMBINE_PATHOGEN_HOST_GFF_FILES_HTSEQ;
    COMBINE_FILES as COMBINE_FILES_PATHOGEN_HOST_GFF;
    COMBINE_FILES as COMBINE_FILES_FASTA;
    COMBINE_FILES as COMBINE_FILES_TRANSCRIPTOME_FILES;
}  from '../../modules/local/combine_files'


include { PREPARE_HOST_TRANSCRIPTOME      } from './prepare_host_transcriptome'
include { PREPARE_PATHOGEN_TRANSCRIPTOME  } from './prepare_pathogen_transcriptome'

include {
    // EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_HOST_HTSEQ;
    EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_HOST_SALMON;
    // EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_PATHOGEN_HTSEQ;
    EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_PATHOGEN_SALMON
} from '../../modules/local/extract_annotations'


workflow PREPARE_REFERENCE_FILES{
  take:
    host_fasta_genome //fasta_host
    host_gff //gff_host
    pathogen_fasta_genome //fasta_pathogen
    pathogen_gff //gff_pathogen

  main:
    ch_transcriptome = Channel.empty()
    ch_host_transcriptome = Channel.empty()
    ch_pathogen_transcriptome = Channel.empty()
    ch_host_pathogen_gff = Channel.empty()

    ch_gene_feature_pathogen = Channel
	    .value(params.gene_feature_gff_to_create_transcriptome_pathogen)
	    .collect()
    
    



    // -------------------
    // uncompress fasta files and gff files
    // -------------------

    // Capture files
    ch_host_fasta_genome = Channel.value(file(params.host_fasta_genome, checkIfExists: true))
    ch_host_gff = Channel.value(file(params.host_gff, checkIfExists: true))
    ch_pathogen_gff = Channel.value(file(params.pathogen_gff, checkIfExists: true))
    ch_pathogen_fasta_genome = Channel.value(file(params.pathogen_fasta_genome, checkIfExists: true)) 

    // Uncompress pathogen (genome) fasta if needed
    if (params.pathogen_fasta_genome.endsWith('.gz') || params.pathogen_fasta_genome.endsWith('.zip')){
        ch_pathogen_fasta_genome_unzipped = UNCOMPRESS_PATHOGEN_FASTA_GENOME(ch_pathogen_fasta_genome)
    } else {
        ch_pathogen_fasta_genome_unzipped = ch_pathogen_fasta_genome
    }

    // Uncompress pathogen gff if needed
    if (params.pathogen_gff.endsWith('.gz') || params.pathogen_gff.endsWith('.zip')){
        ch_pathogen_gff_unzipped = UNCOMPRESS_PATHOGEN_GFF(ch_pathogen_gff)
    } else {
        ch_pathogen_gff_unzipped = ch_pathogen_gff
    }

    // Uncompress host fasta if needed
    if (params.host_fasta_genome.endsWith('.gz') || params.host_fasta_genome.endsWith('.zip')){
        ch_host_fasta_genome_unzipped = UNCOMPRESS_HOST_FASTA_GENOME(ch_host_fasta_genome)
    } else {
        ch_host_fasta_genome_unzipped = ch_host_fasta_genome
    }

    // Uncompress host gff if needed
    if (params.host_gff.endsWith('.gz') || params.host_gff.endsWith('.zip')){
        ch_host_gff_unzipped = UNCOMPRESS_HOST_GFF(ch_host_gff)
    } else {
        ch_host_gff_unzipped = ch_host_gff
    }



    // -------------------
    // Combine pathogen and host fastas
    // -------------------

    // Combine both host and pathogen fasta files and save
    COMBINE_FILES_FASTA(ch_pathogen_fasta_genome_unzipped, ch_host_fasta_genome_unzipped, 'host_pathogen_genome.fasta' )



    //
    // Execute steps specific for mapping-quantification modes
    //


    // Salmon SA or Salmon AB (transcriptome-based)
    if(params.run_salmon_SA | params.run_salmon_AB) {
      
      // HOST - Has a host transcriptome (fasta) been passed? 
      if(params.host_fasta_transcripts){
          
          // Grab the host transcriptome file
          ch_host_fasta_transcripts  = params.host_fasta_transcripts ? Channel.value(file(params.host_fasta_transcripts, checkIfExists: true )) : Channel.empty()

          // Uncompress if needed
          if (params.host_fasta_transcripts.endsWith('.gz') || params.host_fasta_transcripts.endsWith('.zip')){
                  UNCOMPRESS_HOST_TRANSCRIPTOME(ch_host_fasta_transcripts)
                  ch_host_fasta_transcripts_unzipped = UNCOMPRESS_HOST_TRANSCRIPTOME.out.files
          } else {
                  ch_host_fasta_transcripts_unzipped = ch_host_fasta_transcripts
          }

      } else {
      // No host transcriptome passed - prepare      
      PREPARE_HOST_TRANSCRIPTOME(
        ch_host_fasta_genome_unzipped,
        ch_host_gff_unzipped,
      )
      ch_host_fasta_transcripts_unzipped = PREPARE_HOST_TRANSCRIPTOME.out.transcriptome
      }

      // PATHOGEN - Has a pathogen transcriptome (fasta) been passed? 
      if(params.pathogen_fasta_transcripts){
          
          // Grab the transcriptome file
          ch_pathogen_fasta_transcripts  = params.pathogen_fasta_transcripts ? Channel.value(file( params.pathogen_fasta_transcripts, checkIfExists: true )) : Channel.empty()
          
          // Uncompress if needed
          if (params.pathogen_fasta_transcripts.endsWith('.gz') || params.pathogen_fasta_transcripts.endsWith('.zip')){
                    UNCOMPRESS_PATHOGEN_TRANSCRIPTOME(ch_pathogen_fasta_transcripts)
                    ch_pathogen_fasta_transcripts_unzipped = UNCOMPRESS_PATHOGEN_TRANSCRIPTOME.out.files
          } else {
                    ch_pathogen_fasta_transcripts_unzipped = ch_pathogen_fasta_transcripts
          }
        
        } else {
          // No host transcriptome passed - prepare   
          PREPARE_PATHOGEN_TRANSCRIPTOME(
            ch_pathogen_fasta_transcripts_unzipped,
            ch_pathogen_gff_unzipped
          )
          ch_pathogen_fasta_transcripts_unzipped = PREPARE_PATHOGEN_TRANSCRIPTOME.out.transcriptome
        }



      // combine host and pathogen transcriptome
      COMBINE_FILES_TRANSCRIPTOME_FILES(
        ch_host_fasta_transcripts_unzipped,ch_pathogen_fasta_transcripts_unzipped,'host_pathogen_transcripts.fasta'
      )
      ch_combined_fasta_transcripts = COMBINE_FILES_TRANSCRIPTOME_FILES.out



      // --------------
      // Replace attributes in the reference files 
      // --------------
    
      // ---
      // HOST
      // ---

      // Replace selected attributes in gff to parent
      // this changes the gene identifier value such as locus_tag to parent
      REPLACE_ATTRIBUTE_GFF_STAR_SALMON_HOST(
                  ch_host_gff_unzipped,
                  'Parent',
                  'parent'
      )

      // Replace selected gene features in gff - save file with extension "_quant_feature.gff3"
      // this changes the selected values from the third column in the GFF to quant
      // since the input file is from replace attribute, this resulting file now has both
      // gene attribute and gene feature changed. This is useful as the host and pathogen
      // annotations may have different identifiers, such as locus_tag for bacteria and gene_id
      // for human/mouse - and when combining these annotations we want to have common a 
      // naming convention in the final file
      ch_host_genome_gff_salmon_sa = REPLACE_ATTRIBUTE_GFF_STAR_SALMON_HOST.out
      REPLACE_GENE_FEATURE_GFF_HOST_SALMON(
                    ch_host_genome_gff_salmon_sa,
                    params.gene_feature_gff_to_create_transcriptome_host
      )

          

      // ---
      // PATHOGEN
      // ---

      // Replace selected attributes in gff to parent
      // this changes the gene identifier value such as locus_tag to parent
      REPLACE_ATTRIBUTE_GFF_STAR_SALMON_PATHOGEN(
            ch_pathogen_gff_unzipped,
            params.gene_attribute_gff_to_create_transcriptome_pathogen,
            'parent'
      )

      // Replace selected gene features in gff
      // this changes the selected values from the third column in the GFF to quant
      // since the input file is from replace attribute, this resulting file now has both
      // gene attribute and gene feature changed. This is useful as the host and pathogen
      // annotations may have different identifiers, such as locus_tag for bacteria and gene_id
      // for human/mouse - and when combining these annotations we want to have common a 
      // naming convention in the final file
      REPLACE_GENE_FEATURE_GFF_PATHOGEN_SALMON(
            REPLACE_ATTRIBUTE_GFF_STAR_SALMON_PATHOGEN.out,
            ch_gene_feature_pathogen
      )

          


      // ---
      // COMBINED
      // ---
      
      // Combine gff files with replaced features and save
      COMBINE_FILES_PATHOGEN_HOST_GFF(
            REPLACE_GENE_FEATURE_GFF_PATHOGEN_SALMON.out,
            REPLACE_GENE_FEATURE_GFF_HOST_SALMON.out,
            "host_pathogen.gff"
      )


      
      // ---
      // Extracting the GFF annotations into a .tsv file for downstream analysis
      // ---

      // Extract host features
      EXTRACT_ANNOTATIONS_HOST_SALMON (
                REPLACE_GENE_FEATURE_GFF_HOST_SALMON.out,
                'quant',
                'parent',
                params.host_organism,
                'salmon'
      )

      // Extract pathogen featues
      EXTRACT_ANNOTATIONS_PATHOGEN_SALMON (
              REPLACE_ATTRIBUTE_GFF_STAR_SALMON_PATHOGEN.out,
              ch_gene_feature_pathogen,
              "parent",
              params.pathogen_organism,
              'salmon'
        )




    } // end --> if(params.run_salmon_SA | params.run_salmon_AB) {


    emit:
      host_pathogen_fasta_genome = COMBINE_FILES_FASTA.out // 'host_pathogen_genome.fasta'
      host_pathogen_fasta_transcripts = ch_combined_fasta_transcripts // 'host_pathogen_transcripts.fasta'
      host_fasta_transcripts = ch_host_fasta_transcripts_unzipped
      pathogen_fasta_transcripts = ch_pathogen_fasta_transcripts_unzipped
      host_pathogen_gff = COMBINE_FILES_PATHOGEN_HOST_GFF.out // 'host_pathogen.gff'
      annotations_host_salmon = EXTRACT_ANNOTATIONS_HOST_SALMON.out.annotations // extracted_annotations_host_salmon.tsv
      annotations_pathogen_salmon = EXTRACT_ANNOTATIONS_PATHOGEN_SALMON.out.annotations // extracted_annotations_pathogen_salmon.tsv
    }

