include { STAR_GENOMEGENERATE                               } from '../../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN                                        } from '../../modules/nf-core/star/align/main'  
include { SALMON_QUANT                                      } from '../../modules/nf-core/salmon/quant/main'
include { COMBINE_QUANTIFICATION_RESULTS_SALMON             } from '../../modules/local/combine_quantification_results_salmon'
include { SALMON_SPLIT_TABLE as SALMON_SPLIT_TABLE_EACH     } from '../../modules/local/salmon_split_table'
include { SALMON_SPLIT_TABLE as SALMON_SPLIT_TABLE_COMBINED } from '../../modules/local/salmon_split_table'
include { EXTRACT_PROCESSED_READS                           } from '../../modules/local/extract_processed_reads'
include { TXIMPORT                         } from '../../modules/local/tximport/main'
include { COMBINE_QUANTIFICATION_RESULTS_SALMON as COMBINE_QUANTIFICATION_RESULTS_TXIMPORT             } from '../../modules/local/combine_quantification_results_salmon'


workflow SALMON_ALIGNMENT_BASED {

    take:
        ch_reads            // channel: [ val(meta), [ reads ] ]
        ch_host_pathogen_fasta_genome 
        ch_host_pathogen_fasta_transcripts 
        ch_host_pathogen_gff              
        ch_pathogen_fasta_transcripts
        ch_host_fasta_transcripts
        ch_annotations_host_salmon
    main:

        ch_versions = Channel.empty()

        // -------
        // Run create STAR index
        // -------
        STAR_GENOMEGENERATE ( 
            ch_host_pathogen_fasta_genome, 
            ch_host_pathogen_gff 
            )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())


        // -------
        // Run STAR align
        // -------
        STAR_ALIGN ( ch_reads, // reads
                     STAR_GENOMEGENERATE.out.index, // index
                     ch_host_pathogen_gff, // GTF
                     true, //star_ignore_sjdbgtf
                     '', // seq_platform
                     '' // seq_centre
                     )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

      
        // Set to true, as were using alignment-based (with STAR), not selective alignment and Salmon directly
        def alignment_mode = true
        // used to mock the index file, which isnt needed in alignment-mode
        ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true) 

        // -------
        // Run Salmon quant for alignment-based with STAR
        // -------
        SALMON_QUANT(STAR_ALIGN.out.bam_transcript, //reads
                     ch_dummy_file, // dummy file for the index
                     ch_host_pathogen_gff, // GTF
                     ch_host_pathogen_fasta_transcripts, // host / pathogen transcript fasta
                     alignment_mode, // mode
                     params.libtype, // lib type
                     )
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)


        // -------
        // Split the quant table into host and pathogen reads
        // -------
        SALMON_SPLIT_TABLE_EACH(SALMON_QUANT.out.quant, 
                                ch_pathogen_fasta_transcripts, 
                                ch_host_fasta_transcripts
                                )

        // -------
        //  Combine all quant results
        // -------
        // get input file paths (list) of all quant output
        input_files = SALMON_QUANT.out.results.map{it -> it[1]}.collect()

        // Combines all quant results into a file
        COMBINE_QUANTIFICATION_RESULTS_SALMON(
            input_files, 
            Channel.value("both") // combine both host and pathogen results
            )


        // -------
        // Combine all meta data from each datasets
        // -------
        COMBINE_QUANTIFICATION_RESULTS_SALMON.out.combined_quant_data
        .map {it ->
            def meta = [:] // create empty map: meta
            meta.id  = "combined" // sets id with combined
            path_res = it  // assign input (the combined quant path) to path_res
            return [ meta, [ it ] ] // return a tuple containing meta: the metadata map { id: "combined" } and [it]: a list containing the combined quantification file path.
        }.set{ combined_salmon_quant } // set the resulting channel to combined_salmon_quant - containing a tuple [meta, [combined_file_path]]


        // -------
        //  Separate out host and pathogen reads from combined quants
        // -------

        // Generate separate host and pathogen files containing all combined quant results
        // Files: host_quant.sf and pathogen_quant.sf
        SALMON_SPLIT_TABLE_COMBINED( 
            combined_salmon_quant, 
            ch_pathogen_fasta_transcripts, 
            ch_host_fasta_transcripts
            )

        
        // -------
        //  Capture the number of reads processed by Salmon SA and save as output
        // -------
        if (params.mapping_stats) {
            EXTRACT_PROCESSED_READS( 
                SALMON_QUANT.out.json_results, 
                "salmon_alignment" 
                )
        }
        


        // -------
        //  Save gene-level quantifications
        // -------
        TXIMPORT(
            SALMON_SPLIT_TABLE_EACH.out.host, 
            ch_annotations_host_salmon
            )


    emit:
        versions = ch_versions                     // channel: [ versions.yml ]
}


