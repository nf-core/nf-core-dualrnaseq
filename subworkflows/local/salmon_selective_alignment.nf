include { SALMON_INDEX                                      } from '../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT                                      } from '../../modules/nf-core/salmon/quant/main'
include { COMBINE_QUANTIFICATION_RESULTS_SALMON             } from '../../modules/local/combine_quantification_results_salmon'
include { SALMON_SPLIT_TABLE as SALMON_SPLIT_TABLE_EACH     } from '../../modules/local/salmon_split_table'
include { SALMON_SPLIT_TABLE as SALMON_SPLIT_TABLE_COMBINED } from '../../modules/local/salmon_split_table'
include { EXTRACT_PROCESSED_READS                           } from '../../modules/local/extract_processed_reads'
include { COLLATE_PROCESSED_READS                           } from '../../modules/local/collate_processed_reads'
include { TXIMPORT                                          } from '../../modules/local/tximport/main'




workflow SALMON_SELECTIVE_ALIGNMENT {

    take:
        ch_reads            // channel: [ val(meta), [ reads ] ]
        ch_host_pathogen_fasta_genome     // channel: /path/to/host_pathogen_fasta_genome.fasta
        ch_host_pathogen_fasta_transcripts 
        ch_host_pathogen_gff              
        ch_pathogen_fasta_transcripts
        ch_host_fasta_transcripts
        ch_annotations_host_salmon
        
    main:
        ch_versions = Channel.empty()

        // -------
        // Run salmon index
        // -------
        ch_salmon_index = SALMON_INDEX ( ch_host_pathogen_fasta_genome, 
                                         ch_host_pathogen_fasta_transcripts
                                         ).index
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

        // Set to false, as were using selective alignment (not with Star and alignment directly)
        def alignment_mode = false


        // -------
        // Run Salmon quant for selective alignment
        // -------
        SALMON_QUANT(ch_reads, 
                     ch_salmon_index, 
                     ch_host_pathogen_gff, 
                     ch_host_pathogen_fasta_transcripts, 
                     alignment_mode, 
                     params.libtype
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

        // combine all meta data from each datasets
        COMBINE_QUANTIFICATION_RESULTS_SALMON.out.combined_quant_data
        .map {it ->
            def meta = [:] // create empty map: meta
            meta.id  = "combined" // sets id with combined
            path_res = it // assign input (the combined quant path) to path_res
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
            ch_host_fasta_transcripts)




        // -------
        //  Capture the number of reads processed by Salmon SA and save as output
        // -------
        if (params.mapping_stats) {

            // Is saved under: mapping_statistics/salmon_SA
            // Separate files for each sample
            EXTRACT_PROCESSED_READS( 
                SALMON_QUANT.out.json_results, 
                "salmon" // Tell the script to just look at Salmon SA output 
            )

            // Store the read count summary files from each quant run
            EXTRACT_PROCESSED_READS.out.collect_results
                .collect()
                .set { collected_processed_reads_files }

            // Merge all individual results into a single file
            COLLATE_PROCESSED_READS(collected_processed_reads_files)

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
