include { SALMON_INDEX                                      } from '../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT                                      } from '../../modules/nf-core/salmon/quant/main'
include { COMBINE_QUANTIFICATION_RESULTS_SALMON             } from '../../modules/local/combine_quantification_results_salmon'
include { SALMON_SPLIT_TABLE as SALMON_SPLIT_TABLE_EACH     } from '../../modules/local/salmon_split_table'
include { SALMON_SPLIT_TABLE as SALMON_SPLIT_TABLE_COMBINED } from '../../modules/local/salmon_split_table'
include { EXTRACT_PROCESSED_READS                           } from '../../modules/local/extract_processed_reads'
include { TXIMPORT                                          } from '../../modules/local/tximport/main'



workflow SALMON_SELECTIVE_ALIGNMENT {

    take:
        ch_reads            // channel: [ val(meta), [ reads ] ]
        ch_genome_fasta     // channel: /path/to/genome.fasta
        ch_transcript_fasta // channel: /path/to/transcript.fasta
        ch_gtf              // channel: /path/to/genome.gtf
        ch_transcript_fasta_pathogen
        ch_transcript_fasta_host
        ch_annotations_host_salmon
        
    main:
        ch_versions = Channel.empty()
        ch_salmon_index = SALMON_INDEX ( ch_genome_fasta, 
                                         ch_transcript_fasta
                                         ).index
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

        def alignment_mode = false
        SALMON_QUANT(ch_reads, 
                     ch_salmon_index, 
                     ch_gtf, 
                     ch_transcript_fasta, 
                     alignment_mode, 
                     params.libtype
                     )
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

        SALMON_SPLIT_TABLE_EACH(SALMON_QUANT.out.quant, 
                                ch_transcript_fasta_pathogen, 
                                ch_transcript_fasta_host)

        input_files = SALMON_QUANT.out.results.map{it -> it[1]}.collect()

        COMBINE_QUANTIFICATION_RESULTS_SALMON(input_files, Channel.value("both"))


        COMBINE_QUANTIFICATION_RESULTS_SALMON.out.combined_quant_data
        .map {it ->
            def meta = [:]
            meta.id  = "combined"
            path_res = it
            return [ meta, [ it ] ]
        }.set{ combined_salmon_quant }

        SALMON_SPLIT_TABLE_COMBINED( combined_salmon_quant, ch_transcript_fasta_pathogen, ch_transcript_fasta_host)


        EXTRACT_PROCESSED_READS( SALMON_QUANT.out.json_results, "salmon" )

        TXIMPORT(SALMON_SPLIT_TABLE_EACH.out.host,ch_annotations_host_salmon )
        
    emit:
        versions = ch_versions                     // channel: [ versions.yml ]
}
