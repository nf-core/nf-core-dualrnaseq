/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDualrnaseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config,
    // host
    params.host_fasta_genome,
    params.host_gff,
    // params.gff_host_tRNA,
    // pathogen
    params.pathogen_fasta_genome,
    params.pathogen_gff,
]
// loop through and check all params specified above
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// Check mandatory parameters - just samplesheet required at this point
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// configuration channels

//for multiqc
ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

//subworkflow and module inclusion
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPARE_REFERENCE_FILES } from '../subworkflows/local/prepare_reference_files'
include { SALMON_SELECTIVE_ALIGNMENT } from '../subworkflows/local/salmon_selective_alignment'
include { SALMON_ALIGNMENT_BASED } from '../subworkflows/local/salmon_alignment_based'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                            } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_AFTER_TRIMMING   } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                          } from '../modules/nf-core/cutadapt/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DUALRNASEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // if skip_tools passed, but not contain fastqc
    if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
            FASTQC(INPUT_CHECK.out.reads)
            ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    // Store input reads
    // ch_reads = INPUT_CHECK.out.reads // prev code
    ch_reads = INPUT_CHECK.out.reads.map { meta, reads -> tuple(meta, reads) }

    // if skip_tools passed, but not contain cutadapt
    if (!(params.skip_tools && params.skip_tools.split(',').contains('cutadapt'))) {
        CUTADAPT(ch_reads) 
        ch_reads = CUTADAPT.out.reads
        ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
    }

    // if skip_tools passed, but not contain fastqc and cutadapt - so should run fastqc after trimming
    if (!(params.skip_tools && (params.skip_tools.split(',').contains('fastqc') || params.skip_tools.split(',').contains('cutadapt')))) {
            FASTQC_AFTER_TRIMMING(ch_reads)
            ch_versions = ch_versions.mix(FASTQC_AFTER_TRIMMING.out.versions.first())
    }


    // ---------------
    // Prepare reference files
    // ---------------
    // uncompress files, merge host and pathogen files together, 
    // update features in reference files to be compatable with software
    PREPARE_REFERENCE_FILES(
        params.host_fasta_genome,
        params.host_gff,
        params.pathogen_fasta_genome,
        params.pathogen_gff
    )


    // Run Salmon selective alighment
    if ( params.run_salmon_SA ) {
        SALMON_SELECTIVE_ALIGNMENT (
            ch_reads,
            PREPARE_REFERENCE_FILES.out.host_pathogen_fasta_genome,
            PREPARE_REFERENCE_FILES.out.host_pathogen_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.host_pathogen_gff,
            PREPARE_REFERENCE_FILES.out.pathogen_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.host_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.annotations_host_salmon
        )
        ch_versions = ch_versions.mix(SALMON_SELECTIVE_ALIGNMENT.out.versions)
    }

    // Run Salmon alignment based
    if ( params.run_salmon_AB ) {
        SALMON_ALIGNMENT_BASED (
            ch_reads,
            PREPARE_REFERENCE_FILES.out.host_pathogen_fasta_genome,
            PREPARE_REFERENCE_FILES.out.host_pathogen_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.host_pathogen_gff,
            PREPARE_REFERENCE_FILES.out.pathogen_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.host_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.annotations_host_salmon
        )
        ch_versions = ch_versions.mix(SALMON_ALIGNMENT_BASED.out.versions)
    }


    // Capture software versions
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    // MODULE: MultiQC
    workflow_summary    = WorkflowDualrnaseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowDualrnaseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()


}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
