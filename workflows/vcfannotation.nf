/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAGE_INPUT } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { SONG_SCORE_DOWNLOAD } from '../subworkflows/icgc-argo-workflows/song_score_download/main'

include { SNPEFF_SNPEFF } from '../modules/nf-core/snpeff/snpeff/main'

// include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_SNV } from '../modules/local/payload/main'
// include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_INDEL } from '../modules/local/payload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UPLOAD_SNV } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UPLOAD_INDEL } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCFANN {

    println "API token: ${params.api_token}"

    ch_versions = Channel.empty()

    // Validate input, generate metadata, prepare fastq channel
    STAGE_INPUT(
        params.study_id,
        params.analysis_id,
        params.samplesheet
        )

    ch_versions = ch_versions.mix(STAGE_INPUT.out.versions)

    snpeff_db = "GRCh38.105" // latest database in snpeff version 4.1

    // SNPEFF_SNPEFF(
    //     STAGE_INPUT.out.meta_files,
    //     snpeff_db
    // )
    // ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions)

    // XML to VCF conversion

    meta_ch.combine(Channel.fromPath(params.xml)).set{xml_ch}

    XML_VCF (
        xml_ch,
        Channel.fromPath(params.hg19_ref_fa, checkIfExists: true),
        Channel.fromPath(params.hg19_ref_fai, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(XML_VCF.out.versions)

    // VCF lift over
    hg38_ref_ch = Channel.fromPath(params.hg38_ref_fa, checkIfExists: true)
                            .map{ path -> [ [id: 'fasta'], path ] }

    hg38_ref_dict = Channel.fromPath(params.hg38_ref_dict, checkIfExists: true)
                            .map{ path -> [ [id: 'dict'], path ] }

    hg19_to_hg38_chain_ch = Channel.fromPath(params.hg19_to_hg38_chain, checkIfExists: true)
                            .map{ path -> [ [id: 'chain'], path ] }

    // SNV \\
    // lift over
    PICARD_LIFTOVERVCF_SNV (
        XML_VCF.out.snv_vcf,
        hg38_ref_dict,
        hg38_ref_ch,
        hg19_to_hg38_chain_ch
    )
    ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF_SNV.out.versions)

    //Payload Generation
    PICARD_LIFTOVERVCF_SNV.out.vcf_lifted
    .combine(PICARD_LIFTOVERVCF_SNV.out.vcf_lifted_index)
    .combine(meta_ch)
    .combine(PREP_META.out.updated_experiment_info_tsv)
    .map{
        metaA, vcf, metaB, index, meta, metadata_analysis ->
        [
            meta, [vcf, index], metadata_analysis
        ]
    }.set{vcf_and_index_snv}

    PAYLOAD_VARIANT_CALL_SNV (
        vcf_and_index_snv,
        Channel.empty()
        .mix(XML_VCF.out.versions)
        .mix(PICARD_LIFTOVERVCF_SNV.out.versions)
        .collectFile(name: 'collated_versions.yml')
    )
    ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL_SNV.out.versions)

    // Upload
    SONG_SCORE_UPLOAD_SNV(PAYLOAD_VARIANT_CALL_SNV.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
    ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD_SNV.out.versions)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.dump_parameters(workflow, params)
//     NfcoreTemplate.summary(workflow, params, log)
//     if (params.hook_url) {
//         NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//     }
// }

// workflow.onError {
//     if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
//         println("🛑 Default resources exceed availability 🛑 ")
//         println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
