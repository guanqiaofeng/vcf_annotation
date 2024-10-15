/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAGE_INPUT } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { SONG_SCORE_DOWNLOAD } from '../subworkflows/icgc-argo-workflows/song_score_download/main'
include { ENSEMBLVEP_DOWNLOAD } from '../modules/nf-core/ensemblvep/download/main'
include { VCF_ANNOTATE_ENSEMBLVEP_SNPEFF } from '../subworkflows/local/vcf_annotate_ensemblvep_snpeff/main'
include { PAYLOAD_VCFANN } from '../modules/local/payload/main'

// include { SNPEFF_SNPEFF } from '../modules/nf-core/snpeff/snpeff/main'

// include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_SNV } from '../modules/local/payload/main'
// include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_INDEL } from '../modules/local/payload/main'
include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
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

    ch_fasta = Channel.fromPath(params.hg38_ref_fa, checkIfExists: true)
                .map{ path -> [ [id: 'fasta'], path ] }

    // ch_input_all = STAGE_INPUT.out.meta_files.map { entry -> entry + [] }
    ch_input_all = STAGE_INPUT.out.meta_files.map{ entry ->
    // Create a new list that includes the original entry and an empty list
        def newEntry = entry.toList() + [[]] // Convert to a list and append an empty list
        return newEntry
    }

    ch_vep_cache = Channel.fromPath(params.vep_cache, checkIfExists: true)  // from wget https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_merged_vep_112_GRCh38.tar.gz
    ch_vep_extra_files = []
    ch_snpeff_cache = Channel.fromPath(params.snpeff_cache, checkIfExists: true)
                       .map{ path -> [ [id: 'cache'], path ] }
    val_sites_per_chunk = "5000"
    // val_sites_per_chunk = null

    // ch_input_all.subscribe { println(it) }
    // STAGE_INPUT.out.meta_files.subscribe { meta_files ->
    //     println("Data type: ${meta_files.getClass().name}")
    //     println(meta_files)
    // }
    // ch_vep_cache.subscribe { println(it) }
    // ch_snpeff_cache.subscribe { println(it) }

    // ENSEMBLVEP_DOWNLOAD(
    //     STAGE_INPUT.out.meta_files.map { meta, vcf, index ->
    //         // Construct the tuple
    //         tuple(meta, val_vep_genome, val_vep_species, val_vep_cache_version)
    //     }
    // )


    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF(
        ch_input_all,                      // channel: [ val(meta), path(vcf), path(tbi), [path(file1), path(file2)...] ]
        ch_fasta,                    // channel: [ val(meta2), path(fasta) ] (optional)
        params.genome,              //   value: genome to use
        params.species,             //   value: species to use
        params.vep_cache_version,       //   value: cache version to use
        ch_vep_cache,                // channel: [ path(cache) ] (optional)
        ch_vep_extra_files,          // channel: [ path(file1), path(file2)... ] (optional)
        params.snpeff_db,               //   value: the db version to use for snpEff
        ch_snpeff_cache,             // channel: [ val(meta), path(cache) ] (optional)
        params.tools_to_use.split(','),            //   value: a list of tools to use options are: ["ensemblvep", "snpeff"]
        val_sites_per_chunk
    )

    // Combine channels to determine upload status and payload creation
    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF.out.vcf_tbi
    .combine(STAGE_INPUT.out.upRdpc)
    .combine(STAGE_INPUT.out.meta_analysis)
    .map{
        meta,vcf,tbi,upRdpc,metaB,analysis ->
        [
            [
                id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                patient:"${meta.patient}",
                sex:"${meta.sex}",
                sample:"${meta.sample}",
                read_group:"${meta.read_group}",
                data_type:"${meta.data_type}",
                date : "${meta.date}",
                genome_build: "${meta.genome_build}",
                genome_annotation: "${params.genome_annotation}",
                read_groups_count: "${meta.numLanes}",
                study_id : "${meta.study_id}",
                date :"${new Date().format("yyyyMMdd")}",
                upRdpc : upRdpc
            ],[vcf,tbi],analysis
        ]
    }.branch{
        upload : it[0].upRdpc
    }
    .set{ch_payload}

    ch_payload.subscribe { println(it) }

    // Make ALN payload
    PAYLOAD_VCFANN(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
        ch_payload.upload,
        Channel.empty()
        .mix(STAGE_INPUT.out.versions)
        .mix(VCF_ANNOTATE_ENSEMBLVEP_SNPEFF.out.versions)
        .collectFile(name: 'collated_versions.yml')
    )
    ch_versions = ch_versions.mix(PAYLOAD_VCFANN.out.versions)

    // Upload
    SONG_SCORE_UPLOAD(PAYLOAD_VCFANN.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
    ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD.out.versions)


    // SNPEFF_SNPEFF(
    //     STAGE_INPUT.out.meta_files,
    //     snpeff_db,
    //     []
    // )
    // ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions)

    // XML to VCF conversion

    // meta_ch.combine(Channel.fromPath(params.xml)).set{xml_ch}

    // XML_VCF (
    //     xml_ch,
    //     Channel.fromPath(params.hg19_ref_fa, checkIfExists: true),
    //     Channel.fromPath(params.hg19_ref_fai, checkIfExists: true)
    // )
    // ch_versions = ch_versions.mix(XML_VCF.out.versions)

    // // VCF lift over
    // hg38_ref_ch = Channel.fromPath(params.hg38_ref_fa, checkIfExists: true)
    //                         .map{ path -> [ [id: 'fasta'], path ] }

    // hg38_ref_dict = Channel.fromPath(params.hg38_ref_dict, checkIfExists: true)
    //                         .map{ path -> [ [id: 'dict'], path ] }

    // hg19_to_hg38_chain_ch = Channel.fromPath(params.hg19_to_hg38_chain, checkIfExists: true)
    //                         .map{ path -> [ [id: 'chain'], path ] }

    // // SNV \\
    // // lift over
    // PICARD_LIFTOVERVCF_SNV (
    //     XML_VCF.out.snv_vcf,
    //     hg38_ref_dict,
    //     hg38_ref_ch,
    //     hg19_to_hg38_chain_ch
    // )
    // ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF_SNV.out.versions)

    // //Payload Generation
    // PICARD_LIFTOVERVCF_SNV.out.vcf_lifted
    // .combine(PICARD_LIFTOVERVCF_SNV.out.vcf_lifted_index)
    // .combine(meta_ch)
    // .combine(PREP_META.out.updated_experiment_info_tsv)
    // .map{
    //     metaA, vcf, metaB, index, meta, metadata_analysis ->
    //     [
    //         meta, [vcf, index], metadata_analysis
    //     ]
    // }.set{vcf_and_index_snv}

    // PAYLOAD_VARIANT_CALL_SNV (
    //     vcf_and_index_snv,
    //     Channel.empty()
    //     .mix(XML_VCF.out.versions)
    //     .mix(PICARD_LIFTOVERVCF_SNV.out.versions)
    //     .collectFile(name: 'collated_versions.yml')
    // )
    // ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL_SNV.out.versions)

    // // Upload
    // SONG_SCORE_UPLOAD_SNV(PAYLOAD_VARIANT_CALL_SNV.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
    // ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD_SNV.out.versions)

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
