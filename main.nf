#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { download_ncov_tools } from './modules/ncov-tools.nf'
include { download_artic_ncov2019 } from './modules/ncov-tools.nf'
include { download_ncov_watchlists } from './modules/ncov-tools.nf'
include { index_reference_genome } from './modules/ncov-tools.nf'
include { prepare_data_root } from './modules/ncov-tools.nf'
include { create_sample_id_list } from './modules/ncov-tools.nf'
include { find_negative_control } from './modules/ncov-tools.nf'
include { create_config_yaml } from './modules/ncov-tools.nf'
include { ncov_tools } from './modules/ncov-tools.nf'
include { ncov_watch } from './modules/ncov-tools.nf'
include { combine_ncov_watch_variants } from './modules/ncov-tools.nf'
include { ncov_watch_summary } from './modules/ncov-tools.nf'
include { combine_ncov_watch_summaries } from './modules/ncov-tools.nf'

workflow {
  
  ch_primer_scheme_version = Channel.of(params.primer_scheme_version)
  ch_primer_scheme_name = Channel.of('V1200')
  ch_ncov_tools_version = Channel.of(params.ncov_tools_version)
  ch_ncov_watchlists_version = Channel.of(params.ncov_watchlists_version)
  ch_run_name = Channel.of(params.run_name)
  ch_artic_analysis_dir = Channel.fromPath(params.artic_analysis_dir, type: 'dir')
  ch_metadata = Channel.fromPath(params.metadata, type: 'file')
  
  download_ncov_tools(ch_ncov_tools_version)
  download_artic_ncov2019(ch_primer_scheme_version.combine(ch_primer_scheme_name))
  download_ncov_watchlists(ch_ncov_watchlists_version)
  ch_watchlists = download_ncov_watchlists.out.splitCsv().map{ it -> [it[0][0], it[0][1], it[1]] }
  index_reference_genome(download_artic_ncov2019.out)
  
  ch_library_plate_ids = Channel.fromList([194, 195])
  prepare_data_root(ch_artic_analysis_dir.combine(download_artic_ncov2019.out).combine(ch_metadata).combine(ch_library_plate_ids))
  create_sample_id_list(prepare_data_root.out)
  find_negative_control(prepare_data_root.out)
  create_config_yaml(ch_run_name.combine(ch_library_plate_ids).combine(find_negative_control.out).combine(ch_metadata))
  ncov_tools(create_config_yaml.out.join(prepare_data_root.out).combine(index_reference_genome.out).combine(download_ncov_tools.out))
  ncov_watch(prepare_data_root.out.combine(ch_watchlists))
  combine_ncov_watch_variants(ncov_watch.out.map{ it -> [it[0], it[4]] }.groupTuple())
  ncov_watch_summary(create_sample_id_list.out.cross(ncov_watch.out).map{ it -> it[1] + it[0][1].toString() })
  combine_ncov_watch_summaries(ncov_watch_summary.out.groupTuple())
  
}
