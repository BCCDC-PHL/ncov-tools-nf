#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { download_ncov_tools } from './modules/ncov-tools.nf'
include { download_artic_ncov2019 } from './modules/ncov-tools.nf'
include { download_ncov_watchlists } from './modules/ncov-tools.nf'
include { index_reference_genome } from './modules/ncov-tools.nf'
include { prepare_data_root } from './modules/ncov-tools.nf'
include { find_negative_control } from './modules/ncov-tools.nf'
include { create_config_yaml } from './modules/ncov-tools.nf'
include { ncov_tools } from './modules/ncov-tools.nf'
include { ncov_watch } from './modules/ncov-tools.nf'

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
  // closure below is used to transform [[1, 2], 3] into [1, 2, 3]. Is there a better way?
  ch_watchlists = download_ncov_watchlists.out.splitCsv().map{ it -> [it[0][0], it[0][1], it[1]] }
  index_reference_genome(download_artic_ncov2019.out)
  prepare_data_root(ch_artic_analysis_dir.combine(download_artic_ncov2019.out).combine(ch_metadata))
  find_negative_control(prepare_data_root.out)
  create_config_yaml(ch_run_name.combine(find_negative_control.out).combine(ch_metadata))
  ncov_tools(create_config_yaml.out.combine(prepare_data_root.out).combine(index_reference_genome.out).combine(download_ncov_tools.out))
  ncov_watch(prepare_data_root.out.combine(ch_watchlists))
  
}
