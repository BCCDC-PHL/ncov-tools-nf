process update_pangolin {
  executor 'local'

  input:
  val(should_update)

  output:
  val(true)

  script:
  """
  pangolin --update
  """
}

process download_ncov_tools {

  tag { version }

  executor 'local'
  
  input:
  val(version)
  
  output:
  path("ncov-tools", type: 'dir')

  script:
  """
  wget https://github.com/BCCDC-PHL/ncov-tools/archive/v${version}.tar.gz
  tar -xzf v${version}.tar.gz
  mv ncov-tools-${version} ncov-tools
  """
}

process download_artic_ncov2019 {

  tag { version }

  executor 'local'
  
  input:
  tuple val(version), val(primer_scheme)
  
  output:
  path("resources", type: 'dir')

  script:
  """
  wget https://github.com/BCCDC-PHL/artic-ncov2019/archive/v${version}.tar.gz
  tar -xzf v${version}.tar.gz
  mkdir resources
  cp artic-ncov2019-${version}/primer_schemes/nCoV-2019/${primer_scheme}/nCoV-2019.reference.fasta resources
  cp artic-ncov2019-${version}/primer_schemes/nCoV-2019/${primer_scheme}/nCoV-2019.primer.bed resources
  """
}

process download_ncov_watchlists {

  tag { version }

  executor 'local'
  
  input:
  val(version)
  
  output:
  tuple path("watchlists/watchlists.csv"), path("watchlists", type: 'dir')


  script:
  """
  wget https://github.com/BCCDC-PHL/ncov-watchlists/archive/v${version}.tar.gz
  tar -xzf v${version}.tar.gz
  
  cp -r ncov-watchlists-${version}/watchlists watchlists
  """
}

process index_reference_genome {

  tag { "MN908947.3" }

  executor 'local'
  
  input:
  path(resources)
  
  output:
  path("resources", type: 'dir')

  script:
  """
  samtools faidx resources/nCoV-2019.reference.fasta
  """
}

process get_library_plate_ids {

  tag { params.run_name }

  executor 'local'
  
  input:
  path(artic_analysis_dir)
  
  output:
  stdout

  script:
  """
  tail -n+2 ${artic_analysis_dir}/*.qc.csv | \
  grep -v 'POS*' | grep -v 'NEG*' |
  cut -f 1 -d ',' | cut -f 2 -d '-' | sort | uniq
  """
}


process prepare_data_root {

  tag { params.split_by_plate ? params.run_name + " / " + library_plate_id : params.run_name }

  executor 'local'
  
  input:
  tuple path(ncov2019_artic_nf_analysis_dir), path(primer_scheme_dir), path(metadata), val(library_plate_id)
  
  output:
  tuple val(library_plate_id), path("ncov-tools-input", type: 'dir')

  script:
  def metadata = metadata.name != 'NO_FILE' ? "cp ${metadata} ncov-tools-input" : ''
  def filename_glob = params.split_by_plate ? "*-${library_plate_id}-*" : "*"
  def link_downsampled_bams = params.downsampled ? "ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_downsampleAmplicons/${filename_glob} ." : ''
  def link_freebayes_consensus = params.freebayes_consensus ? "ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_callConsensusFreebayes/${filename_glob}.fa ." : ''
  def link_ivar_consensus = params.freebayes_consensus ? '' : "ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_makeConsensus/${filename_glob}.fa ."
  def link_freebayes_variants = params.freebayes_variants ? "ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_callConsensusFreebayes/${filename_glob}.vcf ." : ''
  """
  mkdir ncov-tools-input
  ${metadata}
  pushd ncov-tools-input
  ${link_ivar_consensus}
  ${link_freebayes_consensus}
  ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_readMapping/${filename_glob} .
  ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_trimPrimerSequences/${filename_glob} .
  ${link_downsampled_bams}
  ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_callVariants/${filename_glob}.tsv .
  ${link_freebayes_variants}
  ln -sfn ../${primer_scheme_dir}/nCoV-2019.reference.fasta .
  ln -sfn ../${primer_scheme_dir}/nCoV-2019.primer.bed .
  popd
  """
}

process create_sample_id_list {

  tag { params.split_by_plate ? params.run_name + " / " + library_plate_id : params.run_name }

  executor 'local'

  input:
  tuple val(library_plate_id), path(data_root)

  output:
  tuple val(library_plate_id), path("sample_id_list.tsv")

  script:
  """
  find ${data_root}/ -name '*.variants.tsv' | xargs -n 1 basename | sed 's/\\.variants\\.tsv//' | sort > sample_id_list.tsv
  """
}

process find_negative_control {

  tag { params.split_by_plate ? params.run_name + " / " + library_plate_id : params.run_name }

  executor 'local'
  
  input:
  tuple val(library_plate_id), path(data_root)
  
  output:
  tuple val(library_plate_id), path("neg_control_sample_id.txt")

  script:
  def filename_glob = params.split_by_plate ? "*-${library_plate_id}-*" : "*"
  """
  find ${data_root}/ -name NEG${filename_glob}.consensus.fa -printf "%f" | cut -d '.' -f 1 > neg_control_sample_id.txt
  """
}

process create_config_yaml {

  tag { params.split_by_plate ? params.run_name + " / " + library_plate_id : params.run_name }

  executor 'local'

  input:
  tuple val(library_plate_id), path(negative_control_sample_id), val(metadata)
  
  output:
  tuple val(library_plate_id), path("config.yaml")

  script:
  def metadata = metadata.name != 'NO_FILE' ? "metadata: \\\"{data_root}/metadata.tsv\\\"" : ''
  def bam_pattern = params.downsampled ? "{data_root}/{sample}.mapped.primertrimmed.downsampled.sorted.bam" : "{data_root}/{sample}.mapped.primertrimmed.sorted.bam"
  def consensus_pattern = params.freebayes_consensus ? "{data_root}/{sample}.consensus.fa" : "{data_root}/{sample}.primertrimmed.consensus.fa"
  def variants_pattern = params.freebayes_variants ? "{data_root}/{sample}.variants.norm.vcf" : "{data_root}/{sample}.variants.tsv"
  def run_name_with_plate = params.split_by_plate ? "${params.run_name}_${library_plate_id}" : "${params.run_name}"
  """
  echo "data_root: ncov-tools-input" >> config.yaml
  echo "run_name: ${run_name_with_plate}" >> config.yaml
  if [[ \$( wc -l < ${negative_control_sample_id} ) -ge 1 ]]; then echo "negative_control_samples: [ \\"\$( cat ${negative_control_sample_id} )\\" ]" >> config.yaml; fi
  echo "${metadata}" >> config.yaml
  echo "reference_genome: \\"resources/nCoV-2019.reference.fasta\\"" >> config.yaml
  echo "primer_bed: \\"resources/nCoV-2019.primer.bed\\"" >> config.yaml
  echo "bam_pattern: \\"${bam_pattern}\\"" >> config.yaml
  echo "consensus_pattern: \\"${consensus_pattern}\\"" >> config.yaml
  echo "variants_pattern: \\"${variants_pattern}\\"" >> config.yaml
  echo "platform: illumina" >> config.yaml
  echo "bed_type: unique_amplicons" >> config.yaml
  echo "offset: 0" >> config.yaml
  echo "completeness_threshold: ${params.completeness_threshold}" >> config.yaml
  echo "assign_lineages: true" >> config.yaml
  """
}


process ncov_tools {

  tag { params.split_by_plate ? params.run_name + " / " + library_plate_id : params.run_name }

  executor 'sge'

  penv 'smp'

  queue 'all.q'

  publishDir "${params.outdir}/by_plate/${library_plate_id}", mode: 'copy', pattern: "config.yaml", enabled: params.split_by_plate
  publishDir "${params.outdir}/by_plate/${library_plate_id}", mode: 'copy', pattern: "bed", enabled: params.split_by_plate
  publishDir "${params.outdir}/by_plate/${library_plate_id}", mode: 'copy', pattern: "lineages/${params.run_name}*", enabled: params.split_by_plate
  publishDir "${params.outdir}/by_plate/${library_plate_id}", mode: 'copy', pattern: "plots", enabled: params.split_by_plate
  publishDir "${params.outdir}/by_plate/${library_plate_id}", mode: 'copy', pattern: "qc_analysis", enabled: params.split_by_plate
  publishDir "${params.outdir}/by_plate/${library_plate_id}", mode: 'copy', pattern: "qc_reports/*.tsv", enabled: params.split_by_plate
  publishDir "${params.outdir}/by_plate/${library_plate_id}", mode: 'copy', pattern: "qc_sequencing", enabled: params.split_by_plate
  publishDir "${params.outdir}/by_plate/${library_plate_id}", mode: 'copy', pattern: "qc_annotation", enabled: params.split_by_plate

  publishDir "${params.outdir}", mode: 'copy', pattern: "config.yaml", enabled: !params.split_by_plate
  publishDir "${params.outdir}", mode: 'copy', pattern: "bed", enabled: !params.split_by_plate
  publishDir "${params.outdir}", mode: 'copy', pattern: "lineages/${params.run_name}*", enabled: !params.split_by_plate
  publishDir "${params.outdir}", mode: 'copy', pattern: "plots", enabled: !params.split_by_plate
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_analysis", enabled: !params.split_by_plate
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_reports/*.tsv", enabled: !params.split_by_plate
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_sequencing", enabled: !params.split_by_plate
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_annotation", enabled: !params.split_by_plate

  input:
  tuple val(library_plate_id), path(config_yaml), path(data_root), path(resources), path(ncov_tools), val(pangolin_updated)
  
  output:
  path("config.yaml")
  path("bed")
  path("lineages/${params.run_name}*_lineage_report.csv"), emit: lineage_report
  path("lineages/${params.run_name}*_pangolin_version.txt"), emit: pangolin_version
  path("plots")
  path("qc_analysis")
  path("qc_reports/*.tsv"), emit: qc_reports
  path("qc_sequencing")
  path("qc_annotation")

  script:
  """
  snakemake -s ./ncov-tools/workflow/Snakefile --cores ${task.cpus} all
  snakemake -s ./ncov-tools/workflow/Snakefile --cores 8 all_qc_annotation
  rm qc_reports/${params.run_name}_*ncov_watch_variants.tsv
  """
}

process ncov_watch {

  tag { params.split_by_plate ? params.run_name + " / " + library_plate_id + " / " + mutation_set_id : params.run_name + " / " + mutation_set_id }
  
  cpus 1

  executor 'local'

  publishDir "${params.outdir}/by_plate/${library_plate_id}/ncov_watch", mode: 'copy', pattern: "${params.run_name}*_ncov_watch_variants.tsv", enabled: params.split_by_plate
  publishDir "${params.outdir}/ncov_watch", mode: 'copy', pattern: "${params.run_name}*_ncov_watch_variants.tsv", enabled: !params.split_by_plate

  input:
  tuple val(library_plate_id), path(data_root), val(mutation_set_id), val(watchlist_filename), path(watchlists_dir)
  
  output:
  tuple val(library_plate_id), val(mutation_set_id), val(watchlist_filename), path(watchlists_dir), path("${params.run_name}_*_ncov_watch_variants.tsv")

  script:
  def variants_output_filename = params.split_by_plate == true ? "${params.run_name}_${library_plate_id}_${mutation_set_id}_ncov_watch_variants.tsv" : "${params.run_name}_${mutation_set_id}_ncov_watch_variants.tsv"
  """
  ncov-watch -d ${data_root} --mutation_set ${watchlists_dir}/${watchlist_filename} | sed 's/\\.variants\\.tsv//' > ${variants_output_filename} 2> /dev/null
  """
}

process combine_ncov_watch_variants {

  tag { params.split_by_plate ? params.run_name + " / " + library_plate_id : params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/by_plate/${library_plate_id}/qc_reports", mode: 'copy', pattern: "${params.run_name}*_ncov_watch_variants.tsv", enabled: params.split_by_plate
  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}*_ncov_watch_variants.tsv", enabled: !params.split_by_plate

  input:
  tuple val(library_plate_id), path(variants)

  output:
  path("${params.run_name}*_ncov_watch_variants.tsv")

  script:
  def run_name_with_plate = params.split_by_plate ? "${params.run_name}_${library_plate_id}" : "${params.run_name}"
  """
  head -qn 1 *_variants.tsv | uniq > header.tsv
  tail -qn+2 *_variants.tsv | sort -k1,1 -k4,4n | uniq > data.tsv
  cat header.tsv data.tsv > ${run_name_with_plate}_ncov_watch_variants.tsv
  """
}

process ncov_watch_summary {

  tag { params.split_by_plate ? params.run_name + " / " + library_plate_id + " / " + watchlist_id : params.run_name + " / " + watchlist_id }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/by_plate/${library_plate_id}/ncov_watch", mode: 'copy', pattern: "${params.run_name}*_${watchlist_id}_ncov_watch_summary.tsv", enabled: params.split_by_plate
  publishDir "${params.outdir}/ncov_watch", mode: 'copy', pattern: "${params.run_name}*_${watchlist_id}_ncov_watch_summary.tsv", enabled: !params.split_by_plate

  input:
  tuple val(library_plate_id), val(watchlist_id), val(watchlist_filename), path(watchlists_dir), path(ncov_watch_output), path(sample_ids)

  output:
  tuple val(library_plate_id), path("${params.run_name}*_${watchlist_id}_ncov_watch_summary.tsv")

  script:
  def ncov_watch_summary_filename = params.split_by_plate ? "${params.run_name}_${library_plate_id}_${watchlist_id}_ncov_watch_summary.tsv" : "${params.run_name}_${watchlist_id}_ncov_watch_summary.tsv"
  """
  ncov-watch-summary.py ${ncov_watch_output} --sample-ids ${sample_ids} --watchlist-id ${watchlist_id} --watchlist ${watchlists_dir}/${watchlist_filename} > ncov_watch_summary_tmp.tsv 2> /dev/null
  head -n 1 ncov_watch_summary_tmp.tsv > header.tsv
  tail -n+2 ncov_watch_summary_tmp.tsv | sort -b -k3,3rn -k1,1 > data_sorted.tsv
  cat header.tsv data_sorted.tsv > ${ncov_watch_summary_filename}
  """
}

process combine_ncov_watch_summaries {

  tag { params.split_by_plate ? params.run_name + " / " + library_plate_id : params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/by_plate/${library_plate_id}/qc_reports", mode: 'copy', pattern: "${params.run_name}*_ncov_watch_summary.tsv", enabled: params.split_by_plate
  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}*_ncov_watch_summary.tsv", enabled: !params.split_by_plate

  input:
  tuple val(library_plate_id), path(summaries)

  output:
  path("${params.run_name}*_ncov_watch_summary.tsv")

  script:
  def ncov_watch_summary_filename = params.split_by_plate ? "${params.run_name}_${library_plate_id}_ncov_watch_summary.tsv" : "${params.run_name}_ncov_watch_summary.tsv"
  """
  head -qn 1 *_summary.tsv | uniq > header.tsv
  tail -qn+2 *_summary.tsv | sort -k1,1 -k2,2 > data.tsv
  cat header.tsv data.tsv > ${ncov_watch_summary_filename}
  """
}

process combine_all_qc_summaries_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}_summary_qc.tsv"

  input:
  path(summaries)

  output:
  path("${params.run_name}_summary_qc.tsv")

  script:
  """
  head -qn 1 *_summary_qc.tsv | uniq > header.tsv
  tail -qn+2 *_summary_qc.tsv | sort -k1,1 -k2,2 > data.tsv
  cat header.tsv data.tsv > "${params.run_name}_summary_qc_incorrect_run_name.tsv"
  tail -qn+2 "${params.run_name}_summary_qc_incorrect_run_name.tsv" | awk -F '\t' 'BEGIN {OFS=FS}; {\$2 = substr(\$2, 1, length(\$2) - 4); print}' > data.tsv
  cat header.tsv data.tsv > "${params.run_name}_summary_qc.tsv"
  """
}

process combine_all_mixture_reports_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}_mixture_report.tsv"

  input:
  path(mixture_reports)

  output:
  path("${params.run_name}_mixture_report.tsv")

  script:
  """
  head -qn 1 *_mixture_report.tsv | uniq > header.tsv
  tail -qn+2 *_mixture_report.tsv | sort -k1,1 -k2,2 > data.tsv
  cat header.tsv data.tsv > "${params.run_name}_mixture_report.tsv"
  """
}

process combine_ambiguous_position_reports_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}_ambiguous_position_report.tsv"

  input:
  path(ambiguous_position_reports)

  output:
  path("${params.run_name}_ambiguous_position_report.tsv")

  script:
  """
  head -qn 1 *_ambiguous_position_report.tsv | uniq > header.tsv
  tail -qn+2 *_ambiguous_position_report.tsv | sort -k1,1 -k2,2 > data.tsv
  cat header.tsv data.tsv > "${params.run_name}_ambiguous_position_report.tsv"
  """
}

process combine_all_ncov_watch_summaries_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}_ncov_watch_summary.tsv"

  input:
  path(ncov_watch_summaries)

  output:
  path("${params.run_name}_ncov_watch_summary.tsv")

  script:
  """
  head -qn 1 *_summary.tsv | uniq > header.tsv
  tail -qn+2 *_summary.tsv | sort -k1,1 -k2,2 > data.tsv
  cat header.tsv data.tsv > "${params.run_name}_ncov_watch_summary.tsv"
  """
}

process combine_all_ncov_watch_variants_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}_ncov_watch_variants.tsv"

  input:
  path(ncov_watch_variants)

  output:
  path("${params.run_name}_ncov_watch_variants.tsv")

  script:
  """
  head -qn 1 *_variants.tsv | uniq > header.tsv
  tail -qn+2 *_variants.tsv | sort -k1,1 -k4,4n > data.tsv
  cat header.tsv data.tsv > "${params.run_name}_ncov_watch_variants.tsv"
  """
}

process combine_all_lineage_reports_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/lineages", mode: 'copy', pattern: "${params.run_name}_lineage_report.csv"

  input:
  path(lineage_reports)

  output:
  path("${params.run_name}_lineage_report.csv")

  script:
  """
  head -qn 1 *_lineage_report.csv | uniq > header.csv
  tail -qn+2 *_lineage_report.csv | sort -k1,1 -k2,2 > data.csv
  cat header.csv data.csv > "${params.run_name}_lineage_report.csv"
  """
}

process get_pangolin_version_for_run {

  tag { params.run_name }

  cpus 1

  executor 'local'

  publishDir "${params.outdir}/lineages", mode: 'copy', pattern: "${params.run_name}_pangolin_version.txt"

  input:
  path(lineage_reports)

  output:
  path("${params.run_name}_pangolin_version.txt")

  script:
  """
  cat ${params.run_name}*_pangolin_version.txt > ${params.run_name}_pangolin_version.txt
  """
}
