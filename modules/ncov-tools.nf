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

process prepare_data_root {
  executor 'local'
  
  input:
  tuple path(ncov2019_artic_nf_analysis_dir), path(primer_scheme_dir), path(metadata)
  
  output:
  path("ncov-tools-input", type: 'dir')

  script:
  def metadata = metadata.name != 'NO_FILE' ? "cp ${metadata} ncov-tools-input" : ''
  """
  mkdir ncov-tools-input
  ${metadata}
  pushd ncov-tools-input
  ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_makeConsensus/* .
  ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_readMapping/* .
  ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_trimPrimerSequences/* .
  ln -sfn ../${ncov2019_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_callVariants/* .
  ln -sfn ../${primer_scheme_dir}/nCoV-2019.reference.fasta .
  ln -sfn ../${primer_scheme_dir}/nCoV-2019.primer.bed .
  popd
  """
}

process create_sample_id_list {
  executor 'local'

  input:
  path(data_root)

  output:
  path("sample_id_list.tsv")

  script:
  """
  find ${data_root}/ -name '*.variants.tsv' | xargs -n 1 basename | sed 's/\\.variants\\.tsv//' | sort > sample_id_list.tsv
  """
}

process find_negative_control {
  executor 'local'
  
  input:
  path(data_root)
  
  output:
  stdout

  script:
  """
  find ${data_root}/ -type f -name 'NEG*.primertrimmed.consensus.fa' -printf "%f" | cut -d '.' -f 1 | tr -d \$'\n'
  """
}

process create_config_yaml {
  executor 'local'

  input:
  tuple val(run_name), val(negative_control_sample), val(metadata)
  
  output:
  path("config.yaml")

  script:
  def metadata = metadata.name != 'NO_FILE' ? "metadata: \\\"{data_root}/metadata.tsv\\\"" : ''
  """
  echo "data_root: ncov-tools-input" >> config.yaml
  echo "run_name: ${run_name}" >> config.yaml
  if [ "${negative_control_sample}" != "" ]; then echo "negative_control_samples: [ \\"${negative_control_sample}\\" ]" >> config.yaml; fi
  echo "${metadata}" >> config.yaml
  echo "reference_genome: \\"resources/nCoV-2019.reference.fasta\\"" >> config.yaml
  echo "primer_bed: \\"resources/nCoV-2019.primer.bed\\"" >> config.yaml
  echo "bam_pattern: \\"{data_root}/{sample}.sorted.bam\\"" >> config.yaml
  echo "consensus_pattern: \\"{data_root}/{sample}.primertrimmed.consensus.fa\\"" >> config.yaml
  echo "variants_pattern: \\"{data_root}/{sample}.variants.tsv\\"" >> config.yaml
  echo "platform: illumina" >> config.yaml
  echo "bed_type: unique_amplicons" >> config.yaml
  echo "offset: 0" >> config.yaml
  echo "completeness_threshold: ${params.completeness_threshold}" >> config.yaml
  echo "assign_lineages: true" >> config.yaml
  """
}


process ncov_tools {

  tag { params.run_name }
  
  cpus 16
  executor 'sge'
  penv 'smp'
  queue 'all.q'

  publishDir "${params.outdir}", mode: 'copy', pattern: "config.yaml"
  publishDir "${params.outdir}", mode: 'copy', pattern: "bed"
  publishDir "${params.outdir}", mode: 'copy', pattern: "lineages"
  publishDir "${params.outdir}", mode: 'copy', pattern: "plots"
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_analysis"
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_reports/*.tsv"
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_sequencing"
  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_annotation"

  input:
  tuple path(config_yaml), path(data_root), path(resources), path(ncov_tools)
  
  output:
  path("config.yaml")
  path("bed")
  path("lineages")
  path("plots")
  path("qc_analysis")
  path("qc_reports/*.tsv")
  path("qc_sequencing")
  path("qc_annotation")

  script:
  """
  snakemake -s ./ncov-tools/workflow/Snakefile --cores 16 all
  snakemake -s ./ncov-tools/workflow/Snakefile --cores 2 all_qc_annotation
  rm qc_reports/${params.run_name}_ncov_watch_variants.tsv
  """
}

process ncov_watch {

  tag { mutation_set_id }
  
  cpus 1
  executor 'local'

  publishDir "${params.outdir}/ncov_watch", mode: 'copy', pattern: "${params.run_name}_${mutation_set_id}_ncov_watch_variants.tsv"

  input:
  tuple path(data_root), val(mutation_set_id), val(watchlist_filename), path(watchlists_dir)
  
  output:
  tuple val(mutation_set_id), val(watchlist_filename), path(watchlists_dir), path("${params.run_name}_${mutation_set_id}_ncov_watch_variants.tsv")

  script:
  """
  ncov-watch -d ${data_root} --mutation_set ${watchlists_dir}/${watchlist_filename} | sed 's/\\.variants\\.tsv//' > ${params.run_name}_${mutation_set_id}_ncov_watch_variants.tsv 2> /dev/null
  """
}

process combine_ncov_watch_variants {

  tag { params.run_name }

  cpus 1
  executor 'local'

  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}_ncov_watch_variants.tsv"

  input:
  path(variants)

  output:
  path("${params.run_name}_ncov_watch_variants.tsv")

  script:
  """
  head -qn 1 *_variants.tsv | uniq > header.tsv
  tail -qn+2 *_variants.tsv | sort -k1,1 -k4,4n | uniq > data.tsv
  cat header.tsv data.tsv > ${params.run_name}_ncov_watch_variants.tsv
  """
}

process ncov_watch_summary {

  tag { watchlist_id }

  cpus 1
  executor 'local'

  publishDir "${params.outdir}/ncov_watch", mode: 'copy', pattern: "${params.run_name}_${watchlist_id}_ncov_watch_summary.tsv"

  input:
  tuple val(watchlist_id), val(watchlist_filename), path(watchlists_dir), path(ncov_watch_output), path(sample_ids)

  output:
  path("${params.run_name}_${watchlist_id}_ncov_watch_summary.tsv")

  script:
  """
  ncov-watch-summary.py ${ncov_watch_output} --sample-ids ${sample_ids} --watchlist-id ${watchlist_id} --watchlist ${watchlists_dir}/${watchlist_filename} > ncov_watch_summary_tmp.tsv 2> /dev/null
  head -n 1 ncov_watch_summary_tmp.tsv > header.tsv
  tail -n+2 ncov_watch_summary_tmp.tsv | sort -b -k3,3rn -k1,1 > data_sorted.tsv
  cat header.tsv data_sorted.tsv > ${params.run_name}_${watchlist_id}_ncov_watch_summary.tsv
  """
}

process combine_ncov_watch_summaries {

  tag { params.run_name }

  cpus 1
  executor 'local'

  publishDir "${params.outdir}/qc_reports", mode: 'copy', pattern: "${params.run_name}_ncov_watch_summary.tsv"

  input:
  path(summaries)

  output:
  path("${params.run_name}_ncov_watch_summary.tsv")

  script:
  """
  head -qn 1 *_summary.tsv | uniq > header.tsv
  tail -qn+2 *_summary.tsv | sort -k1,1 -k2,2 > data.tsv
  cat header.tsv data.tsv > ${params.run_name}_ncov_watch_summary.tsv
  """
}
