process download_ncov_tools {
  tag { version }
  executor 'local'
  
  input:
  tuple val(version)
  
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
  tuple path(ncov2010_artic_nf_analysis_dir), path(primer_scheme_dir), path(metadata)
  
  output:
  path("ncov-tools-input", type: 'dir')

  script:
  def metadata = metadata.name != 'NO_FILE' ? "cp ${metadata} ncov-tools-input" : ''
  """
  mkdir ncov-tools-input
  ${metadata}
  cp ${ncov2010_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_makeConsensus/* ncov-tools-input
  cp ${ncov2010_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_readMapping/* ncov-tools-input
  cp ${ncov2010_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_trimPrimerSequences/* ncov-tools-input
  cp ${ncov2010_artic_nf_analysis_dir}/ncovIllumina_sequenceAnalysis_callVariants/* ncov-tools-input
  cp ${primer_scheme_dir}/nCoV-2019.reference.fasta ncov-tools-input
  cp ${primer_scheme_dir}/nCoV-2019.primer.bed ncov-tools-input
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
  publishDir "${params.outdir}", pattern: "config.yaml"
  input:
  tuple val(run_name), val(negative_control_sample), val(metadata)
  
  output:
  path("config.yaml")

  script:
  def metadata = metadata.name != 'NO_FILE' ? "metadata: \\\"{data_root}/metadata.tsv\\\"" : ''
  """
  touch config.yml
  echo "data_root: ncov-tools-input" >> config.yaml
  echo "run_name: ${run_name}" >> config.yaml
  echo "negative_control_samples: [ \\"${negative_control_sample}\\" ]" >> config.yaml
  echo "${metadata}" >> config.yaml
  echo "reference_genome: \\"resources/nCoV-2019.reference.fasta\\"" >> config.yaml
  echo "primer_bed: \\"resources/nCoV-2019.primer.bed\\"" >> config.yaml
  echo "bam_pattern: \\"{data_root}/{sample}.sorted.bam\\"" >> config.yaml
  echo "consensus_pattern: \\"{data_root}/{sample}.primertrimmed.consensus.fa\\"" >> config.yaml
  echo "variants_pattern: \\"{data_root}/{sample}.variants.tsv\\"" >> config.yaml
  echo "platform: illumina" >> config.yaml
  echo "bed_type: unique_amplicons" >> config.yaml
  echo "offset: 0" >> config.yaml
  echo "assign_lineages: true" >> config.yaml
  """
}


process ncov_tools {
  cpus 16
  executor 'sge'
  penv 'smp'
  queue 'all.q'

  publishDir "${params.outdir}", mode: 'copy', pattern: "qc_reports/*_summary_qc.tsv"

  input:
  tuple path(config_yaml), path(data_root), path(resources), path(ncov_tools)
  
  output:
  path("qc_reports/*_summary_qc.tsv")

  script:
  """
  snakemake \
    -s ./ncov-tools/workflow/Snakefile \
    --cores 16 \
    all
  """
}