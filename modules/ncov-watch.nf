process ncov_watch {

    tag { params.run_name }
    
    publishDir "${params.outdir}/variants", mode: 'copy', pattern: "${run_name}_${variant_id}_all_watchlist_mutations.tsv"

    input:
      tuple path(data_root), val(run_name), val(variant_id), val(threshold), path(watchlist)

    output:
      tuple val(variant_id), path(watchlist), path("${run_name}_${variant_id}_all_watchlist_mutations.tsv")

    script:
      """
      ncov-watch.py --directory ${data_root} --watchlist ${watchlist} > ${run_name}_{variant_id}_all_watchlist_mutations.tmp.tsv 2> /dev/null
      sed -e 's/\\.variants\\.tsv//' ${run_name)_{variant_id}_all_mutations.tmp.tsv > ${run_name}_${variant_id}_all_watchlist_mutations.tsv
      """
}