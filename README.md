# Nextflow wrapper for ncov-tools

![Tests](https://github.com/BCCDC-PHL/ncov-tools-nf/actions/workflows/push_master.yml/badge.svg)

This is a wrapper that helps to submit [ncov-tools](https://github.com/jts/ncov-tools) jobs to our compute cluster.

## Usage

This pipeline is intended to be run on the output of the [BCCDC-PHL/ncov2019-artic-nf](https://github.com/BCCDC-PHL/ncov2019-artic-nf) pipeline.
The ncov2019-artic-nf output directory is passed as input to this pipeline using the `--artic_analysis_dir` flag.

The default primer scheme is the [Freed *et al.*](https://academic.oup.com/biomethods/article/5/1/bpaa014/5873518) 1200bp scheme (aka. [`V1200`](https://github.com/BCCDC-PHL/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V1200)). Other primer schemes can be selected using the `--primer_scheme_name` flag. Valid primer scheme names correspond to directory names in our [primer scheme repository](https://github.com/BCCDC-PHL/artic-ncov2019/tree/master/primer_schemes/nCoV-2019).

A `metadata.tsv` file may be optionally provided. If it is provided, it should follow the format as described in [the ncov-tools documentation](https://github.com/jts/ncov-tools#metadata-optional) 

The `completeness_threshold` can optionally be modified. Its default value is the same as that of `ncov-tools`, `0.75`.

The `--update_pangolin` flag allows control over whether or not pangolin will be updated prior to running ncov-tools.

In the command below, `[square brackets]` indicate optional parameters. `<angled brackets>` indicate a placeholder to be replaced by the user in a real pipeline invocation.

```
nextflow run BCCDC-PHL/ncov-tools-nf \
  -profile conda \
  --cache ~/.conda/envs \
  [--primer_scheme_name <primer_scheme_name> ] \
  [--primer_scheme_version <primer_scheme_version> ] \
  [--update_pangolin] \
  [--completeness_threshold <completeness_threshold>] \
  [--metadata <metadata.tsv>] \
  [--pre_downsampled] \
  [--no_split_by_plate] \
  [--ivar_consensus] \
  [--ivar_variants] \
  --run_name <run_name> \
  --artic_analysis_dir <output/from/ncov2019-artic-nf> \
  --outdir <ncov-tools-output> \
```


## Caveats
- We currently assume that the negative control sample is named starting with the letters `NEG`, and that
  no other sample in the run starts with those letters. This is being addressed in [this issue](https://github.com/BCCDC-PHL/ncov-tools-nf/issues/6).
