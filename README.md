# Nextflow wrapper for ncov-tools

This is a wrapper that helps to submit [ncov-tools](https://github.com/jts/ncov-tools) jobs to our compute cluster.

## Usage

This pipeline is intended to be run on the output of the [BCCDC-PHL/ncov2019-artic-nf](https://github.com/BCCDC-PHL/ncov2019-artic-nf) pipeline.
The ncov2019-artic-nf output directory is passed as input to this pipeline using the `--artic_analysis_dir` flag.

A `metadata.tsv` file may be optionally provided. If it is provided, it should follow the format as described in [the ncov-tools documentation](https://github.com/jts/ncov-tools#metadata-optional) 

The `completeness_threshold` can optionally be modified. Its default value is the same as that of `ncov-tools`, `0.75`.

The `negative_control_samples` can optionally be modified to override the default search for a sample starting with 'NEG'.

In the command below, `[square brackets]` indicate optional parameters. `<angled brackets>` indicate a placeholder to be replaced by the user in a real pipeline invocation.

```
nextflow run BCCDC-PHL/ncov-tools-nf \
  -profile conda \
  --cache ~/.conda/envs \
  --run_name <run_name> \
  --artic_analysis_dir <output/from/ncov2019-artic-nf> \
  [--completeness_threshold <completeness_threshold>] \
  [--negative_control_samples <negative_control_samples>] \
  [--metadata <metadata.tsv>] \
  --outdir <ncov-tools-output> \
```

## Caveats
- We currently assume that the negative control sample is named starting with the letters `NEG`, and that
  no other sample in the run starts with those letters. This is being addressed in [this issue](https://github.com/BCCDC-PHL/ncov-tools-nf/issues/6).
