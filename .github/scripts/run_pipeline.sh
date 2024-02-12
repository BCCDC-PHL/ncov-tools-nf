#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH

BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION="1.3.4"

# the github runner only has 4 cpus available
#sed -i s'/cpus = 14/cpus = 4/'g nextflow.config
#sed -i s'/--cores 14/--cores 4/'g modules/ncov-tools.nf
#sed -i s'/--cores 8/--cores 4/'g modules/ncov-tools.nf

nextflow run main.nf \
	 -profile conda \
	 --cache ~/.conda/envs \
	 --artic_analysis_dir test_input/ncov2019-artic-nf-v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION}-output \
	 --metadata .github/data/metadata/metadata.tsv \
	 --run_name test \
	 --outdir results \
	 -resume \
	 -with-trace ncov-tools-nf_trace.tsv

mv .nextflow.log artifacts/nextflow.log
cp -r results artifacts/results
cp -r work artifacts/work
