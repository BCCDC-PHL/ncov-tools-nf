#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH
export PATH=/opt/nextflow/bin:$PATH

BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION="1.3.2"

# write test log as github Action artifact
echo Nextflow run BCCDC-PHL/ncov2019-artic-nf to generate input... >> artifacts/test_artifact.log
NXF_VER=21.04.3 nextflow pull BCCDC-PHL/ncov2019-artic-nf -r v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION} 

mkdir test_input && pushd test_input

NXF_VER=21.04.3 nextflow -quiet run BCCDC-PHL/ncov2019-artic-nf \
       -r v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION} \
       -profile conda \
       --cache ~/.conda/envs \
       --directory $PWD/../.github/data/fastqs/ \
       --ref $PWD/../.github/data/refs/MN908947.3/MN908947.3.fa \
       --bed $PWD/../.github/data/primer_schemes/nCoV-2019_Freed_1200bp.bed \
       --primer_pairs_tsv $PWD/../.github/data/primer_schemes/nCoV-2019_Freed_1200bp_primer_pairs.tsv \
       --gff $PWD/../.github/data/refs/MN908947.3.gff \
       --composite_ref $PWD/../.github/data/refs/mock_composite_ref/mock_composite_ref.fa \
       --illumina \
       --prefix test \
       --outdir ncov2019-artic-nf-v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION}-output

mv .nextflow.log ../artifacts/ncov2019-artic-nf.nextflow.log
cp -r ncov2019-artic-nf-v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION}-output ../artifacts/

popd

# the github runner only has 4 cpus available
sed -i s'/cpus = 14/cpus = 4/'g nextflow.config
sed -i s'/--cores 14/--cores 4/'g modules/ncov-tools.nf
sed -i s'/--cores 8/--cores 4/'g modules/ncov-tools.nf

nextflow -quiet run main.nf \
	 -profile conda \
	 --cache ~/.conda/envs \
	 --artic_analysis_dir test_input/ncov2019-artic-nf-v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION}-output \
	 --run_name test \
	 --outdir results

mv .nextflow.log artifacts/nextflow.log
cp -r results artifacts/results
cp -r work artifacts/work
