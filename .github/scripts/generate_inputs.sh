#!/bin/bash

set -eo pipefail


BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION="1.3.4"

mkdir -p artifacts

# write test log as github Action artifact
echo Nextflow run BCCDC-PHL/ncov2019-artic-nf to generate input... >> artifacts/test_artifact.log
NXF_VER=21.04.3 nextflow pull BCCDC-PHL/ncov2019-artic-nf -r v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION} 

mkdir -p test_input && pushd test_input

NXF_VER=21.04.3 nextflow run BCCDC-PHL/ncov2019-artic-nf \
       -r v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION} \
       -profile conda \
       --cache ${HOME}/.conda/envs \
       --downsampleMinDepth 25 \
       --directory ../.github/data/fastqs/ \
       --ref ${PWD}/../.github/data/refs/MN908947.3/MN908947.3.fa \
       --bed ${PWD}/../.github/data/primer_schemes/nCoV-2019_Freed_1200bp.bed \
       --primer_pairs_tsv ${PWD}/../.github/data/primer_schemes/nCoV-2019_Freed_1200bp_primer_pairs.tsv \
       --gff ${PWD}/../.github/data/refs/MN908947.3.gff \
       --composite_ref ${PWD}/../.github/data/refs/mock_composite_ref/mock_composite_ref.fa \
       --illumina \
       --prefix test \
       --outdir ncov2019-artic-nf-v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION}-output \
       -with-trace artic_nextflow_trace.tsv \
       -with-report artic_nextflow_report.html

mv .nextflow.log ../artifacts/ncov2019-artic-nf.nextflow.log
cp -r ncov2019-artic-nf-v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION}-output ../artifacts/

popd
