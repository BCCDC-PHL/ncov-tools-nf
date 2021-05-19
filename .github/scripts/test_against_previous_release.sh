#!/bin/bash

# set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH
export PATH=/opt/nextflow/bin:$PATH

BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION="1.3.2"

# write test log as github Action artifact
echo Nextflow run BCCDC-PHL/ncov2019-artic-nf to generate input... >> artifacts/test_artifact.log
NXF_VER=20.10.0 nextflow pull BCCDC-PHL/ncov2019-artic-nf -r v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION} 

mkdir ncov2019-artic-nf-v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION}-output && pushd ncov2019-artic-nf-v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION}-output

NXF_VER=20.10.0 nextflow -quiet run BCCDC-PHL/ncov2019-artic-nf \
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
       --outdir .

mv .nextflow.log ../artifacts/ncov2019-artic-nf.nextflow.log

popd

# the github runner only has 2 cpus available
sed -i s'/cpus = 14/cpus = 2/'g nextflow.config
sed -i s'/--cores 14/--cores 2/'g modules/ncov-tools.nf
sed -i s'/--cores 8/--cores 2/'g modules/ncov-tools.nf

echo Nextflow run this pull-request... >> artifacts/test_artifact.log
NXF_VER=20.10.0 nextflow -quiet run main.nf \
       -profile conda \
       --cache ~/.conda/envs \
       --artic_analysis_dir $PWD/ncov2019-artic-nf-v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION}-output \
       --run_name test \
       --outdir results \
       --downsampled \
       --freebayes_consensus \
       --freebayes_variants

mv .nextflow.log artifacts/pull_request.nextflow.log
cp -r results artifacts/pull_request_results
cp -r work artifacts/pull_request_work

# run tests against previous previous_release to compare outputs 
git clone https://github.com/BCCDC-PHL/ncov-tools-nf.git previous_release 
pushd previous_release
git checkout 75a41e7b239fb1f61ab2d1df6c70d35a09c5c745

# the github runner only has 2 cpus available
sed -i s'/cpus = 14/cpus = 2/'g nextflow.config
sed -i s'/--cores 14/--cores 2/'g modules/ncov-tools.nf
sed -i s'/--cores 8/--cores 2/'g modules/ncov-tools.nf

echo Nextflow run previous release >> ../artifacts/test_artifact.log
NXF_VER=20.10.0 nextflow -quiet run main.nf \
       -profile conda \
       --cache ~/.conda/envs \
       --artic_analysis_dir $PWD/../ncov2019-artic-nf-v${BCCDC_NCOV2019_ARTIC_PIPELINE_VERSION}-output \
       --run_name test \
       --outdir results \
       --downsampled \
       --freebayes_consensus \
       --freebayes_variants

mv .nextflow.log ../artifacts/previous_release.nextflow.log
cp -r results ../artifacts/previous_release_results
cp -r work ../artifacts/previous_release_work

popd

# exclude files from comparison
# and list differences
echo "Compare ouputs of current PR vs those of previous release.." >> artifacts/test_artifact.log
find results ./previous_release/results \
     -name "*.fq.gz" \
     -o -name "*.bam" \
     -o -name "*.bam.bai" \
     -o -name "*.vcf" \
    | xargs rm -rf
if ! git diff --stat --no-index results ./previous_release/results > diffs.txt ; then
  echo "test failed: differences found between PR and previous release" >> artifacts/test_artifact.log
  echo "see diffs.txt" >> artifacts/test_artifact.log 
  cp diffs.txt artifacts/  
  exit 1
else
  echo "no differences found between PR and previous release" >> artifacts/test_artifact.log
fi

# clean-up for following tests
rm -rf previous_release && rm -rf results && rm -rf work && rm -rf .nextflow*
