#!/usr/bin/env python3

import argparse
import csv
import glob
import json
import os

def parse_qc_summary(qc_summary):
    """
    Parse the qc_summary file into a dictionary.

    :param qc_summary: Path to the qc_summary file
    :type qc_summary: str
    :return: Dictionary of qc_summary data. Keys are sample IDs.
    :rtype: dict
    """
    qc_summary_by_sample_id = {}
    int_fields = [
        'num_consensus_snvs',
        'num_consensus_n',
        'num_consensus_iupac',
        'num_variants_snvs',
        'num_variants_indel',
        'num_variants_indel_triplet',
        'median_sequencing_depth',
    ]
    float_fields = [
        'mean_sequencing_depth',
        'qpcr_ct',
        'genome_completeness',
    ]
    with open(qc_summary, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample_id = row['sample']
            if sample_id not in qc_summary_by_sample_id:
                for field in int_fields:
                    try:
                        row[field] = int(row[field])
                    except ValueError as e:
                        row[field] = None
                for field in float_fields:
                    try:
                        row[field] = float(row[field])
                    except ValueError as e:
                        row[field] = None
                for k, v in row.items():
                    if v == 'NA' or v == '':
                        row[k] = None
                row['qc_pass'] = row['qc_pass'].split(',')
                qc_summary_by_sample_id[sample_id] = row

    return qc_summary_by_sample_id


def check_expected_summary_qc_pass_values(qc_summary_by_sample_id, expected_qc_pass_by_sample_id):
    """
    Check that the expected QC pass values are present in the summary QC files.

    :param qc_summary_by_sample_id: Summary QC data by sample ID
    :type qc_summary_by_sample_id: dict
    :param expected_qc_pass_by_sample_id: Expected QC pass values by sample ID
    :type expected_qc_pass_by_sample_id: dict
    :return: Whether the expected QC pass values are present in the summary QC files
    :rtype: bool
    """
    for sample_id, expected_qc_pass in expected_qc_pass_by_sample_id.items():
        if sample_id not in qc_summary_by_sample_id:
            return False
        for test_qc_pass_value in qc_summary_by_sample_id[sample_id]['qc_pass']:
            if test_qc_pass_value not in expected_qc_pass:
                return False

    return True


def main(args):
    expected_qc_pass_by_sample_id = {
        'POS-nCoVWGS-1-25x': ['PASS'],
        'SRR27382851-1-25x': ['PASS'],
        'SRR27382852-1-25x': ['PASS'],
        'POS-nCoVWGS-2-25x': ['POSSIBLE_FRAMESHIFT_INDELS'],
        'SRR27503679-2-25x': ['PASS'],
        'SRR27503680-2-25x': ['POSSIBLE_FRAMESHIFT_INDELS', 'EXCESS_AMBIGUITY'],
    }

    summary_qc_glob = os.path.join(
        args.pipeline_outdir,
        'by_plate',
        '*',
        'qc_reports',
        '*_summary_qc.tsv',
    )
    summary_qc_files = glob.glob(summary_qc_glob)
    qc_summary_by_sample_id = {}
    for summary_qc_file in summary_qc_files:
        file_qc_summary_by_sample_id = parse_qc_summary(summary_qc_file)
        qc_summary_by_sample_id.update(file_qc_summary_by_sample_id)
        
    tests = [
        {
            "test_name": "expected_summary_qc_pass_values",
            "test_passed": check_expected_summary_qc_pass_values(qc_summary_by_sample_id, expected_qc_pass_by_sample_id),
        },
    ]

    output_fieldnames = [
        'test_name',
        'test_passed',
    ]

    with open(args.output, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=output_fieldnames, extrasaction='ignore')
        writer.writeheader()
        for test in tests:
            writer.writerow(test)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--pipeline-outdir', type=str, help='Path to the pipeline output directory')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    args = parser.parse_args()
    main(args)
