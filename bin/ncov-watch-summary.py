#!/usr/bin/env python

import argparse
import csv
import json
import os
from pathlib import Path
import sys

import pysam


class Variant:
    def __init__(self, contig, position, reference, alt):
        self.contig = contig
        self.position = position
        self.reference = reference
        self.alt = alt
        self.name = None

    def key(self):
        return ",".join([str(x) for x in [self.contig, self.position, self.reference, self.alt]])


def load_vcf(filename):
    variants = list()
    f = pysam.VariantFile(filename,'r')
    for record in f:
        if len(record.alts) > 1:
            sys.stderr.write("Multi-allelic VCF not supported\n")
            sys.exit(1)

        # Skip any line in the .vcf that starts with '#'
        if not record.chrom.startswith('#'):
            v = Variant(record.chrom, record.pos, record.ref, record.alts[0])
            if "Name" in record.info:
                v.name = record.info["Name"]
            variants.append(v)
        
    return variants


def parse_sample_ids(sample_ids_path):
    """
    input: "path/to/sample_ids.tsv"
    output: ["sample-01", "sample-02", ...]
    """
    with open(sample_ids_path, 'r') as f:
        all_sample_ids = [line.strip() for line in f]

    return all_sample_ids


def parse_ncov_watch_ouput_by_sample_id(ncov_watch_output_path, all_sample_ids):
    """
    input: "/path/to/ncov_watch_output.tsv", ["sample-01", "sample-02", ...]
    output: {
              "sample-01": {
                             "MN908947.3,11287,GTCTGGTTTT,G": "ORF1ab:del3675-3677",
                             "MN908947.3,23063,A,T": "S:N501Y"
                           },
              "sample-02": {
                             ...
                           },
              ...
            }
    """
    ncov_watch_output_by_sample_id = {}

    for sample_id in all_sample_ids:
        ncov_watch_output_by_sample_id[sample_id] = {}

    try:
        with open(ncov_watch_output_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for record in reader:
                sample_id = record['sample'].split('.')[0]
                if sample_id not in ncov_watch_output_by_sample_id.keys():
                    ncov_watch_output_by_sample_id[sample_id] = {}
                concatenated_key = ','.join([record['contig'], record['position'], record['reference'], record['alt']])
                ncov_watch_output_by_sample_id[sample_id][concatenated_key] = record['mutation']
    except:
        pass

    return ncov_watch_output_by_sample_id


def main(args):
    all_sample_ids = parse_sample_ids(args.sample_ids)

    ncov_watch_output_by_sample_id = parse_ncov_watch_ouput_by_sample_id(args.ncov_watch_output, all_sample_ids)

    watch_variants = load_vcf(args.watchlist)

    # Only count watchlist mutations with unique names
    unique_mutation_names_in_watchlist = list({var.name for var in watch_variants})
    num_unique_mutation_names_in_watchlist = len(unique_mutation_names_in_watchlist)

    print("\t".join(["sample_id", "mutation_set_id", "num_observed_mutations", "num_watch_mutations_in_mutation_set", "proportion_watch_mutations_observed"]))
    for sample_id, observed_mutations in ncov_watch_output_by_sample_id.items():
        # Only count observed mutations with unique names
        unique_observed_mutation_names = list({mut for mut in observed_mutations.values()})
        num_unique_observed_mutation_names = len(unique_observed_mutation_names)
        proportion_watchlist_mutations_observed = num_unique_observed_mutation_names / num_unique_mutation_names_in_watchlist
        sample_summary_record = '\t'.join([
            sample_id, args.watchlist_id, str(num_unique_observed_mutation_names), str(num_unique_mutation_names_in_watchlist), "{:.3f}".format(proportion_watchlist_mutations_observed)
        ])
        print(sample_summary_record)


if __name__ == "__main__":
    description = 'Summarize mutations detected per sample and watchlist'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-s', '--sample-ids', help='file containing a single list of all sample IDs')
    parser.add_argument('-w', '--watchlist', help='file containing the mutations to screen for (.vcf format)')
    parser.add_argument('-i', '--watchlist-id', help='identifier for the watchlist')
    parser.add_argument('ncov_watch_output', help='output from ncov-watch.py')
    args = parser.parse_args()
    main(args)
