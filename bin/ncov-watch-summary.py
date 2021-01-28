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

        if not record.chrom.startswith('#'):
            v = Variant(record.chrom, record.pos, record.ref, record.alts[0])
            if "Name" in record.info:
                v.name = record.info["Name"]
            variants.append(v)
        
    return variants


def parse_ncov_watch_ouput_by_sample_id(ncov_watch_output_path):
    ncov_watch_output_by_sample_id = {}
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
    ncov_watch_output_by_sample_id = parse_ncov_watch_ouput_by_sample_id(args.ncov_watch_output)

    watch_variants = load_vcf(args.watchlist)

    watch_dict = {}
    for v in watch_variants:
        watch_dict[v.key()] = v.name
    
    print("\t".join(["sample_id", "variant_id", "num_observed_substitutions", "total_substitutions_for_variant", "proportion_observed"]))
    
    total_substitutions_for_variant = len(watch_dict.keys())
    
    for sample_id, observed_substitutions in ncov_watch_output_by_sample_id.items():
        num_observed_substitutions = len(observed_substitutions.keys())
        proportion_observed = num_observed_substitutions / total_substitutions_for_variant
        sample_summary_record = '\t'.join([
            sample_id, args.variant_id, str(num_observed_substitutions), str(total_substitutions_for_variant), "{:.3f}".format(proportion_observed)
        ])
        print(sample_summary_record)
    

if __name__ == "__main__":
    description = 'Summarize variants detected per sample and watchlist'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-w', '--watchlist', help='file containing the variants to screen for')
    parser.add_argument('-v', '--variant-id', help='identifier for the variant of interest')
    parser.add_argument('ncov_watch_output', help='output from ncov-watch.py')
    args = parser.parse_args()
    main(args)
    exit(0)

    


    
    for f in gen_func():
        if f.find("variants.tsv") >= 0:
            variants = load_ivar_variants(f)
        else:
            variants = load_vcf(f)
        for v in variants:
            if v.key() in watch_dict:
                print("\t".join([os.path.basename(f), watch_dict[v.key()], v.contig, str(v.position), v.reference, v.alt]))
