#!/usr/bin/env python

import argparse
import pysam
import sys
import csv
import os
from pathlib import Path

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

        
        v = Variant(record.chrom, record.pos, record.ref, record.alts[0])
        if "Name" in record.info:
            v.name = record.info["Name"]
        variants.append(v)
        
    return variants

def load_ivar_variants(filename):
    variants = list()
    try:
        with open(filename, 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for record in reader:
                ref = record["REF"]
                alt = record["ALT"]
                if alt[0] == "-":
                    ref += alt[1:]
                    alt = ref[0]
                elif alt[0] == "+":
                    alt = ref + alt[1:]

                variants.append(Variant(record["REGION"], record["POS"], ref, alt))
    except:
        pass
    return variants

def get_from_directory():
    find_files = lambda pattern : [ path for path in Path(args.directory).rglob(pattern) ]
    files = find_files("*pass.vcf") + find_files("*pass.vcf.gz") + find_files("*variants.tsv")
    for f in files:
        yield str(f)

def get_from_stdin():
    for line in sys.stdin:
        yield line.rstrip()

if __name__ == "__main__":
    description = 'Report samples containing a variant in the watchlist'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-w', '--watchlist', help='file containing the variants to screen for')
    parser.add_argument('-d', '--directory', help='root of directories holding variant files')
    args = parser.parse_args()

    watch_variants = load_vcf(args.watchlist)
    watch_dict = dict()
    for v in watch_variants:
        watch_dict[v.key()] = v.name

    gen_func = get_from_stdin
    if args.directory:
        gen_func = get_from_directory

    print("\t".join(["sample", "mutation", "contig", "position", "reference", "alt"]))
    for f in gen_func():
        if f.find("variants.tsv") >= 0:
            variants = load_ivar_variants(f)
        else:
            variants = load_vcf(f)
        for v in variants:
            if v.key() in watch_dict:
                print("\t".join([os.path.basename(f), watch_dict[v.key()], v.contig, str(v.position), v.reference, v.alt]))
