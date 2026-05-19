#!/usr/bin/env python3
"""Merge variant calls from multiple callers."""

import argparse

def main():
    parser = argparse.ArgumentParser(description='Merge VCF files from multiple callers')
    parser.add_argument('--input', nargs='+', required=True, help='Input VCF files')
    parser.add_argument('--output', required=True, help='Output VCF file')
    parser.add_argument('--min-callers', type=int, default=1, help='Minimum callers supporting a variant')
    args = parser.parse_args()

    variants = {}
    for vcf_file in args.input:
        with open(vcf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                key = (parts[0], parts[1], parts[3], parts[4].split(',')[0])
                variants.setdefault(key, 0)
                variants[key] += 1

    with open(args.output, 'w') as out:
        out.write('##fileformat=VCFv4.2\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for (chrom, pos, ref, alt), count in sorted(variants.items()):
            if count >= args.min_callers:
                out.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tCALLERS={count}\n')

if __name__ == '__main__':
    main()
