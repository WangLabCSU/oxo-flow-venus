#!/usr/bin/env python
"""Re-annotate the gene column of an existing CNVkit reference .cnn.

Why: CNVkit tumor-only mode takes gene names from the reference .cnn's gene
column — `cnvkit batch --reference` rejects --annotate (and -n/-t/-a) as a
reference-construction flag (cnvkit >=0.9), so an un-annotated reference
(all genes "-") propagates the gap to every .cnr/.cns and report tables need
coordinate fallbacks. Re-annotating the .cnn itself fixes naming for all
future batch runs against it WITHOUT rebuilding the reference (the numeric
columns — log2, depth, gc, rmask, spread — are untouched, so per-bin
corrections and fingerprints stay valid; the gene column is pure metadata).

Usage:
    annotate_cnv_reference.py <reference.cnn> <genes.bed> [-o OUT.cnn]

The gene BED needs at least 4 columns (chrom, start, end, name), e.g.
GENCODE: bedtools intersect -a gencode.gtf -b /dev/null is NOT needed —
build one directly:
    gtf2bed < gencode.vXX.annotation.gtf | awk '$3=="gene"' \
        | cut -f1-4 > genes.bed
(or any gene-model BED4; bins overlapping multiple genes get comma-joined
names, matching `cnvkit.py target --annotate` behavior exactly — verified
bin-for-bin against cnvlib's into_ranges call in cnvlib/target.py).

Writes the annotated .cnn to stdout path (default: <input>.annotated.cnn);
the output keeps the tabular .cnn header+columns of the input.
"""
import argparse
import sys

from skgenome import tabio


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Add gene names to an existing CNVkit reference .cnn")
    ap.add_argument("reference", help="existing .cnn (tab format) reference")
    ap.add_argument("genes", help="gene model BED4 (chrom, start, end, name)")
    ap.add_argument("-o", "--output",
                    help="output .cnn (default: <reference>.annotated.cnn)")
    args = ap.parse_args(argv)

    out = args.output or (args.reference + ".annotated.cnn")

    ref = tabio.read(args.reference, "tab")
    ann = tabio.read(args.genes, "bed4")

    # Same call as cnvlib/target.py do_target(): any-overlap naming,
    # multi-gene bins comma-joined, non-overlapping bins keep the default.
    ref["gene"] = ann.into_ranges(ref, "gene", "-")
    # into_ranges may append the column out of position; restore the standard
    # .cnn layout (chromosome, start, end, gene, log2, then the rest).
    core = ["chromosome", "start", "end", "gene", "log2"]
    ordered = [c for c in core if c in ref.data.columns] + \
              [c for c in ref.data.columns if c not in core]
    ref.data = ref.data[ordered]
    tabio.write(ref, out)
    named = int((ref["gene"] != "-").sum())
    print(f"Wrote {out}: {len(ref)} bins, {named} gene-named", file=sys.stderr)


if __name__ == "__main__":
    main()
