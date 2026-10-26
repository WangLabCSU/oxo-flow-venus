#!/usr/bin/env python3
"""Fusion-candidate detection from Manta SV VCF for venus clinical reporting.

Maps both breakpoints of each BND record onto gene spans from the GENCODE GTF
and reports candidates matching the curated actionable fusion pair table
(resources/clinical/actionable_fusion_pairs.tsv). A `*` gene column means any
partner matches (single-gene rule, e.g. FGFR2 promoter swaps).

Tumor-only Manta output is candidateSV (unscored); every row is labeled
accordingly. Output: report/clinical/{pair_id}.sv_evidence.tsv
"""
import argparse
import re
import sys

import pandas as pd
import pysam

BND_RE = re.compile(r"[\[\]]([A-Za-z0-9]+):(\d+)[\[\]]")


def gene_spans_from_gtf(gtf_path, wanted):
    wanted_u = {w.upper() for w in wanted if w != "*"}
    spans = {}
    with open(gtf_path) as fh:
        for line in fh:
            if "\tgene\t" not in line or 'gene_name "' not in line:
                continue
            name = line.split('gene_name "')[1].split('"')[0]
            if name.upper() not in wanted_u:
                continue
            fl = line.split("\t")
            spans[name] = (fl[0], int(fl[3]), int(fl[4]))
    return spans


def gene_at(spans, chrom, pos):
    for g, (c, s, e) in spans.items():
        if c == chrom and s <= pos <= e:
            return g
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="Manta candidateSV/somaticSV VCF")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--pairs", required=True, help="actionable_fusion_pairs.tsv")
    ap.add_argument("--tumor-only", action="store_true",
                    help="label rows as unscored candidateSV")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t", dtype=str).fillna("")
    genes = set(pairs["gene_5"]) | set(pairs["gene_3"])
    spans = gene_spans_from_gtf(args.gtf, genes)
    sys.stderr.write(f"[clinical_sv] pair-table genes with GTF spans: "
                     f"{len(spans)}/{len(genes)}\n")

    vcf = pysam.VariantFile(args.vcf)
    # Manta emits each BND twice (record + MATEID mirror). Mate positions are
    # shifted by HOMLEN/CIPOS, so exact coordinate pairs cannot dedupe them —
    # keep the FIRST record of each MATEID pair (second occurrence's mate was
    # already emitted). Coordinate-pair dedupe stays as fallback for callers
    # that emit no MATEID.
    seen_ids = set()
    seen_bnd = set()
    rows = []
    for rec in vcf.fetch():
        if (rec.info.get("SVTYPE") or "") != "BND":
            continue
        if not args.tumor_only:
            filters = set(rec.filter.keys()) if rec.filter else set()
            if filters and not filters <= {"PASS"}:
                continue  # scored somaticSV: honour the caller's filters
        mate = mate[0] if isinstance(mate, (list, tuple)) else mate
        if mate and mate in seen_ids:
            continue  # mirror of an already-kept record
        m = BND_RE.search(rec.alts[0] or "")
        if not m:
            continue
        g1 = gene_at(spans, rec.chrom, rec.pos)
        g2 = gene_at(spans, m.group(1), int(m.group(2)))
        if not g1 and not g2:
            continue
        key = tuple(sorted([f"{rec.chrom}:{rec.pos}", f"{m.group(1)}:{m.group(2)}"]))
        if key in seen_bnd:
            continue
        seen_bnd.add(key)
        if rec.id:
            seen_ids.add(rec.id)
        for _, p in pairs.iterrows():
            a, b = p["gene_5"], p["gene_3"]
            match = (g1 and g2 and {g1.upper(), g2.upper()} == {a.upper(), b.upper()}) \
                or (b == "*" and g1 and g1.upper() == a.upper()) \
                or (b == "*" and g2 and g2.upper() == a.upper()) \
                or (a == "*" and g2 and g2.upper() == b.upper())
            if not match:
                continue
            if g1 and g2:
                fusion = f"{g1}--{g2}"
            else:
                ghit = g1 or g2
                fusion = (f"{ghit} (single breakpoint; partner unresolved)")
            rows.append({
                "fusion_candidate": fusion,
                "gene_5_call": g1, "gene_3_call": g2,
                "breakpoint1": f"{rec.chrom}:{rec.pos}",
                "breakpoint2": f"{m.group(1)}:{m.group(2)}",
                "sv_id": rec.id or "",
                "filter": ";".join(rec.filter.keys()) if rec.filter else "",
                "scoring": "candidateSV (unscored in tumor-only)" if args.tumor_only
                           else "somaticSV",
                "disease_context": p["disease_context"],
                "evidence_context": p["evidence_context"],
                "notes": p.get("notes", ""),
            })

    out = pd.DataFrame(rows, columns=[
        "fusion_candidate", "gene_5_call", "gene_3_call", "breakpoint1",
        "breakpoint2", "sv_id", "filter", "scoring", "disease_context",
        "evidence_context", "notes"])
    out.to_csv(args.out, sep="\t", index=False)
    sys.stderr.write(f"[clinical_sv] wrote {len(out)} fusion candidates -> {args.out}\n")


if __name__ == "__main__":
    main()
