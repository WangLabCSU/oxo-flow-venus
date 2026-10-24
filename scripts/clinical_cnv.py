#!/usr/bin/env python3
"""Gene-level CNV evidence extraction from CNVkit segments for venus.

Reads the pair's call.cns segments, maps a curated actionability gene panel
(resources/clinical/actionable_cnv_genes.tsv) onto gene spans from the GENCODE
GTF (config.gencode_gtf), and classifies focal vs broad amplification/deletion
with conservative thresholds:
  amp:  segment log2 >= --amp-log2 (0.3) and cn >= --amp-cn (4)
        focal when the overlapping segment is <= --focal-mb (5 Mb)
  del:  segment log2 <= --del-log2 (-0.5)
Cross-flags ecDNA membership from the AmpliconSuite gene list when provided.

Output: report/clinical/{pair_id}.cnv_evidence.tsv
"""
import argparse
import sys

import pandas as pd

CANON = {"chr%d" % i for i in range(1, 23)} | {"chrX", "chrY"}


def gene_spans_from_gtf(gtf_path, wanted):
    """gene -> (chrom, start, end) for the wanted gene symbols (tabix-free scan)."""
    wanted_u = {w.upper() for w in wanted}
    spans = {}
    with open(gtf_path) as fh:
        for line in fh:
            if "\tgene\t" not in line:
                continue
            if 'gene_name "' not in line:
                continue
            name = line.split('gene_name "')[1].split('"')[0]
            if name.upper() not in wanted_u:
                continue
            fl = line.split("\t")
            spans[name] = (fl[0], int(fl[3]), int(fl[4]))
    return spans


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cns", required=True, help="CNVkit call.cns")
    ap.add_argument("--gtf", required=True, help="GENCODE GTF (gene records)")
    ap.add_argument("--panel", required=True, help="actionable_cnv_genes.tsv")
    ap.add_argument("--ecdna-genes", default="",
                    help="optional AmpliconSuite gene list TSV for ecDNA flags")
    ap.add_argument("--amp-log2", type=float, default=0.3)
    ap.add_argument("--amp-cn", type=int, default=4)
    ap.add_argument("--del-log2", type=float, default=-0.5)
    ap.add_argument("--focal-mb", type=float, default=5.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    panel = pd.read_csv(args.panel, sep="\t", dtype=str).fillna("")
    spans = gene_spans_from_gtf(args.gtf, set(panel["gene"]))
    sys.stderr.write(f"[clinical_cnv] panel genes with GTF spans: "
                     f"{len(spans)}/{len(panel)}\n")

    ecgenes = set()
    if args.ecdna_genes:
        try:
            ed = pd.read_csv(args.ecdna_genes, sep="\t", dtype=str)
            col = "gene" if "gene" in ed.columns else ed.columns[0]
            ecgenes = {str(g).upper() for g in ed[col].dropna()}
        except Exception as exc:  # tolerate malformed/empty archive
            sys.stderr.write(f"[clinical_cnv] ecdna gene list unreadable: {exc}\n")

    seg = pd.read_csv(args.cns, sep="\t")
    seg = seg[seg["chromosome"].isin(CANON)]
    rows = []
    for _, p in panel.iterrows():
        gene = p["gene"]
        if gene not in spans:
            continue
        chrom, gs, ge = spans[gene]
        hit = seg[(seg.chromosome == chrom) & (seg.end > gs) & (seg.start < ge)]
        if hit.empty:
            continue
        # strongest deviation across overlapping segments
        best_amp = hit.loc[hit.log2.idxmax()] if len(hit) else None
        best_del = hit.loc[hit.log2.idxmin()] if len(hit) else None
        call, log2, cn, segdesc = "", "", "", ""
        if best_amp is not None and best_amp.log2 >= args.amp_log2 \
                and ("cn" not in hit.columns or int(best_amp.get("cn", 0)) >= args.amp_cn):
            length_mb = (best_amp.end - best_amp.start) / 1e6
            call = "focal_amplification" if length_mb <= args.focal_mb else "amplification"
            log2, cn = round(float(best_amp.log2), 3), int(best_amp.get("cn", 0))
            segdesc = f"{best_amp.chromosome}:{best_amp.start}-{best_amp.end}({length_mb:.2f}Mb)"
        if best_del is not None and best_del.log2 <= args.del_log2:
            length_mb = (best_del.end - best_del.start) / 1e6
            call = call or ("focal_deletion" if length_mb <= args.focal_mb else "deletion")
            log2, cn = round(float(best_del.log2), 3), int(best_del.get("cn", 0))
            segdesc = f"{best_del.chromosome}:{best_del.start}-{best_del.end}({length_mb:.2f}Mb)"
        if not call:
            continue
        rows.append({
            "gene": gene, "call": call, "segment_log2": log2, "segment_cn": cn,
            "segment": segdesc, "ecDNA_evidence": "yes" if gene.upper() in ecgenes else "",
            "disease_context": p["disease_context"], "therapy_context": p["therapy_context"],
            "notes": p.get("notes", ""),
        })

    out = pd.DataFrame(rows, columns=[
        "gene", "call", "segment_log2", "segment_cn", "segment",
        "ecDNA_evidence", "disease_context", "therapy_context", "notes"])
    out.to_csv(args.out, sep="\t", index=False)
    sys.stderr.write(f"[clinical_cnv] wrote {len(out)} calls -> {args.out}\n")


if __name__ == "__main__":
    main()
