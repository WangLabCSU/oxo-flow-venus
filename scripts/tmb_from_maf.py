#!/usr/bin/env python3
"""Compute TMB and per-sample variant summaries from PASS-filtered MAFs.

TMB (mutations per Mb) = count of non-silent coding somatic variants divided
by the callable territory: target_mb config for WES/panel, or ~3000 Mb
default for WGS (overridable). Non-silent classes follow the standard TCGA
definition (missense/nonsense/frameshift/inframe/splice/start-lost/nonstop).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

NON_SILENT = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
    "Frame_Shift_Ins", "In_Frame_Ins", "In_Frame_Del",
    "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
}


def read_maf(path: str) -> list[dict[str, str]]:
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    if not lines:
        return []
    header = lines[0].lstrip("#").rstrip("\n").split("\t")
    rows = []
    for line in lines[1:]:
        f = line.rstrip("\n").split("\t")
        rows.append(dict(zip(header, f)))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mafs", nargs="+", required=True)
    ap.add_argument("--target-mb", type=float, default=3000.0)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    out_rows = []
    for path in args.mafs:
        mafs = read_maf(path)
        sample = ""
        for r in mafs:
            sample = r.get("Tumor_Sample_Barcode", "")
            if sample:
                break
        if not sample:
            sample = Path(path).name.replace(".somatic.pass.maf", "").replace(".maf", "")
        total = len(mafs)
        coding = sum(1 for r in mafs if r.get("Variant_Classification") in NON_SILENT)
        snps = sum(1 for r in mafs if r.get("Variant_Type") == "SNP")
        indels = sum(1 for r in mafs if r.get("Variant_Type") in ("INS", "DEL"))
        tmb = (coding / args.target_mb) if args.target_mb > 0 else 0.0
        out_rows.append({
            "sample": sample,
            "total_pass_variants": total,
            "coding_non_silent": coding,
            "snps": snps,
            "indels": indels,
            "target_mb": args.target_mb,
            "tmb_mut_per_mb": f"{tmb:.3f}",
        })

    columns = ["sample", "total_pass_variants", "coding_non_silent", "snps",
               "indels", "target_mb", "tmb_mut_per_mb"]
    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"[tmb_from_maf] {len(out_rows)} samples -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
