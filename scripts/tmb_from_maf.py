#!/usr/bin/env python3
"""Per-sample TMB from one or more MAF files, with two clearly separated
denominators:

- tmb_mut_per_mb = coding non-silent mutations / coding Mb. This is the
  clinically comparable TMB (the >=10 mut/Mb high-TMB cutoff is defined on
  coding territory). --coding-mb defaults to 34.0 Mb (approx. CCDS coding
  length), which is a reasonable denominator for BOTH WGS and WES; for WES
  you can pass --target-bed and the coding territory inside the capture is
  used instead (only when the BED total is smaller).
- all_somatic_mut_per_mb = all PASS mutations / target Mb. WGS research
  convention. --target-mb defaults to 3000 (approx. callable WGS); with
  --target-bed the Mb is computed from the merged BED intervals, which is
  the honest denominator for WES.

Counting rules: 'count' = unique (Chromosome, Start_Position, End_Position,
Reference_Allele, Tumor_Seq_Allele2); 'coding' keeps only coding classes
(non-silent) after dropping 'Precision'-qualified rows. # comment lines and
version banners are skipped.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict

CODING_CLASSES = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Ins",
    "In_Frame_Del",
    "Splice_Site",
    "Translation_Start_Site",
    "Nonstop_Mutation",
}

COLUMNS = [
    "sample",
    "total_pass_variants",
    "coding_non_silent",
    "snps",
    "indels",
    "coding_mb",
    "target_mb",
    "tmb_mut_per_mb",
    "all_somatic_mut_per_mb",
]


def merge_bed(bed_path: str) -> tuple[int, int]:
    """Sum merged BED interval bases.

    Returns (merged_bases, original_bases). Overlapping/touching intervals
    are merged per chromosome so the reported Mb is not inflated.
    """
    per_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    total_bases = 0
    opener = gzip.open if bed_path.endswith(".gz") else open
    with opener(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            start, end = int(f[1]), int(f[2])
            if end > start:
                per_chrom[f[0]].append((start, end))
                total_bases += end - start
    merged_bases = 0
    for ivs in per_chrom.values():
        ivs.sort()
        cur_s, cur_e = ivs[0]
        for s, e in ivs[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged_bases += cur_e - cur_s
                cur_s, cur_e = s, e
        merged_bases += cur_e - cur_s
    return merged_bases, total_bases


def read_maf(path: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    if not lines:
        return rows
    header = lines[0].lstrip("#").rstrip("\n").split("\t")
    for ln in lines[1:]:
        f = ln.rstrip("\n").split("\t")
        if len(f) != len(header):
            print(f"[tmb] skipping malformed row in {path}", file=sys.stderr)
            continue
        rows.append(dict(zip(header, f)))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mafs", nargs="+", required=True)
    ap.add_argument(
        "--target-bed",
        default="",
        help="capture/target interval BED; target Mb computed from merged intervals",
    )
    ap.add_argument(
        "--target-mb",
        type=float,
        default=3000.0,
        help="target territory Mb when no BED given (WGS default ~3000)",
    )
    ap.add_argument(
        "--coding-mb",
        type=float,
        default=34.0,
        help="coding territory Mb for the clinical TMB (default ~34 Mb CCDS)",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    coding_mb = args.coding_mb
    target_mb = args.target_mb
    if args.target_bed:
        merged, raw = merge_bed(args.target_bed)
        bed_mb = merged / 1e6
        if raw != merged:
            print(
                f"[tmb] target BED: {raw} raw bases merged to {merged} "
                f"({bed_mb:.2f} Mb)",
                file=sys.stderr,
            )
        else:
            print(f"[tmb] target BED: {merged} bases ({bed_mb:.2f} Mb)", file=sys.stderr)
        target_mb = bed_mb
        if bed_mb < coding_mb:
            coding_mb = bed_mb  # WES: coding territory cannot exceed captured

    per_sample: dict[str, dict] = defaultdict(
        lambda: {
            "total": 0,
            "coding": 0,
            "snps": 0,
            "indels": 0,
            "keys": set(),
        }
    )
    for path in args.mafs:
        for row in read_maf(path):
            key = (
                row.get("Chromosome", ""),
                row.get("Start_Position", ""),
                row.get("End_Position", ""),
                row.get("Reference_Allele", ""),
                row.get("Tumor_Seq_Allele2", ""),
            )
            s = per_sample[row.get("Tumor_Sample_Barcode", "UNKNOWN")]
            if key in s["keys"]:
                continue
            s["keys"].add(key)
            s["total"] += 1
            vtype = row.get("Variant_Type", "")
            if vtype == "SNP":
                s["snps"] += 1
            elif vtype in ("INS", "DEL"):
                s["indels"] += 1
            vc = row.get("Variant_Classification", "")
            if vc in CODING_CLASSES:
                s["coding"] += 1

    with open(args.output, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for sample in sorted(per_sample):
            s = per_sample[sample]
            tmb = s["coding"] / coding_mb if coding_mb else 0.0
            tmb_all = s["total"] / target_mb if target_mb else 0.0
            w.writerow([
                sample,
                s["total"],
                s["coding"],
                s["snps"],
                s["indels"],
                f"{coding_mb:.2f}",
                f"{target_mb:.2f}",
                f"{tmb:.3f}",
                f"{tmb_all:.3f}",
            ])
    print(f"[tmb] wrote {len(per_sample)} samples -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
