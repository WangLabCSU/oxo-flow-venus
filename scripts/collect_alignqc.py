#!/usr/bin/env python3
"""Per-sample WGS alignment QC summary: flagstat + stats + mosdepth (N10/N11).

Consumes the four raw QC products of rules/qc.oxoflow (wgs_flagstat,
wgs_bamstats, wgs_mosdepth) and emits one row per DNA sample with mapping
rate, duplicate rate, insert size, GC, depth, coverage fractions, and the
sex-chromosome coverage check (N11). Input files are routed by SUFFIX (like
cohort_tables.py), so argument order does not matter:
  *.flagstat                 samtools flagstat text
  *.stats                    samtools stats text (SN fields + GC section)
  *.mosdepth.summary.txt     per-chromosome mean depth (+ 'total' row)
  *.mosdepth.thresholds.bed.gz / .bed
                             per-window bases >= 10x/20x (mosdepth --thresholds)

Sex check (N11, disclosed constants): compares mosdepth per-chromosome means
against the autosome average. diploid-normalised expectation is X/A ~= 1.0 (F)
or ~= 0.5 (M); Y/A ~= 0 (F) or >= 0.05 (M, mappability-reduced). This is a
screening call on COVERAGE ONLY — the cohort has no recorded sex metadata, so
every row is a hypothesis to review, not a sample-sheet correction.

Cross-sample contamination (N11): NOT ASSESSED. VerifyBamID2/freemix needs a
1000G genotype reference panel (SVDP bundle, GBs) and meaningful freemix needs
cohort density; deferring to P2 (see docs/20261104-qc-bundle.md). The
contamination column is an explicit not_assessed marker so downstream tables
never silently present a blank as a measurement.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

# --- sex-call thresholds (disclosed; tune only with cohort evidence) ---------
SEX_MALE_X_MAX = 0.75      # M: X/A clearly below diploid expectation
SEX_MALE_Y_MIN = 0.05      # M: Y/A above noise (cross-mapping background)
SEX_FEMALE_X_MIN = 0.85    # F: X/A at/near diploid expectation
SEX_FEMALE_Y_MAX = 0.02    # F: Y/A at cross-mapping background level

CONTAMINATION_MARK = "not_assessed(no_genotype_panel_p2)"

COLUMNS = ["sample", "total_reads", "mapped_pct", "duplicate_pct",
           "insert_size_median", "insert_size_average", "insert_size_stdev",
           "gc_pct", "mean_depth", "coverage_10x_pct", "coverage_20x_pct",
           "sex_x_ratio", "sex_y_ratio", "predicted_sex",
           "contamination_freemix"]

SUFFIXES = [".mosdepth.thresholds.bed.gz", ".mosdepth.thresholds.bed",
            ".mosdepth.summary.txt", ".flagstat", ".stats"]


def route(path: str) -> str:
    name = Path(path).name
    for suf in SUFFIXES:
        if name.endswith(suf):
            return suf
    raise SystemExit(f"[collect_alignqc] unrecognised QC file suffix: {name}")


def sample_of(path: str) -> str:
    name = Path(path).name
    for suf in SUFFIXES:
        if name.endswith(suf):
            return name[: -len(suf)]
    return name


def open_maybe_gzip(path: str):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_flagstat(path: str) -> dict[str, str]:
    total = mapped = dups = None
    mapped_pct = None
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^(\d+) \+ \d+ in total", line)
            if m:
                total = int(m.group(1))
                continue
            m = re.match(r"^(\d+) \+ \d+ duplicates", line)
            if m:
                dups = int(m.group(1))
                continue
            m = re.match(r"^(\d+) \+ \d+ mapped(?: \(([\d.]+)%)?", line)
            if m:
                mapped = int(m.group(1))
                mapped_pct = m.group(2)  # None on pre-1.14 flagstat (no pct)
    out = {}
    if total:
        out["total_reads"] = str(total)
        if mapped_pct is None and mapped is not None:
            mapped_pct = f"{mapped / total * 100:.4f}"
        if mapped_pct is not None:
            out["mapped_pct"] = mapped_pct
        if dups is not None:
            out["duplicate_pct"] = f"{dups / total * 100:.4f}"
    return out


def parse_stats(path: str) -> dict[str, str]:
    sn = {}
    gc_num = gc_den = 0.0
    with open(path) as fh:
        for line in fh:
            if line.startswith("SN\t"):
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3:
                    sn[f[1].rstrip(":")] = f[2]
            elif line.startswith("GC\t"):
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3:
                    gc = float(f[1])
                    n = int(f[2])
                    gc_num += gc * n
                    gc_den += n
    out = {}
    med = sn.get("insert size: median")
    avg = sn.get("insert size: average")
    std = sn.get("insert size: standard deviation")
    if med:
        out["insert_size_median"] = med
    if avg:
        out["insert_size_average"] = avg
    if std:
        out["insert_size_stdev"] = std
    if gc_den:
        out["gc_pct"] = f"{gc_num / gc_den:.2f}"
    return out


def parse_mosdepth_summary(path: str) -> dict[str, dict[str, float]]:
    """{chrom: {length, bases, mean}} rows ('total' included)."""
    rows = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[0] == "chrom":
                continue  # header
            try:
                rows[f[0]] = {"length": int(f[1]), "bases": int(f[2]),
                              "mean": float(f[3])}
            except ValueError:
                continue
    return rows


def autosome_mean(rows) -> float | None:
    """Coverage mean over chr1..chr22 (length-weighted), name-agnostic."""
    num = den = 0
    found = False
    for chrom, r in rows.items():
        c = chrom.removeprefix("chr")
        if c.isdigit() and 1 <= int(c) <= 22:
            num += r["bases"]
            den += r["length"]
            found = True
    return num / den if found and den else None


def parse_thresholds(path: str) -> dict[str, str]:
    """Fractions of reference bases covered >= 10x / 20x (mosdepth --thresholds
    over 1Mb windows tiling the reference; denominator = summed window length)."""
    idx10 = idx20 = None
    den = num10 = num20 = 0
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                cols = line.lstrip("#").split("\t")
                if "10X" in [c.upper() for c in cols]:
                    up = [c.upper() for c in cols]
                    idx10 = up.index("10X")
                    idx20 = up.index("20X") if "20X" in up else None
                continue
            f = line.split("\t")
            if len(f) < 5:
                continue
            start, end = int(f[1]), int(f[2])
            den += end - start
            if idx10 is not None and len(f) > idx10:
                num10 += int(f[idx10])
            if idx20 is not None and len(f) > idx20:
                num20 += int(f[idx20])
    out = {}
    if den:
        out["coverage_10x_pct"] = f"{num10 / den * 100:.4f}"
        if idx20 is not None:
            out["coverage_20x_pct"] = f"{num20 / den * 100:.4f}"
    return out


def call_sex(x_ratio, y_ratio):
    """Returns (predicted_sex, x_ratio_str, y_ratio_str); screening-level."""
    if x_ratio is None:
        return "ambiguous(no_chrX)", "", ""
    xs = f"{x_ratio:.3f}"
    if y_ratio is None:
        # no chrY contig in the reference/BAM: F still decidable, M is not
        if x_ratio >= SEX_FEMALE_X_MIN:
            return "F(no_chrY_contig)", xs, ""
        return "ambiguous(no_chrY)", xs, ""
    ys = f"{y_ratio:.3f}"
    if x_ratio < SEX_MALE_X_MAX and y_ratio >= SEX_MALE_Y_MIN:
        return "M", xs, ys
    if x_ratio >= SEX_FEMALE_X_MIN and y_ratio <= SEX_FEMALE_Y_MAX:
        return "F", xs, ys
    return "ambiguous", xs, ys


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="*", default=[],
                    help="flagstat/stats/mosdepth outputs, any order "
                         "(empty list -> header-only table)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    samples: dict[str, dict[str, str]] = {}
    summaries: dict[str, dict] = {}
    for path in args.inputs:
        kind = route(path)
        sample = sample_of(path)
        rec = samples.setdefault(sample, {"sample": sample,
                                          "contamination_freemix":
                                          CONTAMINATION_MARK})
        if kind == ".flagstat":
            rec.update(parse_flagstat(path))
        elif kind == ".stats":
            rec.update(parse_stats(path))
        elif kind == ".mosdepth.summary.txt":
            summaries[sample] = parse_mosdepth_summary(path)
        elif kind.startswith(".mosdepth.thresholds"):
            rec.update(parse_thresholds(path))

    # depth + sex after all summaries are in
    for sample, rec in samples.items():
        rows = summaries.get(sample)
        if not rows:
            continue
        total = rows.get("total")
        if total:
            rec["mean_depth"] = f"{total['mean']:.2f}"
        am = autosome_mean(rows)
        if am:
            # robust lookup regardless of 'chr' prefix
            xrow = rows.get("chrX") or rows.get("X")
            yrow = rows.get("chrY") or rows.get("Y")
            xr = xrow["mean"] / am if xrow else None
            yr = yrow["mean"] / am if yrow else None
            sex, xs, ys = call_sex(xr, yr)
            rec["predicted_sex"], rec["sex_x_ratio"], rec["sex_y_ratio"] = \
                sex, xs, ys

    ordered = [samples[s] for s in sorted(samples)]
    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(ordered)

    sexes = {}
    for rec in ordered:
        sexes.setdefault(rec.get("predicted_sex", ""), []).append(rec["sample"])
    print("[collect_alignqc] sex-call thresholds (coverage-based screening): "
          f"M: X/A<{SEX_MALE_X_MAX} and Y/A>={SEX_MALE_Y_MIN}; "
          f"F: X/A>={SEX_FEMALE_X_MIN} and Y/A<={SEX_FEMALE_Y_MAX}; "
          "else ambiguous", file=sys.stderr)
    for s, lst in sexes.items():
        print(f"[collect_alignqc] predicted_sex={s or 'n/a'}: "
              f"{len(lst)} samples", file=sys.stderr)
    print("[collect_alignqc] contamination: " + CONTAMINATION_MARK +
          " (VerifyBamID2 needs a 1000G SVDP genotype panel — P2)",
          file=sys.stderr)
    print(f"[collect_alignqc] {len(ordered)} samples -> {args.output}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
