#!/usr/bin/env python3
"""SpecHLA built-in HLA-LOH caller driver (venus#31 P1 item 4).

DASH (the publication-standard HLA-LOH caller, Shukla et al. 2015) is NOT
feasible under venus's pixi lock-only env policy: it has no bioconda/pypi
package (verified 2026-11 against both channels) and would need a
hand-built toolchain outside the lock. The honest alternative shipped INSIDE
the existing SpecHLA 1.0.12 env (envs/hla) is `spechla-loh`
(share/spechla/script/cal.hla.copy.pl — verified from the actual bioconda
package, NOT from memory):

  spechla-loh -S <sample> -purity <f> -ploidy <f> -C <het_cutoff=5>
              -F <filelist of *_freq.txt> -T <hla.result.txt> -O <outdir>
  -> <outdir>/merge.hla.copy.txt with columns
     Sample HLA Allele1 Allele2 copyratio KeptHLA LossHLA Freq1 Freq2
     Purity Het_num LOH
  LOH = "Y" iff the inferred allele copy ratios differ AND the heterozygous
  SNP count >= -C (default 5). One-allele-zero-frequency ("homogeneous")
  loci are tagged LOH = "N" by the tool itself — copy-neutral loss of one
  allele's phasing is NOT LOH evidence and is never re-interpreted here.

This driver is thin glue only:
  1. purity/ploidy come from the ASCAT fit (ascat/{pair}/purity_ploidy.tsv,
     run_ascat.R schema: pair, sample, purity, ploidy, psi, ...) — the
     cal.hla.copy.pl usage text itself says tumor purity and ploidy "can be
     inferred by software like ABSOLUTE and ASCAT". Missing/failed ASCAT fit
     fails LOUDLY (no silent default purity).
  2. writes the freq filelist (the 8 per-gene {sample}/HLA_<gene>_freq.txt
     files produced by SpecHLA.sh's phase_variants.py and NOT deleted by
     clear_output.sh — verified), then invokes spechla-loh.
  3. passes merge.hla.copy.txt through verbatim + pair/note columns.

Stdlib only; spechla-loh must be on PATH (envs/hla pixi env).
"""
import argparse
import csv
import os
import subprocess
import sys

MERGE_HEADER = ["Sample", "HLA", "Allele1", "Allele2", "copyratio", "KeptHLA",
                "LossHLA", "Freq1", "Freq2", "Purity", "Het_num", "LOH"]
OUT_HEADER = ["pair"] + MERGE_HEADER + ["note"]
FREQ_GENES = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


def load_purity_ploidy(path, pair_id, sample):
    """run_ascat.R purity_ploidy.tsv -> (purity, ploidy); loud on failure."""
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")]
    if not rows:
        sys.exit(f"[hla_loh] {path}: no rows — ASCAT fit missing/failed; "
                 "spechla-loh REQUIRES purity and ploidy, refusing to guess")
    row = next((r for r in rows
                if r.get("pair") == pair_id or r.get("sample") == sample),
               None)
    if row is None:
        sys.exit(f"[hla_loh] {path}: no row for pair={pair_id} "
                 f"sample={sample}")
    try:
        purity = float(row["purity"])
        ploidy = float(row["ploidy"])
    except (TypeError, ValueError, KeyError) as exc:
        sys.exit(f"[hla_loh] {path}: unusable purity/ploidy ({exc}) — ASCAT "
                 "fit failed for this pair; not fabricating defaults")
    if not (0.0 < purity <= 1.0) or ploidy <= 0:
        sys.exit(f"[hla_loh] {path}: out-of-range purity={purity} "
                 f"ploidy={ploidy} — refusing to run LOH on a broken fit")
    return purity, ploidy


def write_freq_list(path, freq_files):
    with open(path, "w") as fh:
        for f in freq_files:
            fh.write(f"{f}\n")


def parse_merge(path, pair_id):
    """merge.hla.copy.txt -> list of row dicts (verbatim pass-through).

    Adds a note column only for the tool's own "homogeneous" tag (LossHLA
    column == 'homogeneous', LOH = 'N'): one allele had zero haplotype
    frequency — the tool calls this LOH=N and we do not second-guess it.
    """
    with open(path) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r]
    if not rows:
        sys.exit(f"[hla_loh] {path}: empty merge.hla.copy.txt")
    if rows[0] != MERGE_HEADER:
        sys.exit(f"[hla_loh] {path}: unexpected header {rows[0]} — "
                 f"expected {MERGE_HEADER}")
    out = []
    for cols in rows[1:]:
        if len(cols) != len(MERGE_HEADER):
            sys.exit(f"[hla_loh] {path}: row with {len(cols)} cols "
                     f"(expected {len(MERGE_HEADER)}): {cols}")
        rec = dict(zip(MERGE_HEADER, cols))
        rec["pair"] = pair_id
        rec["note"] = ("tool-tagged homogeneous (one allele zero-frequency; "
                       "LOH=N by definition)" if rec["LossHLA"].startswith(
                           "homogeneous") else "")
        out.append(rec)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pair", required=True, help="pair_id")
    ap.add_argument("--sample", required=True,
                    help="tumor DNA sample label (pairs.tsv experiment)")
    ap.add_argument("--result", required=True,
                    help="SpecHLA hla.result.txt (allele-name source for "
                         "spechla-loh)")
    ap.add_argument("--freq", nargs="+", required=True,
                    help="the 8 HLA_<gene>_freq.txt files (exact, declared)")
    ap.add_argument("--purity-tsv", required=True,
                    help="ASCAT ascat/{pair}/purity_ploidy.tsv")
    ap.add_argument("--het-cutoff", type=float, default=5.0,
                    help="spechla-loh -C heterogeneous-SNP cutoff (default 5)")
    ap.add_argument("--outdir", required=True,
                    help="working dir for spechla-loh (-O; holds "
                         "merge.hla.copy.txt)")
    ap.add_argument("--out", required=True, help="tidy LOH TSV")
    args = ap.parse_args()

    purity, ploidy = load_purity_ploidy(args.purity_tsv, args.pair,
                                        args.sample)
    for f in args.freq:
        if not os.path.exists(f):
            sys.exit(f"[hla_loh] missing freq file: {f} — SpecHLA did not "
                     "phase this gene; declare the failure, do not skip it")
    os.makedirs(args.outdir, exist_ok=True)

    freq_list = os.path.join(args.outdir, "freq.list")
    write_freq_list(freq_list, args.freq)

    cmd = ["spechla-loh",
           "-S", args.sample,
           "-purity", f"{purity:.6f}",
           "-ploidy", f"{ploidy:.6f}",
           "-C", f"{args.het_cutoff:g}",
           "-F", freq_list,
           "-T", args.result,
           "-O", args.outdir]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"[hla_loh] spechla-loh failed rc={proc.returncode}\n"
                 f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")

    merge = os.path.join(args.outdir, "merge.hla.copy.txt")
    rows = parse_merge(merge, args.pair)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_HEADER, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    n_loh = sum(1 for r in rows if r["LOH"] == "Y")
    sys.stderr.write(
        f"[hla_loh] {args.pair}: {n_loh}/{len(rows)} HLA loci LOH=Y "
        f"(purity={purity:.3f}, ploidy={ploidy:.3f}, het_cutoff="
        f"{args.het_cutoff:g}; caller=spechla-loh/SpecHLA 1.0.12; DASH "
        f"not packaged on bioconda/pypi — disclosed in method rationale)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
