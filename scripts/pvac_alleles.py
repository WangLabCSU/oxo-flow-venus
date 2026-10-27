#!/usr/bin/env python3
"""Build the pVACseq allele list: SpecHLA alleles ∩ MHCflurry-supported set.

- Reads the hla_parse.py long TSV (2-field "HLA-A*02:01" alleles).
- Intersects with the predictor's supported-allele list so pvacseq never
  dies on a rare typed allele (venus#31 P0 runs MHCflurry only).
- Refuses to emit fewer than --min-alleles typed+supported alleles: an
  HLA typing that collapses to nothing is a QC failure, not a runnable
  input (fail loud, keep the report honest).

Stdlib only (supported-allele discovery is via the pvac/mhcflurry env CLI).
"""
import argparse
import csv
import subprocess
import sys


def mhcflurry_supported():
    try:
        out = subprocess.run(
            ["mhcflurry-predict", "--supported-alleles"],
            capture_output=True, text=True, check=True, timeout=300)
        return {ln.strip() for ln in out.stdout.splitlines() if ln.strip()}
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) \
            as exc:
        sys.exit(f"[pvac_alleles] cannot enumerate MHCflurry supported "
                 f"alleles: {exc}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alleles-tsv", required=True)
    ap.add_argument("--min-alleles", type=int, default=2)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    typed = set()
    with open(args.alleles_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in ("allele_1", "allele_2"):
                if r.get(k, "").strip():
                    typed.add(r[k].strip())
    supported = mhcflurry_supported()
    usable = sorted(typed & supported)
    unsupported = sorted(typed - supported)
    if unsupported:
        sys.stderr.write(f"[pvac_alleles] typed but not supported by "
                         f"MHCflurry ({len(unsupported)}): "
                         f"{', '.join(unsupported)}\n")
    if len(usable) < args.min_alleles:
        sys.exit(f"[pvac_alleles] only {len(usable)} typed allele(s) usable "
                 f"(< --min-alleles {args.min_alleles}) — treat as HLA-typing "
                 f"QC failure, not an empty result")
    with open(args.out, "w") as fh:
        fh.write(",".join(usable) + "\n")
    sys.stderr.write(f"[pvac_alleles] {len(usable)} usable alleles "
                     f"-> {args.out}\n")


if __name__ == "__main__":
    main()
