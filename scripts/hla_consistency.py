#!/usr/bin/env python3
"""SpecHLA (DNA) vs arcasHLA (RNA) HLA typing consistency QC (venus#31 P1 item 3).

Joins the two typing layers of one patient and emits a per-gene/per-allele
consistency table plus an honest disclosure row for every one-sided absence:

  - SpecHLA side : hla/{experiment}.alleles.tsv (scripts/hla_parse.py 12-col
    schema; 2-field alleles in allele_1/allele_2, tumor-only caveat column)
  - arcasHLA side: hla_rna/{rna_sample}/{rna_sample}.genotype.json — the raw
    arcasHLA 0.6.0 output ({gene: [allele1, allele2]}, 3-field alleles).
    Loci below arcasHLA's min_count (75 reads) are silently OMITTED from the
    JSON (verified in scripts/genotype.py:748); the script renders that as
    "absent_rna", never as a mismatch.

Comparison happens at 2-field resolution (HLA-A*02:01) using the SAME
normalization as hla_parse.two_field — the G-group/P-group text arcas emits is
stripped the same way. Matches/mismatches are reported per gene with both
alleles considered as a set, so heterozygote orientation differences (arcas
swapping allele1/allele2) do not fake a mismatch.

Agreement is a QC signal only: DNA typing from tumor (organoid) DNA can show
HLA-LOH that RNA typing does not share (or vice versa); discordance is
disclosed, not auto-resolved. Either side missing renders class "absent_*",
an empty/malformed input renders "not_run"/"error" rows with a note.

Stdlib only.
"""
import argparse
import csv
import json
import os
import re
import sys

GENES = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
ALLELE_RE = re.compile(r"([A-Z0-9]+)\*(\d+):(\d+)")

HEADER = ["patient", "gene", "dna_alleles", "rna_alleles", "concordance",
          "detail"]


def two_field(allele):
    """Same normalization as scripts/hla_parse.py: two_field()."""
    s = str(allele or "").strip().replace("HLA_", "").replace("HLA-", "")
    m = ALLELE_RE.search(s)
    return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}" if m else ""


def load_dna(path, sample):
    """SpecHLA alleles.tsv -> {gene: set(2-field alleles)}; '' gene = not typed."""
    out = {}
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r.get("sample") == sample or sample is None]
    # the table is per-sample already; fall back to all rows if the sample
    # label differs (single-sample table written by hla_parse.py)
    if not rows:
        with open(path) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        g = (r.get("gene") or "").strip()
        if g not in GENES:
            continue
        alleles = {two_field(r.get("allele_1")), two_field(r.get("allele_2"))}
        out[g] = {a for a in alleles if a}
    return out


def load_rna(path):
    """arcasHLA genotype.json -> {gene: set(2-field alleles)}.

    Genes omitted by arcasHLA (low coverage, < min_count) are simply absent
    from the JSON — the caller renders that as absent_rna.
    """
    with open(path) as fh:
        raw = json.load(fh)
    if not isinstance(raw, dict):
        raise ValueError(f"{path}: expected {{gene: [a1, a2]}} JSON object")
    out = {}
    for g, alleles in raw.items():
        gene = str(g).replace("HLA-", "").strip()
        if gene not in GENES:
            continue
        if isinstance(alleles, str):
            alleles = [alleles]
        out[gene] = {two_field(a) for a in alleles if two_field(a)}
    return out


def gene_row(patient, gene, dna_set, rna_set):
    """One gene -> one consistency row (honest class labels)."""
    dna = dna_set.get(gene) or set()
    rna = rna_set.get(gene) or set()
    if not dna and not rna:
        return {"patient": patient, "gene": gene, "dna_alleles": "",
                "rna_alleles": "", "concordance": "not_run",
                "detail": "neither layer typed this locus"}
    if not dna:
        return {"patient": patient, "gene": gene,
                "dna_alleles": "", "rna_alleles": ";".join(sorted(rna)),
                "concordance": "absent_dna",
                "detail": "SpecHLA did not type this locus (tumor DNA typing "
                          "gap / LOH); RNA-side call shown"}
    if not rna:
        return {"patient": patient, "gene": gene,
                "dna_alleles": ";".join(sorted(dna)), "rna_alleles": "",
                "concordance": "absent_rna",
                "detail": "arcasHLA omitted this locus (below min_count 75 "
                          "extracted reads, or not requested); DNA-side call "
                          "shown"}
    if dna == rna:
        return {"patient": patient, "gene": gene,
                "dna_alleles": ";".join(sorted(dna)),
                "rna_alleles": ";".join(sorted(rna)),
                "concordance": "match",
                "detail": "2-field alleles agree (order-insensitive)"}
    shared = dna & rna
    detail = "2-field allele sets differ"
    if shared:
        detail += (f"; shared: {';'.join(sorted(shared))}; "
                   f"dna-only: {';'.join(sorted(dna - rna))}; "
                   f"rna-only: {';'.join(sorted(rna - dna))}")
    return {"patient": patient, "gene": gene,
            "dna_alleles": ";".join(sorted(dna)),
            "rna_alleles": ";".join(sorted(rna)),
            "concordance": "mismatch", "detail": detail}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dna-tsv", required=True,
                    help="SpecHLA alleles.tsv (hla_parse.py schema)")
    ap.add_argument("--dna-sample", required=True,
                    help="DNA sample label (pairs.tsv experiment, e.g. Pt01-TD)")
    ap.add_argument("--rna-json", required=True,
                    help="arcasHLA {sample}.genotype.json")
    ap.add_argument("--rna-sample", required=True,
                    help="RNA sample label (e.g. Pt01-TR)")
    ap.add_argument("--patient", required=True, help="pair_id")
    ap.add_argument("--out", required=True, help="consistency TSV (8 rows)")
    args = ap.parse_args()

    if not os.path.exists(args.dna_tsv):
        sys.exit(f"[hla_consistency] DNA table missing: {args.dna_tsv}")
    if not os.path.exists(args.rna_json):
        sys.exit(f"[hla_consistency] arcasHLA JSON missing: {args.rna_json}")
    try:
        dna = load_dna(args.dna_tsv, args.dna_sample)
        rna = load_rna(args.rna_json)
    except (ValueError, KeyError) as exc:
        sys.exit(f"[hla_consistency] unreadable input: {exc}")

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for g in GENES:
            w.writerow(gene_row(args.patient, g, dna, rna))

    n_match = sum(1 for g in GENES
                  if (dna.get(g) or set()) and (rna.get(g) or set())
                  and dna[g] == rna[g])
    n_cmp = sum(1 for g in GENES
                if (dna.get(g) or set()) and (rna.get(g) or set()))
    sys.stderr.write(
        f"[hla_consistency] {args.patient}: {n_match}/{n_cmp} "
        f"comparable loci match at 2-field resolution "
        f"(DNA={args.dna_sample}, RNA={args.rna_sample})\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
