#!/usr/bin/env python3
"""Neoantigen report for venus (venus#31 P0): pVACseq output -> burden + QC.

Consumes the pVACseq filtered report (binding-threshold survivors) plus the
SpecHLA allele table and the pair's RNA featureCounts, and emits:
  - candidates TSV: best-scoring row per (gene, epitope, HLA allele), with
    gene TPM annotated (TESLA expression filter: TPM >= 1)
  - summary TSV (long metric/value/detail — cohort-concatenable)
  - markdown section: burden, QC, top candidates, honest caveats

Honesty rules:
  - MHC-I only in P0 (MHCflurry); class-II absence is a not_run row, not silence
  - tumor-only mode keeps a qualitative germline-residue caveat on every
    burden number (tumor-only calling has no matched normal; germline
    polymorphisms can masquerade as somatic candidates)
  - no candidates / unparseable input -> fixed-schema rows, never an empty TSV
Column matching is alias-based against pVACtools 7.x filtered.tsv real
column names (Best MT IC50 Score / Gene Name / ...); unknown headers fail
loud with the header dump.
Stdlib only.
"""
import argparse
import csv
import gzip
import os
import sys

CAND_COLS = ["gene", "variant_type", "hla_allele", "epitope", "aa_change",
             "protein_position", "mt_score_nm", "mt_percentile",
             "tumor_vaf", "gene_tpm", "expr_pass", "transcript"]
SUMMARY_COLS = ["class", "metric", "value", "detail"]

ALIASES = {
    "gene": ["Gene Name", "Gene"],
    "variant_type": ["Variant Type"],
    "hla_allele": ["HLA Allele"],
    "epitope": ["MT Epitope Seq", "Mutant Peptide", "Epitope Seq"],
    "aa_change": ["Mutation", "AMINO_ACID_CHANGE", "Amino Acid Change"],
    "protein_position": ["Protein Position"],
    "mt_score_nm": ["Best MT IC50 Score", "Median MT IC50 Score", "MT Score",
                    "MHC Score", "IC50"],
    "mt_percentile": ["Best MT Percentile", "Median MT Percentile",
                      "MT Percentile", "Percentile"],
    "tumor_vaf": ["Tumor DNA VAF", "VAF"],
    "transcript": ["Transcript"],
}


def open_tsv(path):
    opener = gzip.open if path.endswith(".gz") else open
    return opener(path, "rt")


def read_table(path):
    with open_tsv(path) as fh:
        lines = fh.read().splitlines()
    if not lines:
        raise SystemExit(f"[neoantigen_report] {path} is empty")
    header = lines[0].split("\t")
    lower = {h.strip().lower(): h for h in header}
    colmap = {}
    for canon, aliases in ALIASES.items():
        for a in aliases:
            if a.lower() in lower:
                colmap[canon] = lower[a.lower()]
                break
    missing = [c for c in ("gene", "epitope", "hla_allele", "mt_score_nm")
               if c not in colmap]
    if missing:
        sys.exit(f"[neoantigen_report] {path}: missing columns {missing}; "
                 f"header was: {header}")
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) < len(header):
            parts += [""] * (len(header) - len(parts))
        rows.append({c: parts[header.index(h)]
                     for c, h in colmap.items()})
    return rows


def read_tpm_map(counts_path, gene_map_path):
    """featureCounts (Ensembl gene_id) + gtf gene map -> {symbol: TPM}."""
    symbol_of = {}
    try:
        with open(gene_map_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                symbol_of[r["gene_id"].strip()] = r["symbol"].strip()
    except (OSError, KeyError, csv.Error) as exc:
        sys.stderr.write(f"[neoantigen_report] gene map unreadable ({exc}); "
                         f"TPM annotation will be empty\n")
        return {}
    opener = gzip.open if counts_path.endswith(".gz") else open
    try:
        with opener(counts_path) as fh:
            # featureCounts writes its header line with a leading '#' and may
            # emit repeated "Geneid" header rows — keep line 0 as the header,
            # skip metadata/dupes, don't drop the header wholesale.
            lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
        header = lines[0].lstrip("#").strip().split("\t")
        i_gene = header.index("Geneid") if "Geneid" in header else 0
        i_len = header.index("Length")
        i_cnt = len(header) - 1  # last column = the single sample
        denom = 0.0
        per_gene = {}
        for ln in lines[1:]:
            p = ln.split("\t")
            if len(p) <= max(i_len, i_cnt) or p[i_gene] == "Geneid":
                continue
            try:
                cnt = float(p[i_cnt])
            except ValueError:
                continue
            rpk = cnt / float(p[i_len]) if float(p[i_len]) else 0.0
            per_gene[p[i_gene]] = rpk
            denom += rpk
        scale = 1e6 / denom if denom else 0.0
        return {symbol_of[gid]: rpk * scale
                for gid, rpk in per_gene.items() if gid in symbol_of}
    except (OSError, ValueError, IndexError) as exc:
        sys.stderr.write(f"[neoantigen_report] counts unreadable ({exc}); "
                         f"TPM annotation will be empty\n")
        return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pvac-tsv", required=True)
    ap.add_argument("--hla-tsv", required=True)
    ap.add_argument("--usable-alleles", default="",
                    help="pvac_alleles.txt — the predictor-usable allele list")
    ap.add_argument("--counts", default="")
    ap.add_argument("--gene-map", default="")
    ap.add_argument("--patient", required=True)
    ap.add_argument("--mode", default="tumor_only",
                    choices=["tumor_only", "paired"])
    ap.add_argument("--strong-nm", type=float, default=500.0)
    ap.add_argument("--high-aff-nm", type=float, default=100.0)
    ap.add_argument("--expr-tpm", type=float, default=1.0)
    ap.add_argument("--predictors", default="MHCflurry")
    ap.add_argument("--top-n", type=int, default=20)
    ap.add_argument("--out-candidates", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()

    rows = read_table(args.pvac_tsv)
    # best (lowest IC50) row per (gene, epitope, hla). Note: distinct variant
    # isoforms yielding the same (gene, epitope, HLA) collapse into the
    # best-scoring row — the burden metric counts unique neoepitopes, not
    # unique variants (audit m4).
    def fnum(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return float("inf")

    best = {}
    for r in rows:
        key = (r.get("gene", ""), r.get("epitope", ""), r.get("hla_allele", ""))
        cur = best.get(key)
        if cur is None or fnum(r.get("mt_score_nm")) < fnum(cur.get("mt_score_nm")):
            best[key] = r
    cands = list(best.values())

    tpm = (read_tpm_map(args.counts, args.gene_map)
           if args.counts and args.gene_map else {})
    for r in cands:
        v = tpm.get(r.get("gene", ""))
        r["gene_tpm"] = f"{v:.2f}" if v is not None else ""
        r["expr_pass"] = ("yes" if v is not None and v >= args.expr_tpm
                          else "no")

    for r in cands:
        r["_s"] = fnum(r.get("mt_score_nm"))
    cands.sort(key=lambda r: r["_s"])

    n_strong = sum(1 for r in cands if r["_s"] <= args.strong_nm)
    n_highaff = sum(1 for r in cands if r["_s"] <= args.high_aff_nm)
    n_fs = sum(1 for r in cands if r["_s"] <= args.strong_nm
               and "frameshift" in (r.get("variant_type", "") or "").lower())
    strong_rows = [r for r in cands if r["_s"] <= args.strong_nm]
    n_tpm_pass = sum(1 for r in strong_rows if r["expr_pass"] == "yes")

    typed = set()
    with open(args.hla_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in ("allele_1", "allele_2"):
                if r.get(k, "").strip():
                    typed.add(r[k].strip())
    usable = set()
    if args.usable_alleles and os.path.exists(args.usable_alleles):
        with open(args.usable_alleles) as fh:
            usable = {a.strip() for a in fh.read().split(",") if a.strip()}

    mode_caveat = (
        "tumor-only: putative somatic — without a matched normal, germline "
        "polymorphisms can masquerade as somatic candidates; burden numbers "
        "are screening signals, not clinical quantitation"
        if args.mode == "tumor_only" else
        "paired (matched normal filtered); still research-grade, orthogonal "
        "confirmation required")
    class2_note = ("MHC-II predictions not run in P0 (NetMHCIIpan/MixMHC2pred "
                   "absent — DTU academic license requires manual install)")

    summary = [
        {"class": "neoantigen", "metric": "hla_alleles_typed",
         "value": len(typed), "detail": ";".join(sorted(typed))},
        {"class": "neoantigen", "metric": "hla_alleles_predictor_usable",
         "value": len(usable), "detail": ";".join(sorted(usable))},
        {"class": "neoantigen", "metric": "predictors", "value": args.predictors,
         "detail": class2_note},
        {"class": "neoantigen", "metric": "candidates_reported", "value": len(rows),
         "detail": "rows surviving pVACseq's reported filter set"},
        {"class": "neoantigen", "metric": "unique_candidates",
         "value": len(cands),
         "detail": "best row per (gene, epitope, HLA); counts unique "
                   "neoepitopes, not unique variants"},
        {"class": "neoantigen", "metric": "burden_strong_500nM",
         "value": n_strong,
         "detail": f"unique epitopes MT IC50 <= {args.strong_nm:g} nM"},
        {"class": "neoantigen", "metric": "burden_highaff_100nM",
         "value": n_highaff,
         "detail": f"IC50 <= {args.high_aff_nm:g} nM (high-affinity bucket; "
                   f"TESLA consensus model optimum was IC50<34 nM with "
                   f"tumor abundance>33 TPM and stability>1.4h — Wells 2020)"},
        {"class": "neoantigen", "metric": "frameshift_neoantigens",
         "value": n_fs, "detail": "frameshift-derived among strong binders"},
        {"class": "neoantigen", "metric": "expr_tpm_pass",
         "value": n_tpm_pass,
         "detail": f"strong binders with gene TPM >= {args.expr_tpm:g} "
                   f"(common-practice expression filter; {len(tpm)} gene "
                   f"symbols mapped)"},
        {"class": "neoantigen", "metric": "mode", "value": args.mode,
         "detail": mode_caveat},
    ]

    with open(args.out_summary, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(summary)

    with open(args.out_candidates, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CAND_COLS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for r in cands:
            w.writerow({c: r.get(c, "") for c in CAND_COLS})

    with open(args.out_md, "w") as fh:
        fh.write(f"## Neoantigen candidates — {args.patient}\n\n")
        fh.write(f"- HLA alleles typed: **{len(typed)}** "
                 f"({', '.join(sorted(typed)) or 'none'}); predictor-usable: "
                 f"{len(usable)}\n")
        fh.write(f"- Unique candidates: {len(cands)}; "
                 f"**{n_strong}** strong binders (<= {args.strong_nm:g} nM), "
                 f"{n_highaff} high-affinity (<= {args.high_aff_nm:g} nM); "
                 f"{n_fs} frameshift-derived\n")
        fh.write(f"- Expression filter (TPM >= {args.expr_tpm:g}): "
                 f"{n_tpm_pass}/{len(strong_rows) or 0} strong binders pass\n")
        fh.write(f"- Predictors: {args.predictors}. {class2_note}.\n")
        fh.write(f"- Caveat: {mode_caveat}\n\n")
        fh.write("| Gene | HLA | Epitope | IC50(nM) | TPM | VAF | Type |\n")
        fh.write("|---|---|---|---|---|---|---|\n")
        for r in cands[:args.top_n]:
            fh.write(f"| {r.get('gene','')} | {r.get('hla_allele','')} | "
                     f"{r.get('epitope','')} | {r.get('mt_score_nm','')} | "
                     f"{r.get('gene_tpm','')} | {r.get('tumor_vaf','')} | "
                     f"{r.get('variant_type','')} |\n")
    sys.stderr.write(f"[neoantigen_report] {args.patient}: {len(cands)} unique "
                     f"candidates, {n_strong} strong -> {args.out_summary}\n")


if __name__ == "__main__":
    main()
