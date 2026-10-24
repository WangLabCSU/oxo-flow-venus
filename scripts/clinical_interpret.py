#!/usr/bin/env python3
"""Per-patient clinical interpretation: merge evidence streams for venus.

Combines SNV tiering (clinical_tier.py), CNV gene calls (clinical_cnv.py),
fusion candidates (clinical_sv.py), and cohort-level biomarkers (TMB, MSI,
signatures) into:
  {pair_id}.biomarkers.tsv  — actionable biomarker list, AMP tier sorted
  {pair_id}.therapies.tsv   — therapy options per biomarker, with provenance
                               (CIViC assertion/evidence + NMPA/CSCO table)

All tumor-only outputs carry the putative-somatic caveat. Missing inputs
(e.g. MSI in tumor-only cohorts) degrade to honest not_run rows.
"""
import argparse
import os
import sys

import pandas as pd

TIER_ORDER = {"I": 0, "II": 1, "III": 2, "IV": 3, "": 4}
TMB_HIGH_MUT_MB = 10.0  # AMP/ASCO/CAP & FDA tissue-TMB threshold (mut/Mb, coding)


def read_tsv(path):
    if path and os.path.exists(path):
        try:
            return pd.read_csv(path, sep="\t", dtype=str).fillna("")
        except Exception as exc:
            sys.stderr.write(f"[clinical_interpret] WARN unreadable {path}: {exc}\n")
    return pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", required=True)
    ap.add_argument("--snv-evidence", required=True)
    ap.add_argument("--cnv-evidence", required=True)
    ap.add_argument("--sv-evidence", required=True)
    ap.add_argument("--tmb", default="", help="report/tmb/{pair_id}.tmb.tsv")
    ap.add_argument("--msi", default="", help="deliver/{pair_id}.msi.tsv (paired only)")
    ap.add_argument("--signatures", default="", help="SBS exposures TSV")
    ap.add_argument("--nmpa", required=True, help="nmpa_oncology_approvals.tsv")
    ap.add_argument("--kb-tag", default="")
    ap.add_argument("--out-biomarkers", required=True)
    ap.add_argument("--out-therapies", required=True)
    args = ap.parse_args()

    snv = read_tsv(args.snv_evidence)
    cnv = read_tsv(args.cnv_evidence)
    sv = read_tsv(args.sv_evidence)
    nmpa = read_tsv(args.nmpa)

    biomarkers = []

    # --- SNV/indel evidence ---
    if snv.empty:
        biomarkers.append({"class": "SNV/indel", "biomarker": "not_run",
                           "amp_tier": "", "detail": "no SNV evidence table"})
    else:
        keep = snv[snv["amp_tier"].isin(["I", "II", "IV"]) | snv["civic_variant_id"].astype(bool)]
        for _, r in keep.iterrows():
            detail = (f"VAF={r['tumor_vaf']} DP={r['tumor_dp']} "
                      f"gnomAD_AF={r['gnomad_af']} cons={r['consequence']}")
            biomarkers.append({
                "class": "SNV/indel",
                "biomarker": f"{r['gene']} {r['hgvsp'] or r['alt']} "
                             f"({r['chrom']}:{r['pos']}{r['ref']}>{r['alt']})",
                "amp_tier": r["amp_tier"], "detail": detail,
                "evidence": r["evidence_source"],
                "therapies_civic": r["civic_therapies"],
                "clinvar_sig": r["clinvar_sig"],
                "caveat": "putative somatic (tumor-only)" if not r["clinvar_sig"] else
                          f"germline context: {r['clinvar_sig']}",
            })

    # --- CNV gene calls ---
    if cnv.empty:
        biomarkers.append({"class": "CNV", "biomarker": "no_panel_calls",
                           "amp_tier": "", "detail": "no actionable CNV gene call"})
    else:
        for _, r in cnv.iterrows():
            ec = f" + ecDNA evidence" if r.get("ecDNA_evidence") == "yes" else ""
            biomarkers.append({
                "class": "CNV",
                "biomarker": f"{r['gene']} {r['call']}",
                "amp_tier": "II" if "amplification" in r["call"] else "II",
                "detail": f"{r['segment']} log2={r['segment_log2']} cn={r['segment_cn']}{ec}",
                "evidence": "CNVkit+panel",
                "therapies_civic": "",
                "clinvar_sig": "", "caveat": "putative somatic (tumor-only)",
            })

    # --- Fusion candidates ---
    if sv.empty:
        biomarkers.append({"class": "SV/fusion", "biomarker": "no_panel_candidates",
                           "amp_tier": "", "detail": "no actionable fusion candidate"})
    else:
        for _, r in sv.iterrows():
            biomarkers.append({
                "class": "SV/fusion",
                "biomarker": f"{r['fusion_candidate']} fusion candidate",
                "amp_tier": "II",
                "detail": f"{r['breakpoint1']} / {r['breakpoint2']} [{r['scoring']}]",
                "evidence": r.get("evidence_context", ""),
                "therapies_civic": "",
                "clinvar_sig": "", "caveat": "putative somatic (tumor-only); unscored SV",
            })

    # --- Cohort-level biomarkers ---
    tmb = read_tsv(args.tmb)
    if tmb.empty:
        biomarkers.append({"class": "TMB", "biomarker": "not_run", "amp_tier": "",
                           "detail": "TMB table missing", "evidence": "",
                           "therapies_civic": "", "clinvar_sig": "", "caveat": ""})
    else:
        row = tmb.iloc[0]
        val = float(row["tmb_mut_per_mb"])
        flag = "TMB-H" if val >= TMB_HIGH_MUT_MB else "not TMB-H"
        biomarkers.append({
            "class": "TMB", "biomarker": f"{flag} ({val} mut/Mb coding)",
            "amp_tier": "I" if val >= TMB_HIGH_MUT_MB else "",
            "detail": f"coding {row['coding_mb']}Mb denominator; "
                      f"all-somatic {row['all_somatic_mut_per_mb']} mut/Mb (research)",
            "evidence": "AMP/ASCO/CAP >=10 mut/Mb",
            "therapies_civic": "immune checkpoint inhibitor context" if val >= TMB_HIGH_MUT_MB else "",
            "clinvar_sig": "", "caveat": "tumor-only TMB overestimates vs paired",
        })

    msi = read_tsv(args.msi)
    if msi.empty:
        biomarkers.append({"class": "MSI", "biomarker": "not_run", "amp_tier": "",
                           "detail": "MSI requires paired normal (tumor-only cohort)",
                           "evidence": "", "therapies_civic": "", "clinvar_sig": "",
                           "caveat": ""})
    else:
        cols = {c.lower(): c for c in msi.columns}
        status = msi.iloc[0].get(cols.get("msi_status", "msi_status"), "") \
            if "msi_status" in cols else str(msi.iloc[0].to_dict())
        biomarkers.append({
            "class": "MSI", "biomarker": f"MSI status: {status}",
            "amp_tier": "I" if "msi-h" in str(status).lower() else "",
            "detail": "MSIsensor-pro", "evidence": "",
            "therapies_civic": "ICI (pembro/nivo±ipi/envolimab)" if "msi-h" in str(status).lower() else "",
            "clinvar_sig": "", "caveat": "",
        })

    sig = read_tsv(args.signatures)
    if not sig.empty:
        if "Signature" in sig.columns:  # long format: Signature|Attributed_Mutations|Exposure_Fraction
            frac = pd.to_numeric(sig.get("Exposure_Fraction"), errors="coerce")
            top2 = sig.assign(_f=frac.fillna(0)).sort_values("_f", ascending=False).head(2)
            detail = ", ".join(f"{r['Signature']}={r['_f']:.2f}"
                               for _, r in top2.iterrows())
        else:  # wide fallback: one row, signature-named columns
            top = sig.iloc[0].to_dict()
            pairs_ = sorted(((k, v) for k, v in top.items() if v),
                            key=lambda kv: -float(kv[1] or 0))[:2]
            detail = ", ".join(f"{k}={float(v):.2f}" for k, v in pairs_)
        biomarkers.append({
            "class": "SBS signature", "biomarker": "exposure context",
            "amp_tier": "", "detail": detail, "evidence": "COSMIC v3.1 NNLS",
            "therapies_civic": "", "clinvar_sig": "",
            "caveat": "context only (e.g. SBS17 = oxidative/gastric)" if "SBS17" in detail else "",
        })

    bm = pd.DataFrame(biomarkers)
    if not bm.empty:
        bm["_o"] = bm["amp_tier"].map(TIER_ORDER)
        bm = bm.sort_values(["_o", "class", "biomarker"]).drop(columns="_o")
    bm.to_csv(args.out_biomarkers, sep="\t", index=False)

    # --- Therapy options (CIViC therapies + NMPA/CSCO curated context) ---
    therapies = []
    for _, r in bm.iterrows():
        if r.get("therapies_civic"):
            for t in str(r["therapies_civic"]).split(";"):
                if t.strip():
                    therapies.append({
                        "biomarker": r["biomarker"], "therapy": t.strip(),
                        "target_or_context": r["class"], "region_or_approval": "CIViC evidence",
                        "line_or_level": "", "source": "CIViC", "kb_tag": args.kb_tag,
                    })
        gene_key = str(r["biomarker"]).split()[0]
        hits = nmpa[(nmpa["gene"] == gene_key) |
                    (nmpa["biomarker"].str.startswith(gene_key + " "))]
        for _, n in hits.iterrows():
            therapies.append({
                "biomarker": r["biomarker"], "therapy": n["therapy"],
                "target_or_context": n["biomarker"], "region_or_approval": n["approval_region"],
                "line_or_level": f"{n['line']} | {n['guideline_context']}",
                "source": "NMPA/CSCO curated (" + n["as_of"] + ")",
                "kb_tag": args.kb_tag,
            })

    th = pd.DataFrame(therapies)
    th.to_csv(args.out_therapies, sep="\t", index=False)
    sys.stderr.write(f"[clinical_interpret] {len(bm)} biomarkers, "
                     f"{len(th)} therapy rows -> {args.out_biomarkers}\n")


if __name__ == "__main__":
    main()
