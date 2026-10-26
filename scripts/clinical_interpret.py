#!/usr/bin/env python3
"""Per-patient clinical interpretation: merge evidence streams for venus.

Combines SNV tiering (clinical_tier.py), CNV gene calls (clinical_cnv.py),
fusion candidates (clinical_sv.py), and cohort-level biomarkers (TMB, MSI,
signatures) into:
  {pair_id}.biomarkers.tsv  — actionable biomarker list, AMP tier sorted
  {pair_id}.therapies.tsv   — therapy options per biomarker, with provenance

Therapy matching is BIOMARKER-KEYED, not gene-keyed: the curated NMPA/CSCO
table (resources/clinical/nmpa_oncology_approvals.tsv) carries a
`biomarker_key` (allele/fusion/amp + required confirmation modality) and a
`requires` column that is surfaced on every therapy row. Guard rails:
  - allele-specific keys (KRAS G12C, BRAF V600E, PIK3CA hotspots) match only
    the named protein change, never other variants in the same gene
  - MSI/TMB therapies attach only to the POSITIVE status (never MSS / not-H)
  - fusion therapies require a RESOLVED two-gene call; single-breakpoint
    candidateSV hits are "screening signals" (un-tiered, no therapy mapping)
  - ERBB2 amplification carries an IHC/FISH-concordance requirement
  - context_only rows (trial-only / historical-threshold) are emitted with
    context_only=true and never as Tier I/II evidence

All tumor-only outputs carry the putative-somatic caveat. Missing inputs
(e.g. MSI in tumor-only cohorts) degrade to honest not_run rows.
"""
import argparse
import os
import sys

import pandas as pd

TIER_ORDER = {"I": 0, "II": 1, "III": 2, "IV": 3, "": 4}
TMB_HIGH_FDA = 13.0        # FDA label since 2026-10-16 traditional approval
TMB_HIGH_LEGACY = 10.0     # historical AMP/ASCO/CAP context (no longer label)
VAF_REPORT_FLOOR = 0.05
DP_REPORT_FLOOR = 20
ALLELE_NMPA_KEYS = {
    ("KRAS", "G12C"): "KRAS G12C snv",
    ("BRAF", "V600E"): "BRAF V600E snv",
    ("PIK3CA", "E545K"): "PIK3CA hotspot snv",
    ("PIK3CA", "H1047R"): "PIK3CA hotspot snv",
}


def read_tsv(path):
    if path and os.path.exists(path):
        try:
            return pd.read_csv(path, sep="\t", dtype=str).fillna("")
        except Exception as exc:
            sys.stderr.write(f"[clinical_interpret] WARN unreadable {path}: {exc}\n")
    return pd.DataFrame()


def _num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def nmpa_match(key, nmpa):
    if not key:
        return pd.DataFrame()
    return nmpa[nmpa["biomarker_key"] == key]


def fusion_nmpa_key(gene):
    """Resolved-partner gene -> NMPA therapy key (audit M3): the curated
    table folds NTRK1/2/3 under one tissue-agnostic row."""
    return ("NTRK fusion resolved" if str(gene).upper().startswith("NTRK")
            else f"{gene} fusion resolved")


def emit_therapies(rows_out, biomarker, civic_therapies, nmpa, nmpa_key,
                   kb_tag, klass):
    for t in str(civic_therapies or "").split(";"):
        t = t.strip()
        if t:
            rows_out.append({
                "biomarker": biomarker, "therapy": t, "target_or_context": klass,
                "region_or_approval": "CIViC evidence", "line_or_level": "",
                "requires": "", "source": "CIViC", "kb_tag": kb_tag})
    if nmpa_key:
        for _, n in nmpa_match(nmpa_key, nmpa).iterrows():
            if str(n.get("context_only", "")).lower() == "true":
                region = f"{n['approval_region']} [context only]"
            else:
                region = n["approval_region"]
            rows_out.append({
                "biomarker": biomarker, "therapy": n["therapy"],
                "target_or_context": n["biomarker_key"],
                "region_or_approval": region,
                "line_or_level": f"{n['line']} | {n['guideline_context']}",
                "requires": n.get("requires", ""),
                "source": f"NMPA/CSCO curated ({n.get('as_of', '')})",
                "kb_tag": kb_tag})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", required=True)
    ap.add_argument("--snv-evidence", required=True)
    ap.add_argument("--cnv-evidence", required=True)
    ap.add_argument("--sv-evidence", required=True)
    ap.add_argument("--tmb", default="")
    ap.add_argument("--msi", default="")
    ap.add_argument("--signatures", default="")
    ap.add_argument("--nmpa", required=True)
    ap.add_argument("--cohort-mode", default="tumor_only", choices=["tumor_only", "paired"],
                    help="affects somatic-quality caveats (audit m1)")
    ap.add_argument("--kb-tag", default="")
    ap.add_argument("--out-biomarkers", required=True)
    ap.add_argument("--out-therapies", required=True)
    args = ap.parse_args()

    snv = read_tsv(args.snv_evidence)
    cnv = read_tsv(args.cnv_evidence)
    sv = read_tsv(args.sv_evidence)
    nmpa = read_tsv(args.nmpa)

    biomarkers = []
    therapies = []

    # --- SNV/indel evidence ---
    n_snv = n_matched = 0
    if snv.empty:
        biomarkers.append({"class": "SNV/indel", "biomarker": "not_run",
                           "amp_tier": "", "detail": "no SNV evidence table",
                           "evidence": "", "therapies_civic": "",
                           "clinvar_sig": "", "caveat": "",
                           "civic_diseases": "", "evidence_types": ""})
    else:
        for _, r in snv.iterrows():
            n_snv += 1
            if r.get("civic_variant_id"):
                n_matched += 1
        keep = snv[snv["amp_tier"].isin(["I", "II", "IV"]) | snv["civic_variant_id"].astype(bool)]
        for _, r in keep.iterrows():
            detail = (f"VAF={r['tumor_vaf']} DP={r['tumor_dp']} "
                      f"gnomAD_AF={r['gnomad_af']} cons={r['consequence']}")
            if r["clinvar_sig"]:
                caveat = f"germline context: {r['clinvar_sig']}"
            else:
                caveat = ("putative somatic (tumor-only)" if args.cohort_mode == "tumor_only"
                          else "somatic (matched-normal filtered)")
            vaf, dp = _num(r["tumor_vaf"]), _num(r["tumor_dp"])
            if (vaf == vaf and vaf < VAF_REPORT_FLOOR) or (dp == dp and dp < DP_REPORT_FLOOR):
                caveat += (f" [below common reporting floor VAF>={VAF_REPORT_FLOOR}, "
                           f"DP>={DP_REPORT_FLOOR}]")
            allele_key = ALLELE_NMPA_KEYS.get((r["gene"], r["hgvsp"]), "")
            biomarkers.append({
                "class": "SNV/indel",
                "biomarker": f"{r['gene']} {r['hgvsp'] or r['alt']} "
                             f"({r['chrom']}:{r['pos']}{r['ref']}>{r['alt']})",
                "amp_tier": r["amp_tier"], "detail": detail,
                "evidence": r["evidence_source"],
                "therapies_civic": r["civic_therapies"],
                "clinvar_sig": r["clinvar_sig"],
                "caveat": caveat,
                "civic_diseases": r.get("civic_diseases", ""),
                "evidence_types": r.get("civic_evidence_types", ""),
                "nmpa_key": allele_key,
            })
            emit_therapies(therapies, biomarkers[-1]["biomarker"],
                           r["civic_therapies"], nmpa, allele_key,
                           args.kb_tag, "SNV/indel")
        if keep.empty:
            biomarkers.append({
                "class": "SNV/indel",
                "biomarker": f"no retained calls ({n_snv} tiered, {n_matched} CIViC-matched)",
                "amp_tier": "", "detail": "no Tier I/II and no CIViC-matched variant",
                "evidence": "", "therapies_civic": "", "clinvar_sig": "",
                "caveat": "", "civic_diseases": "", "evidence_types": ""})

    # --- CNV gene calls (amplification-type rows carry "<gene> amp" keys) ---
    if cnv.empty:
        biomarkers.append({"class": "CNV", "biomarker": "no_panel_calls",
                           "amp_tier": "", "detail": "no actionable CNV gene call",
                           "evidence": "", "therapies_civic": "",
                           "clinvar_sig": "", "caveat": "",
                           "civic_diseases": "", "evidence_types": ""})
    else:
        for _, r in cnv.iterrows():
            ec = f" + ecDNA evidence" if r.get("ecDNA_evidence") == "yes" else ""
            is_amp = "amplification" in r["call"]
            is_del = "deletion" in r["call"]
            nmpa_key = f"{r['gene']} amp" if is_amp else ""
            # audit M4: BRCA1/2 biallelic loss surfaces the PARP trial-context
            # row as context only — never as germline-grade pathogenicity.
            if is_del and r["gene"] in ("BRCA1", "BRCA2"):
                nmpa_key = "BRCA1/2 loss (somatic, context)"
            caveat = ("putative somatic (tumor-only)" if args.cohort_mode == "tumor_only"
                      else "somatic (matched-normal filtered)")
            if r["gene"] == "ERBB2" and is_amp:
                caveat += "; HER2 positivity requires IHC3+/ISH+ concordance (DNA amp alone insufficient)"
            # audit M5: tier is per-gene evidence from the panel table, not a
            # blanket Tier II.
            tier = str(r.get("amp_tier", "") or "III")
            biomarkers.append({
                "class": "CNV",
                "biomarker": f"{r['gene']} {r['call']}",
                "amp_tier": tier,
                "detail": f"{r['segment']} log2={r['segment_log2']} cn={r['segment_cn']}{ec}",
                "evidence": "CNVkit+panel",
                "therapies_civic": "",
                "clinvar_sig": "", "caveat": caveat,
                "civic_diseases": "", "evidence_types": "",
                "nmpa_key": nmpa_key,
            })
            emit_therapies(therapies, biomarkers[-1]["biomarker"], "", nmpa,
                           nmpa_key, args.kb_tag, "CNV")

    # --- Fusion candidates: resolved = Tier II + therapy mapping; ---
    # unresolved single-breakpoint candidateSV = un-tiered screening signal.
    if sv.empty:
        biomarkers.append({"class": "SV/fusion", "biomarker": "no_panel_candidates",
                           "amp_tier": "", "detail": "no actionable fusion candidate",
                           "evidence": "", "therapies_civic": "",
                           "clinvar_sig": "", "caveat": "",
                           "civic_diseases": "", "evidence_types": ""})
    else:
        for _, r in sv.iterrows():
            resolved = bool(r.get("gene_5_call")) and bool(r.get("gene_3_call"))
            if resolved:
                # audit M3/M5: fold NTRK1/2/3 onto the curated NTRK row; tier
                # II only when an approved (non-context) therapy row matched.
                genes = (r["gene_5_call"], r["gene_3_call"])
                keys = sorted({fusion_nmpa_key(g) for g in genes})
                has_approved = any(
                    str(n.get("context_only", "")).lower() != "true"
                    for k in keys for _, n in nmpa_match(k, nmpa).iterrows())
                fmode = ("tumor-only; " if args.cohort_mode == "tumor_only" else "")
                biomarkers.append({
                    "class": "SV/fusion",
                    "biomarker": f"{r['fusion_candidate']} fusion",
                    "amp_tier": "II" if has_approved else "",
                    "detail": f"{r['breakpoint1']} / {r['breakpoint2']} [{r['scoring']}]",
                    "evidence": "curated fusion panel",
                    "therapies_civic": "", "clinvar_sig": "",
                    "caveat": f"{fmode}confirm by RNA/FISH",
                    "civic_diseases": "", "evidence_types": "",
                    "nmpa_key": ";".join(keys),
                })
                for k in keys:
                    emit_therapies(therapies, biomarkers[-1]["biomarker"], "",
                                   nmpa, k, args.kb_tag, "SV/fusion")
            else:
                ghit = r.get("gene_5_call") or r.get("gene_3_call")
                biomarkers.append({
                    "class": "SV/fusion screening",
                    "biomarker": f"{r['fusion_candidate']} (screening signal)",
                    "amp_tier": "",
                    "detail": (f"{r['breakpoint1']} / {r['breakpoint2']} "
                               f"[{r['scoring']}] — un-resolved single breakpoint "
                               f"in {ghit}; NOT evidence of an actionable fusion"),
                    "evidence": "candidateSV coordinate screen",
                    "therapies_civic": "", "clinvar_sig": "",
                    "caveat": "no therapy mapping; confirm partner by RNA/long-read",
                    "civic_diseases": "", "evidence_types": "",
                })

    # --- TMB: dual-threshold (FDA >=13 label since 2026-10-16; >=10 legacy),
    # with germline-contamination quantification from the SNV evidence table.
    germline_flagged = germline_af = ""
    if not snv.empty:
        flag = snv["clinvar_sig"].str.contains(
            "Pathogenic", case=False, na=False)
        germline_flagged = int(flag.sum())
        af = pd.to_numeric(snv["gnomad_af"], errors="coerce")
        germline_af = int((af > 0.001).sum())
    tmb = read_tsv(args.tmb)
    if tmb.empty:
        biomarkers.append({"class": "TMB", "biomarker": "not_run", "amp_tier": "",
                           "detail": "TMB table missing", "evidence": "",
                           "therapies_civic": "", "clinvar_sig": "", "caveat": "",
                           "civic_diseases": "", "evidence_types": ""})
    else:
        row = tmb.iloc[0]
        val = _num(row.get("tmb_mut_per_mb"))
        if val != val:
            biomarkers.append({"class": "TMB", "biomarker": "unreadable TMB value",
                               "amp_tier": "",
                               "detail": f"raw row: {row.to_dict()}",
                               "evidence": "", "therapies_civic": "",
                               "clinvar_sig": "", "caveat": "",
                               "civic_diseases": "", "evidence_types": ""})
        else:
            germline_note = ""
            if germline_flagged != "":
                gmode = ("tumor-only, no PoN" if args.cohort_mode == "tumor_only"
                         else "matched-normal cohort")
                germline_note = (f"; germline screen: {germline_flagged} ClinVar-"
                                 f"flagged + {germline_af} gnomAD_AF>0.001 of "
                                 f"{n_snv} coding rows ({gmode})")
            if val >= TMB_HIGH_FDA:
                tier, bio = "I", f"TMB-H (FDA >=13: {val} mut/Mb coding)"
                nmpa_key = "TMB-H (FDA >=13 mut/Mb)"
            elif val >= TMB_HIGH_LEGACY:
                tier, bio = "", (f"TMB borderline (10-<13 legacy: {val} mut/Mb)")
                nmpa_key = "TMB borderline (10-<13 mut/Mb)"
            else:
                tier, bio, nmpa_key = "", f"not TMB-H ({val} mut/Mb coding)", ""
            biomarkers.append({
                "class": "TMB", "biomarker": bio, "amp_tier": tier,
                "detail": (f"coding {row.get('coding_mb', '')}Mb denominator; "
                           f"all-somatic {row.get('all_somatic_mut_per_mb', '')} "
                           f"mut/Mb (research){germline_note}"),
                "evidence": f"FDA label >=13 (2026-10-16); legacy >=10 (AMP/ASCO/CAP)",
                "therapies_civic": "immune checkpoint inhibitor context" if tier else "",
                "clinvar_sig": "",
                "caveat": ("tumor-only TMB overestimates vs paired"
                           if args.cohort_mode == "tumor_only"
                           else "matched-normal cohort TMB"),
                "civic_diseases": "", "evidence_types": "",
            })
            emit_therapies(therapies, bio,
                           "immune checkpoint inhibitor" if tier else "",
                           nmpa, nmpa_key, args.kb_tag, "TMB")

    # --- MSI ---
    msi = read_tsv(args.msi)
    if msi.empty:
        biomarkers.append({"class": "MSI", "biomarker": "not_run", "amp_tier": "",
                           "detail": "MSI requires paired normal (tumor-only cohort)",
                           "evidence": "", "therapies_civic": "", "clinvar_sig": "",
                           "caveat": "", "civic_diseases": "", "evidence_types": ""})
    else:
        cols = {c.lower(): c for c in msi.columns}
        status = (msi.iloc[0].get(cols.get("msi_status", ""), "")
                  if "msi_status" in cols else str(msi.iloc[0].to_dict()))
        is_msih = "msi-h" in str(status).lower()
        biomarkers.append({
            "class": "MSI", "biomarker": f"MSI status: {status}",
            "amp_tier": "I" if is_msih else "",
            "detail": "MSIsensor-pro", "evidence": "",
            "therapies_civic": "ICI (pembro/nivo±ipi/envafolimab)" if is_msih else "",
            "clinvar_sig": "", "caveat": "",
            "civic_diseases": "", "evidence_types": "",
        })
        if is_msih:
            emit_therapies(therapies, f"MSI status: {status}", "ICI", nmpa,
                           "MSI-H positive", args.kb_tag, "MSI")

    # --- SBS signature context ---
    sig = read_tsv(args.signatures)
    if not sig.empty:
        if "Signature" in sig.columns:
            frac = pd.to_numeric(sig.get("Exposure_Fraction"), errors="coerce")
            top2 = sig.assign(_f=frac.fillna(0)).sort_values("_f", ascending=False).head(2)
            detail = ", ".join(f"{r['Signature']}={r['_f']:.2f}"
                               for _, r in top2.iterrows())
        else:
            top = sig.iloc[0].to_dict()
            pairs_ = sorted(((k, v) for k, v in top.items() if v),
                            key=lambda kv: -float(kv[1] or 0))[:2]
            detail = ", ".join(f"{k}={float(v):.2f}" for k, v in pairs_)
        biomarkers.append({
            "class": "SBS signature", "biomarker": "exposure context",
            "amp_tier": "", "detail": detail, "evidence": "COSMIC v3.1 NNLS",
            "therapies_civic": "", "clinvar_sig": "",
            "caveat": "context only (e.g. SBS17 = oxidative/gastric)" if "SBS17" in detail else "",
            "civic_diseases": "", "evidence_types": "",
        })

    bm = pd.DataFrame(biomarkers).fillna("")
    if not bm.empty:
        bm["_o"] = bm["amp_tier"].map(TIER_ORDER)
        bm = bm.sort_values(["_o", "class", "biomarker"]).drop(columns="_o")
    bm.to_csv(args.out_biomarkers, sep="\t", index=False)

    # dedupe identical therapy rows (same biomarker+therapy+region+requires)
    th = pd.DataFrame(therapies, columns=[
        "biomarker", "therapy", "target_or_context", "region_or_approval",
        "line_or_level", "requires", "source", "kb_tag"])
    if not th.empty:
        th = th.drop_duplicates()
    th.to_csv(args.out_therapies, sep="\t", index=False)
    sys.stderr.write(f"[clinical_interpret] {len(bm)} biomarkers, "
                     f"{len(th)} therapy rows -> {args.out_biomarkers}\n")


if __name__ == "__main__":
    main()
