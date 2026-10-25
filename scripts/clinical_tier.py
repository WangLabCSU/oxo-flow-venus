#!/usr/bin/env python3
"""SNV/indel clinical evidence extraction + AMP/ASCO/CAP tiering for venus.

Input: the rule's bcftools-annotated PASS+VEP VCF (ClinVar ID/CLNSIG joined by
coordinate, gnomAD AF joined as INFO/AF) plus CIViC monthly dump TSVs.

CIViC coordinates are GRCh37 while the pipeline is GRCh38, so matching is
build-independent (same strategy as the OncoKB annotator):
  1. gene + normalized protein change (HGVSp core, 1-letter aa)  [primary]
  2. ClinVar VariationID bridge (VCF ID joined from clinvar.vcf.gz,
     compared against CIViC clinvar_ids)                          [secondary]
  3. GRCh38 coordinate equality                                   [tertiary]

AMP tiers follow the VICC meta-knowledgebase mapping (Wagner et al., Nat Genet
2020, PMC7127986): CIViC assertion amp_category passthrough when present,
otherwise evidence level A->Tier I-A, B->I-B, C->II-C, D->II-D, E->III;
benign-only assertion -> Tier IV; unmatched protein-altering -> Tier III.

Output: long TSV, one row per protein-altering or KB-matched PASS variant.
Matching is substitution-centric: missense/nonsense protein changes get exact
keys; in-frame indels/delins (e.g. EGFR E746_A750del, MET exon-14 skips) have
no portable protein key yet and currently match only via ClinVar VariationID —
a known recall gap, disclosed here and in docs/refdata.md §10.
Stdlib + pandas + pysam; CSQ subfield order read from the VCF header
(same contract as scripts/vcf2maf_light.py).
"""
import argparse
import re
import sys

import pandas as pd
import pysam

AA3TO1 = {
    "ala": "A", "arg": "R", "asn": "N", "asp": "D", "cys": "C", "gln": "Q",
    "glu": "E", "gly": "G", "his": "H", "ile": "I", "leu": "L", "lys": "K",
    "met": "M", "phe": "F", "pro": "P", "ser": "S", "thr": "T", "trp": "W",
    "tyr": "Y", "val": "V", "ter": "*",
}
CODE_KEEP = {"missense_variant", "synonymous_variant", "stop_gained", "stop_lost",
             "start_lost", "frameshift_variant", "inframe_insertion",
             "inframe_deletion", "protein_altering_variant",
             "splice_acceptor_variant", "splice_donor_variant",
             "transcript_ablation"}
SUB_RE = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")
P_RE = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")
CSQ_FORMAT_RE = re.compile(r"Format[:=]\s*([^\">]+)")


def _aa1(three):
    return AA3TO1.get(three.lower(), three)


def norm_hgvsp_core(hgvsp):
    """p.Thr315Ile -> T315I ; p.Glu384Ter -> E384* ; '' if unusual form."""
    if not hgvsp:
        return ""
    s = hgvsp.split(":")[-1].strip()
    if s.startswith("p."):
        s = s[2:]
    s = s.split("(")[0].strip()
    m = P_RE.match("p." + s)
    if m:
        return f"{_aa1(m.group(1))}{m.group(2)}{_aa1(m.group(3))}"
    return ""


def civic_protein_keys(name, hgvs_field):
    keys = set()
    m = SUB_RE.match(name or "")
    if m:
        keys.add(f"{m.group(1)}{m.group(2)}{m.group(3)}")
    for h in (hgvs_field or "").split(","):
        m = P_RE.search(h.strip())
        if m:
            keys.add(f"{_aa1(m.group(1))}{m.group(2)}{_aa1(m.group(3))}")
    return keys


def parse_amp_category(raw):
    """'Tier I - Level A' / 'TIER_II_LEVEL_C' / 'Tier IV' -> (tier, level)."""
    if not raw:
        return "", ""
    s = str(raw).upper().replace("_", " ")
    mt = re.search(r"TIER\s*-?\s*(III|IV|II|I)\b", s)
    ml = re.search(r"LEVEL\s*-?\s*([ABCD])\b", s)
    return (mt.group(1) if mt else "", ml.group(1) if ml else "")


def tier_from_evidence_level(lvl):
    return {"A": ("I", "A"), "B": ("I", "B"), "C": ("II", "C"),
            "D": ("II", "D"), "E": ("III", "")}.get(lvl, ("", ""))


def load_civic(variant_tsv, evidence_tsv, assertion_tsv):
    tv = pd.read_csv(variant_tsv, sep="\t", dtype=str).fillna("")
    ev = pd.read_csv(evidence_tsv, sep="\t", dtype=str).fillna("")
    asex = pd.read_csv(assertion_tsv, sep="\t", dtype=str).fillna("")

    ev_by_mp = {str(k): v for k, v in ev.groupby("molecular_profile_id")}
    as_by_mp = {str(k): v for k, v in asex.groupby("molecular_profile_id")}

    variants, by_protein, by_clinvar, by_coord = {}, {}, {}, {}
    for _, r in tv.iterrows():
        vid = r["variant_id"]
        mp = str(r.get("single_variant_molecular_profile_id", ""))
        ev_levels, ev_types, therapies, diseases = set(), set(), set(), set()
        amp_tier, amp_level, amp_raw, reg_approval, nccn, benign = "", "", "", "", "", False
        for sub in (ev_by_mp.get(mp), as_by_mp.get(mp)):
            if sub is None:
                continue
            therapies |= set(t for t in sub["therapies"] if t)
            diseases |= set(d for d in sub["disease"] if d)
            if "evidence_level" in sub.columns:
                ev_levels |= set(sub["evidence_level"])
                ev_types |= set(sub["evidence_type"])
            if "amp_category" in sub.columns:
                for _, a in sub.iterrows():
                    t, l = parse_amp_category(a.get("amp_category", ""))
                    if t and (not amp_tier or (t, l) < (amp_tier, amp_level)):
                        amp_tier, amp_level, amp_raw = t, l, a.get("amp_category", "")
                        reg_approval = a.get("regulatory_approval", "")
                        nccn = a.get("nccn_guideline", "")
                    if "benign" in str(a.get("significance", "")).lower():
                        benign = True
        if not amp_tier and ev_levels:
            amp_tier, amp_level = tier_from_evidence_level(sorted(ev_levels)[0])
        clinv = {x.strip() for x in re.split(r"[,;]", r.get("clinvar_ids", "")) if x.strip()}
        build = r.get("reference_build", "")
        coord = (r.get("chromosome", ""), r.get("start", ""), r.get("stop", ""),
                 r.get("reference_bases", ""), r.get("variant_bases", ""))
        variants[vid] = {
            "gene": r["feature_name"], "civic_variant": r["variant"],
            "protein_keys": civic_protein_keys(r["variant"], r.get("hgvs_descriptions", "")),
            "clinvar_ids": clinv,
            "evidence_levels": ";".join(sorted(ev_levels)),
            "evidence_types": ";".join(sorted(ev_types)),
            "therapies": ";".join(sorted(therapies)),
            "diseases": ";".join(sorted(diseases)),
            "amp_tier": amp_tier, "amp_level": amp_level, "amp_raw": amp_raw,
            "regulatory_approval": reg_approval, "nccn_guideline": nccn,
            "benign": benign,
        }
        for k in variants[vid]["protein_keys"]:
            by_protein.setdefault((r["feature_name"].upper(), k), vid)
        for cid in clinv:
            by_clinvar.setdefault(cid, vid)
        if build == "GRCh38" and coord[0] and coord[1] and coord[3] and coord[4]:
            # normalize away chr prefixes on both sides (CIViC rows use "1")
            by_coord[(coord[0].removeprefix("chr"), coord[1], coord[3], coord[4])] = vid
    return variants, by_protein, by_clinvar, by_coord


def csq_field_names(vcf):
    """Ordered CSQ subfield names from the header (vcf2maf_light contract)."""
    for rec in vcf.header.records:
        if rec.key == "INFO" and rec.get("ID") == "CSQ":
            m = CSQ_FORMAT_RE.search(rec.get("Description", ""))
            if m:
                return [f.strip() for f in m.group(1).split("|")]
    return []


def csq_pick(csq_raw, fields):
    """First CSQ entry (--pick) as a dict keyed by subfield name."""
    if csq_raw is None:
        return None
    first = csq_raw[0] if isinstance(csq_raw, (list, tuple)) else csq_raw
    parts = str(first).split("|")
    return dict(zip(fields, parts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--civic-variant-tsv", required=True)
    ap.add_argument("--civic-evidence-tsv", required=True)
    ap.add_argument("--civic-assertion-tsv", required=True)
    ap.add_argument("--kb-tag", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    variants, by_protein, by_clinvar, by_coord = load_civic(
        args.civic_variant_tsv, args.civic_evidence_tsv, args.civic_assertion_tsv)
    sys.stderr.write(f"[clinical_tier] CIViC variants indexed: {len(variants)}\n")

    vcf = pysam.VariantFile(args.vcf)
    fields = csq_field_names(vcf)
    if not fields:
        sys.exit("[clinical_tier] CSQ Format header not found — not a vep-rs VCF?")
    sample = list(vcf.header.samples)[0] if vcf.header.samples else ""

    rows = []
    for rec in vcf.fetch():
        c = csq_pick(rec.info.get("CSQ"), fields)
        if not c:
            continue
        cons = c.get("Consequence", "")
        if not cons or not CODE_KEEP.intersection(cons.split("&")):
            continue
        gene = c.get("SYMBOL", "")
        hgvsp = norm_hgvsp_core(c.get("HGVSp", ""))
        vaf, dp = "", ""
        smp = rec.samples.get(sample, {}) if sample else {}
        af = smp.get("AF")
        if af not in (None, ()):
            vaf = round(float(af[0]), 4)
        else:
            ad = smp.get("AD")
            if ad and len(ad) == 2 and sum(ad) > 0:
                vaf = round(ad[1] / sum(ad), 4)
        if smp.get("DP"):
            dp = int(smp["DP"])
        gaf = rec.info.get("AF")
        gnomad_af = round(float(gaf[0]), 6) if gaf not in (None, ()) else ""
        clinvar_id = str(rec.id) if rec.id else ""
        clinvar_sig = ";".join(rec.info.get("CLNSIG", []) or [])

        vid = basis = None
        hit = by_protein.get((gene.upper(), hgvsp)) if hgvsp else None
        if hit:
            vid, basis = hit, f"protein:{gene}:{hgvsp}"
        elif clinvar_id and clinvar_id in by_clinvar:
            vid, basis = by_clinvar[clinvar_id], f"clinvar:{clinvar_id}"
        elif (rec.chrom.removeprefix("chr"), str(rec.pos), rec.ref, rec.alts[0]) in by_coord:
            vid = by_coord[(rec.chrom.removeprefix("chr"), str(rec.pos), rec.ref, rec.alts[0])]
            basis = "grch38_coord"
        base = {
            "chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": rec.alts[0],
            "gene": gene, "hgvsp": hgvsp, "consequence": cons, "tumor_vaf": vaf,
            "tumor_dp": dp, "gnomad_af": gnomad_af, "clinvar_id": clinvar_id,
            "clinvar_sig": clinvar_sig, "kb_tag": args.kb_tag,
        }
        if vid is None:
            if not hgvsp:
                if not cons.startswith("splice") and cons != "transcript_ablation":
                    continue  # no protein change, no KB match, no mechanism note
                rows.append({**base, "civic_variant_id": "", "civic_variant": "",
                             "civic_evidence_levels": "", "civic_evidence_types": "",
                             "civic_therapies": "", "civic_diseases": "",
                             "amp_tier": "", "amp_level": "", "evidence_source": "",
                             "tier_basis": "non_protein_mechanism"})
                continue
            rows.append({**base, "civic_variant_id": "", "civic_variant": "",
                         "civic_evidence_levels": "", "civic_evidence_types": "",
                         "civic_therapies": "", "civic_diseases": "",
                         "amp_tier": "III", "amp_level": "", "evidence_source": "",
                         "tier_basis": "no_kb_match"})
            continue
        v = variants[vid]
        # a benign-only assertion must not erase an actionable tier from the
        # same variant; flip to IV only when nothing actionable was derived
        tier = v["amp_tier"]
        if v["benign"] and not tier:
            tier = "IV"
        rows.append({**base, "civic_variant_id": vid, "civic_variant": v["civic_variant"],
                     "civic_evidence_levels": v["evidence_levels"],
                     "civic_evidence_types": v["evidence_types"],
                     "civic_therapies": v["therapies"], "civic_diseases": v["diseases"],
                     "amp_tier": tier,
                     "amp_level": v["amp_level"],
                     "evidence_source": "CIViC:" + (v["amp_raw"] or v["evidence_levels"]),
                     "tier_basis": basis})

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    sys.stderr.write(f"[clinical_tier] wrote {len(out)} rows -> {args.out}\n")


if __name__ == "__main__":
    main()
