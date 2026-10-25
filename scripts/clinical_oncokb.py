#!/usr/bin/env python3
"""Optional OncoKB enrichment for venus clinical interpretation (token-gated).

Reads the clinical_tier.py SNV evidence table and queries the OncoKB annotate
API (https://api.oncokb.org/api/v1/variants/annotate) for each matched/
protein-altering variant, appending OncoKB levels (therapeutic 1-4, resistance
R1/R2) when returned. Requires a free academic token (register at
https://www.oncokb.org/api-access with an institutional email).

Failures degrade to empty oncokb_* columns with a stderr note — the CIViC
layer remains the backbone when OncoKB is unavailable.

Output: same table + oncokb_level/oncokb_treatment/oncokb_mutation_effect columns.
"""
import argparse
import json
import os
import sys
import urllib.error
import urllib.request

import pandas as pd

API = "https://api.oncokb.org/api/v1/variants/annotate"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snv-evidence", required=True)
    ap.add_argument("--token", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.snv_evidence, sep="\t", dtype=str).fillna("")
    for col in ("oncokb_level", "oncokb_treatment", "oncokb_mutation_effect"):
        if col not in df.columns:
            df[col] = ""

    todo = df[(df["hgvsp"] != "") | (df["civic_variant_id"] != "")]
    if todo.empty:
        df.to_csv(args.out, sep="\t", index=False)
        sys.exit(0)

    queries = []
    for _, r in todo.iterrows():
        variant = r["hgvsp"] if r["hgvsp"] else str(r["civic_variant"])
        queries.append({"gene": {"hugoSymbol": r["gene"]}, "variant": variant,
                        "type": "variant", "sampleId": "query"})
    body = json.dumps(queries).encode()
    req = urllib.request.Request(API, data=body, headers={
        "Content-Type": "application/json", "Authorization": f"Bearer {args.token}"})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            results = json.load(resp)
    except urllib.error.HTTPError as exc:
        sys.stderr.write(f"[clinical_oncokb] API HTTP {exc.code} — "
                         f"check token validity; leaving columns empty\n")
        df.to_csv(args.out, sep="\t", index=False)
        sys.exit(0)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        sys.stderr.write(f"[clinical_oncokb] API unreachable ({exc}) — "
                         f"leaving columns empty; CIViC layer unaffected\n")
        df.to_csv(args.out, sep="\t", index=False)
        sys.exit(0)

    for (idx, _), res in zip(todo.iterrows(), results):
        if not isinstance(res, dict):
            continue
        oncogenic = res.get("oncogenic", "") or ""
        df.at[idx, "oncokb_mutation_effect"] = oncogenic
        # OncoKB returns treatment evidence inside higher-level payloads;
        # best-effort extraction of the strongest therapeutic level.
        lvl = res.get("highestLevelOfEvidence") or ""
        treatments = []
        for t in (res.get("treatments") or []):
            treatments.append(f"{t.get('drugs','') or t.get('drug','')}"
                              f"[{t.get('level','')}]")
        df.at[idx, "oncokb_level"] = str(lvl)
        df.at[idx, "oncokb_treatment"] = ";".join(x for x in treatments if x.strip("[]"))

    df.to_csv(args.out, sep="\t", index=False)
    sys.stderr.write(f"[clinical_oncokb] enriched {len(todo)} variants -> {args.out}\n")


if __name__ == "__main__":
    main()
