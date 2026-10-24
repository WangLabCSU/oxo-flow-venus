#!/usr/bin/env bash
# Host-side knowledgebase preparation for the venus clinical module.
#
# Downloads the LATEST CIViC monthly dump (variant/evidence/assertion TSVs)
# and the ClinVar weekly GRCh38 VCF into /home/data/reference/clinical_kb/
# (override with CLINICAL_KB_DIR). Not a pipeline rule: KB refreshes are
# deliberate host operations so rule fingerprints stay stable. After running,
# update the KB paths/tag in venus.oxoflow [config] and docs/refdata.md.
#
# CIViC release file names are discovered via the public GraphQL API
# (dataReleases { name variantTsv { path } ... }); ClinVar comes from the NCBI
# FTP. Re-running is idempotent (skips complete files).
set -euo pipefail

KB_DIR="${CLINICAL_KB_DIR:-/home/data/reference/clinical_kb}"
API="https://civicdb.org/api/graphql"
CLINVAR_BASE="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"

mkdir -p "$KB_DIR"
cd "$KB_DIR"

echo "[prepare_clinical_kb] discovering latest CIViC monthly release..."
CIVIC_JSON=$(curl -sSL --retry 3 "$API" -H 'Content-Type: application/json' \
  -d '{"query":"{ dataReleases { name variantTsv { path } acceptedEvidenceTsv { path } acceptedAssertionTsv { path } } }"}')
# dataReleases[0] is nightly; pick the first non-nightly (most recent monthly).
RELEASE=$(python3 - "$CIVIC_JSON" <<'PY'
import json, sys
d = json.load(sys.stdin)
rel = [r for r in d["data"]["dataReleases"] if r["name"] != "nightly"][0]
print(rel["name"])
for k in ("variantTsv", "acceptedEvidenceTsv", "acceptedAssertionTsv"):
    print("https://civicdb.org" + rel[k]["path"])
PY
)
mapfile -t URLS <<< "$RELEASE"
TAG="${URLS[0]}"
echo "[prepare_clinical_kb] CIViC release: $TAG"
i=0
for kind in variant evidence assertion; do
  i=$((i+1))
  url="${URLS[$i]}"
  fname="$(basename "$url")"
  if [ ! -s "$fname" ]; then
    echo "  downloading $fname"
    curl -sSL --retry 3 -o "$fname" "$url"
  else
    echo "  present: $fname"
  fi
done

for f in clinvar.vcf.gz clinvar.vcf.gz.tbi; do
  if [ ! -s "$f" ]; then
    echo "[prepare_clinical_kb] downloading $f"
    curl -sSL --retry 3 -o "$f" "$CLINVAR_BASE/$f"
  else
    echo "[prepare_clinical_kb] present: $f"
  fi
done

echo "[prepare_clinical_kb] done. Contents of $KB_DIR:"
ls -la
echo
echo "Update venus.oxoflow [config]:"
echo "  clinical_kb_tag = \"CIViC_$TAG\""
echo "  civic_variant_tsv   = \"$KB_DIR/${TAG}-VariantSummaries.tsv\""
echo "  civic_evidence_tsv  = \"$KB_DIR/${TAG}-AcceptedClinicalEvidenceSummaries.tsv\""
echo "  civic_assertion_tsv = \"$KB_DIR/${TAG}-AcceptedAssertionSummaries.tsv\""
echo "  clinvar_vcf         = \"$KB_DIR/clinvar.vcf.gz\""
