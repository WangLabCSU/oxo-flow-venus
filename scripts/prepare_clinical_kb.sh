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

print_config() {
  echo "[prepare_clinical_kb] done. Contents of $KB_DIR:"
  ls -la "$KB_DIR"
  echo
  echo "Update venus.oxoflow [config]:"
  echo "  clinical_kb_tag = \"CIViC_$TAG\""
  echo "  civic_variant_tsv   = \"$KB_DIR/${TAG}-VariantSummaries.tsv\""
  echo "  civic_evidence_tsv  = \"$KB_DIR/${TAG}-AcceptedClinicalEvidenceSummaries.tsv\""
  echo "  civic_assertion_tsv = \"$KB_DIR/${TAG}-AcceptedAssertionSummaries.tsv\""
  echo "  clinvar_vcf         = \"$KB_DIR/clinvar.chr.vcf.gz\""
}

mkdir -p "$KB_DIR"
cd "$KB_DIR"

# Fast path: everything already deployed for the CURRENT tag — skip discovery
# and downloads entirely (also makes offline re-runs and cron no-ops safe).
CIVIC_TAG_FILE="$KB_DIR/.civic_tag"
if [ -s "$CIVIC_TAG_FILE" ] && grep -q . "$CIVIC_TAG_FILE"; then
  EXISTING_TAG=$(cat "$CIVIC_TAG_FILE")
  if [ -s "${EXISTING_TAG}-VariantSummaries.tsv" ] \
     && [ -s "${EXISTING_TAG}-AcceptedClinicalEvidenceSummaries.tsv" ] \
     && [ -s "${EXISTING_TAG}-AcceptedAssertionSummaries.tsv" ] \
     && [ -s clinvar.vcf.gz ] && [ -s clinvar.chr.vcf.gz ]; then
    TAG="$EXISTING_TAG"
    echo "[prepare_clinical_kb] KB already complete for tag $TAG — nothing to do"
    print_config
    exit 0
  fi
fi

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

# Normalize ClinVar contigs to chr-prefixed names (NCBI ships "1"; the
# pipeline's GRCh38.d1.vd1 uses "chr1" — without this the bcftools join is
# silently empty). One-time rename + index into clinvar.chr.vcf.gz.
if [ ! -s clinvar.chr.vcf.gz ]; then
  echo "[prepare_clinical_kb] renaming ClinVar contigs to chr-prefix..."
  printf '1 chr1\n2 chr2\n3 chr3\n4 chr4\n5 chr5\n6 chr6\n7 chr7\n8 chr8\n9 chr9\n10 chr10\n11 chr11\n12 chr12\n13 chr13\n14 chr14\n15 chr15\n16 chr16\n17 chr17\n18 chr18\n19 chr19\n20 chr20\n21 chr21\n22 chr22\nX chrX\nY chrY\nMT chrM\n' > /tmp/clinvar_chr_map.txt
  bcftools annotate --rename-chrs /tmp/clinvar_chr_map.txt clinvar.vcf.gz -Oz -o clinvar.chr.vcf.gz
  bcftools index -f -t clinvar.chr.vcf.gz
fi

# Record which CIViC tag this directory is complete for (fast-path key).
echo -n "$TAG" > "$CIVIC_TAG_FILE"
print_config

