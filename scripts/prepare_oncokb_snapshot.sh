#!/usr/bin/env bash
# OncoKB snapshot via the academic REST API (the only official data channel;
# OncoKB content must NOT be committed into any repository — snapshots live
# under the server reference dir, versioned by date, refreshed monthly).
#   ONCOKB_TOKEN_FILE  file containing the bearer token (required for the
#                      production host; the demo host needs none)
#   ONCOKB_BASE        API base (default https://api.oncokb.org/api/v1;
#                      https://demo.oncokb.org/api/v1 for connectivity checks)
#   ONCOKB_OUT         snapshot dir (default /home/data/reference/oncokb/<date>)
# Produces: genes.json, drugs.json, alterations per gene, a flat TSV view and
# an md5 manifest; prints the [config] lines for venus.oxoflow.
set -euo pipefail

BASE="${ONCOKB_BASE:-https://api.oncokb.org/api/v1}"
OUT="${ONCOKB_OUT:-/home/data/reference/oncokb/$(date +%Y%m%d)}"
TOKEN_FILE="${ONCOKB_TOKEN_FILE:-}"
AUTH=()
if [ -n "$TOKEN_FILE" ] && [ -s "$TOKEN_FILE" ]; then
    AUTH=(-H "Authorization: Bearer $(cat "$TOKEN_FILE")")
elif [ "$BASE" != "https://demo.oncokb.org/api/v1" ]; then
    echo "ERROR: production host requires ONCOKB_TOKEN_FILE" >&2
    exit 1
fi

mkdir -p "$OUT"
get() { curl -sL --max-time 120 --retry 3 --retry-delay 10 ${AUTH[@]+"${AUTH[@]}"} "$BASE/$1" -o "$OUT/$2"; }

echo "== /genes =="
get genes genes.json
gene_ids=$(python3 - << PY
import json
genes = json.load(open("$OUT/genes.json"))
print(" ".join(str(g["entrezGeneId"]) for g in genes))
PY
)
n=$(echo $gene_ids | wc -w)
echo "genes: $n"

echo "== per-gene alterations + evidence =="
: > "$OUT/alterations.tsv"
printf "entrezGeneId\thugoSymbol\talteration\toncogenicity\tdescription\n" > "$OUT/alterations.tsv"
for gid in $gene_ids; do
    get "genes/$gid/alterations" "alt_$gid.json" || continue
    python3 - "$OUT/alt_$gid.json" "$gid" >> "$OUT/alterations.tsv" << 'PY'
import json, sys, csv
gid = sys.argv[2]
try:
    alts = json.load(open(sys.argv[1]))
except Exception:
    sys.exit(0)
for a in alts if isinstance(alts, list) else []:
    csv.writer(sys.stdout, delimiter="\t", lineterminator="\n").writerow(
        [gid, a.get("gene", {}).get("hugoSymbol", ""), a.get("name", ""),
         a.get("oncogenic", ""), (a.get("description", "") or "").replace("\t", " ")])
PY
done

echo "== /drugs =="
get drugs drugs.json

echo "== manifest =="
( cd "$OUT" && md5sum genes.json drugs.json alterations.tsv alt_*.json > md5_manifest.txt 2>/dev/null || true )
echo "== [config] lines =="
cat << EOF
oncokb_snapshot_dir = "$OUT"
EOF
echo "ONCOKB_SNAPSHOT_DONE ($OUT)"
