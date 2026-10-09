#!/bin/bash
# Download + verify the official ASCAT G1000 hg38 WGS reference set
# (Zenodo record 14008443) used by rules/ascat.oxoflow.
#
#   G1000_alleles_WGS_hg38.zip    10,173,770 B   e8aebc542222c3e58f1168ba47b299fa
#   G1000_loci_WGS_hg38.zip        8,463,901 B   80c04afd5c7c75c0593591d6be22fa6f
#   GC_G1000_WGS_hg38.zip        104,766,910 B   4e5ab5437aab7dbe097592b0c7a025fb
#   RT_G1000_WGS_hg38.zip        106,368,554 B   6ebda3b1abe2879adaee5fd16cea3fe4
#
# NOTE on names: the GC/RT zips ship INNER files named GC_G1000_hg38.txt /
# RT_G1000_hg38.txt (no "WGS"); this script extracts and copies them to the
# *_WGS_* names that rules/ascat.oxoflow passes via --gc/--rt.
#
# Every zip is size-checked AND md5-checked before extraction, and the unpacked
# text tables are validated afterwards (22 alleles + 22 loci files; GC/RT are
# tab-separated with a header, every row 18 fields; RT spans chr1..X, so chr22
# presence is asserted explicitly). Background: a corrupted RESUMED download of
# the RT zip (12 MB of appended garbage that still passed a naive size glance)
# killed ASCAT mid-run with "Stopped early on line 1813274. Expected 18 fields
# but found 51". The md5 gate below makes that failure mode impossible to
# reach the refdir.
#
# Usage:
#   scripts/download_ascat_refdata.sh -o /path/to/refdir [-r]
#     -o  target directory (what config.ascat_refdir points at)
#     -r  resume mode for flaky links: curl -C - in a loop until complete.
#         Safe here because the md5 gate runs on the final file either way;
#         curl also refuses to resume against a server that ignores Range.
set -euo pipefail

RESUME=0
OUT=""
while getopts "ro:" opt; do
  case $opt in
    r) RESUME=1 ;;
    o) OUT=$OPTARG ;;
    *) echo "usage: $0 -o <refdir> [-r]" >&2; exit 64 ;;
  esac
done
[ -n "$OUT" ] || { echo "usage: $0 -o <refdir> [-r]" >&2; exit 64; }
mkdir -p "$OUT"

RECORD="https://zenodo.org/records/14008443/files"

fetch() { # <name> <expected-size> <expected-md5>
  local name=$1 size=$2 md5=$3
  local zip="$OUT/$name"
  if [ -f "$zip" ] \
     && [ "$(stat -c %s "$zip")" = "$size" ] \
     && [ "$(md5sum "$zip" | cut -d' ' -f1)" = "$md5" ]; then
    echo "[skip] $name already present and verified"
  else
    rm -f "$zip"
    if [ "$RESUME" = 1 ]; then
      local i
      for i in $(seq 1 60); do
        [ "$(stat -c %s "$zip" 2>/dev/null || echo 0)" = "$size" ] && break
        echo "[fetch $i] $name ($(stat -c %s "$zip" 2>/dev/null || echo 0)/$size bytes)"
        curl -fsSL -C - --speed-time 30 --speed-limit 1024 \
          -o "$zip" "$RECORD/$name?download=1" || true
        sleep 3
      done
    else
      curl -fsSL --retry 5 --retry-delay 10 -o "$zip" "$RECORD/$name?download=1"
    fi
    local got_sz got_md5
    got_sz=$(stat -c %s "$zip")
    got_md5=$(md5sum "$zip" | cut -d' ' -f1)
    [ "$got_sz" = "$size" ] || { echo "FATAL: $name size $got_sz != $size" >&2; exit 1; }
    [ "$got_md5" = "$md5" ] || { echo "FATAL: $name md5 $got_md5 != $md5" >&2; exit 1; }
    unzip -t "$zip" > /dev/null || { echo "FATAL: $name fails unzip -t" >&2; exit 1; }
    echo "[ok] $name verified, extracting"
    unzip -o -q "$zip" -d "$OUT"
  fi
}

fetch G1000_alleles_WGS_hg38.zip  10173770 e8aebc542222c3e58f1168ba47b299fa
fetch G1000_loci_WGS_hg38.zip      8463901 80c04afd5c7c75c0593591d6be22fa6f
fetch GC_G1000_WGS_hg38.zip      104766910 4e5ab5437aab7dbe097592b0c7a025fb
fetch RT_G1000_WGS_hg38.zip      106368554 6ebda3b1abe2879adaee5fd16cea3fe4

# --- rule-expected names: GC/RT zips ship *_hg38.txt, the rules consume
# --- *_WGS_hg38.txt. Copy (don't move) so re-extraction stays idempotent.
cp -f "$OUT/GC_G1000_hg38.txt" "$OUT/GC_G1000_WGS_hg38.txt"
cp -f "$OUT/RT_G1000_hg38.txt" "$OUT/RT_G1000_WGS_hg38.txt"

# --- post-extraction validation of the text tables --------------------------
n_alleles=$(ls "$OUT"/G1000_alleles_hg38_*.txt 2>/dev/null | wc -l)
n_loci=$(ls "$OUT"/G1000_loci_hg38_*.txt 2>/dev/null | wc -l)
[ "$n_alleles" -eq 22 ] || { echo "FATAL: expected 22 alleles files, got $n_alleles" >&2; exit 1; }
[ "$n_loci" -eq 22 ]    || { echo "FATAL: expected 22 loci files, got $n_loci" >&2; exit 1; }

for t in GC RT; do
  bad=$(awk -F '\t' 'NF!=18' "$OUT/${t}_G1000_WGS_hg38.txt" | wc -l)
  [ "$bad" -eq 0 ] || { echo "FATAL: $t table has $bad malformed rows (expect 18 tab fields)" >&2; exit 1; }
done
# RT spans chr1..X; make sure chr22 is present (truncation sentinel)
awk -F '\t' '$2=="22"{found=1} END{exit !found}' "$OUT/RT_G1000_WGS_hg38.txt" \
  || { echo "FATAL: RT table does not cover chr22 (truncated download?)" >&2; exit 1; }

echo "ALL OK: ASCAT hg38 WGS reference set verified in $OUT"
