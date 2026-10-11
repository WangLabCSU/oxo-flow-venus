#!/bin/bash
# Verify every host reference path declared in venus.oxoflow [config].
#
#   scripts/verify_refdata.sh [path/to/venus.oxoflow]
#
# Checks existence/type of each configured path (genome FASTA + BWA-MEM2
# index, STAR index, GENCODE GTF, GATK resources, vep-rs binary + cache,
# dnbc4tools reference/binary, ASCAT tables) and re-runs the ASCAT table
# validation from scripts/download_ascat_refdata.sh without re-downloading.
#
# Layout note (2026-10-11): all assets physically live under the canonical
# root /home/data/reference/; a few config strings still point at historical
# locations that are now compat symlinks into that root. This script checks
# the CONFIGURED strings, so symlink breakage is exactly what it catches.
# cnv_reference is WARN-only: it is documented as configured-but-not-deployed
# (docs/refdata.md §7) and only tumor-only CNV would need it.
#
# Linux-only (stat -c / readlink -f); run on the reference host.
set -euo pipefail

WF="${1:-venus.oxoflow}"
[ -f "$WF" ] || { echo "FATAL: workflow file not found: $WF" >&2; exit 1; }
WF=$(readlink -f "$WF")

FAIL=0
ok()   { echo "[ok]   $*"; }
bad()  { echo "[FAIL] $*"; FAIL=1; }
warn() { echo "[WARN] $*"; }
note() { echo "[note] $*"; }

cfg() { # cfg <key> — string value of a [config] key in $WF (keys are unique)
  sed -n 's/^[[:space:]]*'"$1"'[[:space:]]*=[[:space:]]*"\(.*\)".*/\1/p' "$WF"
}

# check <label> <kind: f|d|x|ne> <path>  (ne = non-empty dir)
check() {
  local label=$1 kind=$2 path=$3
  if [ ! -e "$path" ]; then bad "$label missing: $path"; return 0; fi
  case $kind in
    f)  [ -f "$path" ] || { bad "$label not a regular file: $path"; return 0; } ;;
    d)  [ -d "$path" ] || { bad "$label not a directory: $path"; return 0; } ;;
    x)  [ -x "$path" ] || { bad "$label not executable: $path"; return 0; } ;;
    ne) [ -d "$path" ] && [ -n "$(ls -A "$path" 2>/dev/null)" ] \
          || { bad "$label empty or not a directory: $path"; return 0; } ;;
  esac
  local phys
  phys=$(readlink -f "$path")
  if [ "$phys" != "$path" ]; then
    note "$label: $path -> $phys (compat symlink)"
  fi
  ok "$label"
}

echo "== verify_refdata: $WF =="

# --- genome + aligner index -------------------------------------------------
GENOME_FA=$(cfg genome_fa);        GENOME_FAI=$(cfg genome_fai)
GENOME_DICT=$(cfg genome_dict);    BWA_PREFIX=$(cfg bwa_prefix)
check genome_fa f "$GENOME_FA"
check genome_fai f "$GENOME_FAI"
check genome_dict f "$GENOME_DICT"
for sfx in 0123 amb ann bwt.2bit.64 pac; do
  check "bwa_index.$sfx" f "$BWA_PREFIX.$sfx"
done

# --- bulk RNA ---------------------------------------------------------------
STAR_INDEX=$(cfg star_index); GENCODE_GTF=$(cfg gencode_gtf)
check star_index d "$STAR_INDEX"
check "star_index.genomeParameters" f "$STAR_INDEX/genomeParameters.txt"
check gencode_gtf f "$GENCODE_GTF"

# --- GATK resources ---------------------------------------------------------
for key in dbsnp mills gnomad pon; do
  p=$(cfg "$key")
  check "$key" f "$p"
  check "$key.tbi" f "$p.tbi"
done

# --- VEP --------------------------------------------------------------------
VEP_BIN=$(cfg vep_bin); VEP_CACHE=$(cfg vep_cache_dir)
check vep_bin x "$VEP_BIN"
check vep_cache_dir d "$VEP_CACHE"
check "vep_cache.info.json" f "$VEP_CACHE/info.json"

# --- single cell (dnbc4tools) ----------------------------------------------
SCRNA_REF=$(cfg scrna_reference); SCRNA_BIN=$(cfg scrna_bin)
check scrna_reference ne "$SCRNA_REF"
check scrna_bin x "$SCRNA_BIN"

# --- CNVkit (configured, NOT deployed — docs/refdata.md §7) -----------------
CNV_REF=$(cfg cnv_reference)
if [ -e "$CNV_REF" ]; then
  check cnv_reference f "$CNV_REF"
else
  warn "cnv_reference not deployed: $CNV_REF (expected; only tumor-only CNV needs it)"
fi

# --- ASCAT (Zenodo 14008443, WGS hg38) --------------------------------------
ASDIR=$(cfg ascat_refdir)
check ascat_refdir d "$ASDIR"
if [ -d "$ASDIR" ]; then
  # run_ascat.R probes the PLAIN-numbered names first — that is the contract.
  n_all=$(ls "$ASDIR"/G1000_alleles_hg38_[0-9]*.txt 2>/dev/null | wc -l || true)
  n_loc=$(ls "$ASDIR"/G1000_loci_hg38_[0-9]*.txt 2>/dev/null | wc -l || true)
  [ "$n_all" -eq 22 ] || bad "ASCAT plain alleles files: $n_all (expect 22)"
  [ "$n_loc" -eq 22 ] || bad "ASCAT plain loci files: $n_loc (expect 22)"
  [ "$FAIL" -eq 0 ] && ok "ASCAT plain allele/loci sets (22+22)"
  # chr-prefixed copies are run_ascat.R's fallback naming — report, don't gate.
  n_allc=$(ls "$ASDIR"/G1000_alleles_hg38_chr*.txt 2>/dev/null | wc -l || true)
  n_locc=$(ls "$ASDIR"/G1000_loci_hg38_chr*.txt 2>/dev/null | wc -l || true)
  note "ASCAT chr-prefixed fallback sets: alleles=$n_allc loci=$n_locc (optional)"
  for t in GC RT; do
    f="$ASDIR/${t}_G1000_WGS_hg38.txt"
    if [ -f "$f" ]; then
      badrows=$(awk -F '\t' 'NF!=18' "$f" | wc -l || true)
      [ "$badrows" -eq 0 ] || bad "$t table has $badrows malformed rows (expect 18 tab fields)"
      [ "$FAIL" -eq 0 ] && ok "ASCAT ${t}_G1000_WGS_hg38.txt (18-field)"
    else
      bad "missing $f"
    fi
  done
  if [ -f "$ASDIR/RT_G1000_WGS_hg38.txt" ]; then
    if awk -F '\t' '$2=="22"{f=1} END{exit !f}' "$ASDIR/RT_G1000_WGS_hg38.txt"; then
      ok "ASCAT RT table covers chr22 (truncation sentinel)"
    else
      bad "ASCAT RT table does not cover chr22 (truncated?)"
    fi
  fi
fi

# --- pipeline helpers -------------------------------------------------------
VT=$(cfg venus_tools)
[ -d "$VT" ] && ok "venus_tools dir: $VT" || note "venus_tools relative: $VT (resolve from the run cwd)"

echo "== done =="
if [ "$FAIL" -eq 0 ]; then
  echo "ALL REFERENCE PATHS OK"
else
  echo "REFERENCE VERIFICATION FAILED" >&2
  exit 1
fi
