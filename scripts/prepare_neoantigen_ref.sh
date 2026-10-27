#!/usr/bin/env bash
# Neoantigen reference/tools bring-up (venus#31 P0; audit-corrected).
#  1. Verify the SpecHLA conda deployment and locate its bundled IMGT db
#  2. Verify the MHCflurry model cache (models are pre-fetched ONCE with
#     MHCFLURRY_DATA_DIR pointed at the shared dir; see envs/pvac rule shell)
#  3. Verify the Ensembl VEP plugin dir used by the re-annotation rule
#     (Wildtype + Frameshift perl modules — operator-copied, bundled from
#     the Ensembl/VEP_plugins repo at the VEP release tag)
#  4. Audit DTU tools (NetMHCpan/NetMHCIIpan) — academic license forbids
#     redistribution; USED ONLY IF installed manually at DTU_DIR
# Idempotent; prints the [config] lines to paste into venus.oxoflow.
set -euo pipefail

MHCFLURRY_DATA="${MHCFLURRY_DATA:-/home/data/reference/mhcflurry_data}"
VEP_PLUGINS_DIR="${VEP_PLUGINS_DIR:-/home/data/reference/vep_plugins}"
DTU_DIR="${DTU_DIR:-/home/data/reference/dtu_mhc}"

echo "== 1. SpecHLA deployment =="
if ! command -v spechla >/dev/null 2>&1; then
    echo "ERROR: spechla not on PATH — install envs/hla (pixi) and retry" >&2
    exit 1
fi
# bioconda spechla bundles the IMGT db at $CONDA_PREFIX/share/spechla/db and
# its wrapper exports SPECHLA_DB; verify both.
CONDA_PREFIX_HIT="$(cd "$(dirname "$(command -v spechla)")/.." && pwd)"
DB_DIR="$CONDA_PREFIX_HIT/share/spechla/db"
if [ -d "$DB_DIR/HLA" ] && [ -d "$DB_DIR/ref" ]; then
    echo "SpecHLA db OK: $DB_DIR"
    ls "$DB_DIR" | head -3
else
    echo "WARN: expected bundled db missing at $DB_DIR — check the conda pkg"
    echo "      (fallback: clone deepomicslab/SpecHLA and set SPECHLA_DB)"
fi

echo "== 2. MHCflurry model cache =="
mkdir -p "$MHCFLURRY_DATA"
MHCFLURRY_DATA_DIR="$MHCFLURRY_DATA" mhcflurry-downloads fetch models 2>&1 | tail -2 || \
    echo "WARN: model prefetch failed here — run it inside envs/pvac:"
echo "      MHCFLURRY_DATA_DIR=$MHCFLURRY_DATA pixi run --manifest-path envs/pvac/pixi.toml mhcflurry-downloads fetch models"
MHCFLURRY_DATA_DIR="$MHCFLURRY_DATA" mhcflurry-downloads info 2>&1 | head -5 || true

echo "== 3. VEP plugin dir (Wildtype/Frameshift) =="
mkdir -p "$VEP_PLUGINS_DIR"
for pm in Wildtype.pm Frameshift.pm; do
    if [ -f "$VEP_PLUGINS_DIR/$pm" ]; then
        echo "OK  $VEP_PLUGINS_DIR/$pm"
    else
        echo "MISSING $pm — copy it from the Ensembl/VEP_plugins repo (match"
        echo "        your VEP release) into $VEP_PLUGINS_DIR and re-run"
    fi
done

echo "== 4. DTU tools audit =="
printf "tool\tpresent\tpath\n" > /tmp/neo_tool_status.tsv
for t in NetMHCpan NetMHCIIpan; do
    p=""
    [ -x "$DTU_DIR/$t" ] && p="$DTU_DIR/$t"
    command -v "$t" >/dev/null 2>&1 && p="$(command -v $t)"
    printf "%s\t%s\t%s\n" "$t" "$([ -n "$p" ] && echo yes || echo no)" "$p" \
        >> /tmp/neo_tool_status.tsv
done
cat /tmp/neo_tool_status.tsv
if grep -q "no	" /tmp/neo_tool_status.tsv; then
    echo "NOTE: absent DTU tools stay excluded from P0 predictions (class-II"
    echo "      predictions will be reported as not_run, honestly). To enable:"
    echo "      download NetMHCpan-4.1b / NetMHCIIpan (academic license) from"
    echo "      https://services.healthtech.dtu.dk and place executables in $DTU_DIR"
fi

echo "== [config] lines =="
cat << EOF
neo_predictors = "MHCflurry"
neo_binding_nm = 500
neo_min_alleles = 2
neo_tumor_purity = ""
mhcflurry_data_dir = "$MHCFLURRY_DATA"
vep_plugins_dir = "$VEP_PLUGINS_DIR"
EOF
echo "PREPARE_NEOANTIGEN_DONE"
