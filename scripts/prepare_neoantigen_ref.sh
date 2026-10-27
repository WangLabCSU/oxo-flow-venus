#!/usr/bin/env bash
# Neoantigen reference/tools bring-up (venus#31 P0).
#  1. Verify the SpecHLA deployment (conda package + bundled IMGT db index)
#  2. Prefetch MHCflurry model data into a shared dir (open license, safe to
#     cache) so rules never hit the network at run time
#  3. Audit DTU tools (NetMHCpan/NetMHCIIpan) — academic license forbids
#     redistribution; they are USED ONLY IF the operator installed them
#     manually at the documented location
# Idempotent; prints the [config] lines to paste into venus.oxoflow.
set -euo pipefail

SPECHLA_REF="${SPECHLA_REF:-/home/data/reference/spechla}"
MHCFLURRY_DATA="${MHCFLURRY_DATA:-/home/data/reference/mhcflurry_data}"
DTU_DIR="${DTU_DIR:-/home/data/reference/dtu_mhc}"
export MHCFLURRY_DATA_DIR="$MHCFLURRY_DATA"

echo "== 1. SpecHLA deployment =="
if ! command -v spechla >/dev/null 2>&1; then
    echo "ERROR: spechla not on PATH — install envs/hla (pixi) and retry" >&2
    exit 1
fi
spechla --version 2>/dev/null || true
# The conda package carries the IMGT db; locate it and verify the index.
SPECHLA_DB="${SPECHLA_DB:-$(python - << 'PY' 2>/dev/null || true
import glob, os
for pat in (os.path.expanduser("~/.pixi/envs/hla/share/spechla*"),
            "/opt/conda/share/spechla*", os.environ.get("CONDA_PREFIX", "") + "/share/spechla*"):
    for d in sorted(glob.glob(pat)):
        if os.path.isdir(d):
            print(d)
            raise SystemExit
PY
)}"
echo "SPECHLA_DB=${SPECHLA_DB:-<unset — rules pass it via config>}"
if [ -n "${SPECHLA_DB:-}" ] && [ -d "$SPECHLA_DB" ]; then
    ls "$SPECHLA_DB" | head -5
else
    echo "WARN: SPECHLA_DB not auto-located; clone deepomicslab/SpecHLA (pinned) "
    echo "      to $SPECHLA_REF and index it (bash index.sh) — then set config.spechla_db"
fi

echo "== 2. MHCflurry models =="
mkdir -p "$MHCFLURRY_DATA"
pip install -q mhcflurry 2>/dev/null || true
mhcflurry-downloads fetch models 2>&1 | tail -2
mhcflurry-downloads info 2>&1 | head -5

echo "== 3. DTU tools audit =="
: > /tmp/neo_tool_status.tsv
printf "tool\tpresent\tpath\n" > /tmp/neo_tool_status.tsv
for t in NetMHCpan NetMHCIIpan; do
    p=""
    [ -x "$DTU_DIR/$t" ] && p="$DTU_DIR/$t"
    command -v "$t" >/dev/null 2>&1 && p="$(command -v $t)"
    printf "%s\t%s\t%s\n" "$t" "$([ -n "$p" ] && echo yes || echo no)" "$p" \
        >> /tmp/neo_tool_status.tsv
done
cat /tmp/neo_tool_status.tsv
if grep -qP "no\t$" /tmp/neo_tool_status.tsv; then
    echo "NOTE: absent DTU tools stay excluded from P0 predictions (class-II"
    echo "      predictions will be reported as not_run, honestly). To enable:"
    echo "      download NetMHCpan-4.1b / NetMHCIIpan (academic license) from"
    echo "      https://services.healthtech.dtu.dk and place the executables in $DTU_DIR"
fi

echo "== [config] lines =="
cat << EOF
neo_predictors = "MHCflurry"
neo_binding_nm = 500
neo_tumor_purity = ""
neo_tumor_ploidy = ""
EOF
echo "PREPARE_NEOANTIGEN_DONE"
