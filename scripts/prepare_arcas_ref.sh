#!/usr/bin/env bash
# One-time arcasHLA IMGT reference bring-up (venus#31 P1 item 3).
#
# WHY THIS EXISTS
# arcasHLA 0.6.0 (bioconda hdfd78af_2, verified from the actual package) does
# NOT ship its dat/ reference; `arcasHLA reference` populates
#   <env-prefix>/dat/IMGTHLA  (git clone of https://github.com/ANHIG/IMGTHLA)
#   <env-prefix>/dat/ref/{hla.p.json,hla.idx,hla.fa}
# INSIDE the pixi environment prefix. That means the reference build mutates
# the pixi env directory and must run ONCE per machine, after `pixi install`
# of envs/arcas, BEFORE the hla_type_rna rule runs. It is an ops step, not a
# rule step: a rule cannot depend on its own environment's state.
#
# USAGE (on lab-wsx-in, inside a checked-out venus workdir)
#   cd <workdir> && pixi run --environment arcas bash scripts/prepare_arcas_ref.sh [IMGT_VERSION]
#
# ARGUMENT
#   IMGT_VERSION  arcasHLA `reference --version` name, e.g. 3.46.0 (default).
#                 Must be a name in arcasHLA's shipped parameters.json versions
#                 map (3.3.0..3.46.0; verified map in the 0.6.0 package).
#                 Newer-than-map IMGT releases need `--commit <hash>` instead —
#                 edit DEFAULT below only after checking the map.
#
# NOTES
#   - IMGT hla.dat is stored with Git LFS for versions >= 3.35.0; the arcas
#     env pulls in git + git-lfs for exactly this step (verified dep list).
#   - Re-runnable: `arcasHLA reference` skips completed stages; pass --rebuild
#     via env var REBUILD=1 to force a fresh build.
set -euo pipefail

DEFAULT_VERSION="3.46.0"
VERSION="${1:-$DEFAULT_VERSION}"

command -v arcasHLA >/dev/null 2>&1 || {
  echo "[prepare_arcas_ref] arcasHLA not on PATH — run inside the arcas env" >&2
  exit 1
}

ARCASHLA_ROOT="${ARCASHLA_ROOT_DIR:-$(realpath "$(dirname "$(command -v arcasHLA)")/..")}"
DAT_DIR="${ARCASHLA_ROOT}/dat"
echo "[prepare_arcas_ref] env root : ${ARCASHLA_ROOT}"
echo "[prepare_arcas_ref] dat dir  : ${DAT_DIR}"
echo "[prepare_arcas_ref] version  : ${VERSION}"

if [ "${REBUILD:-0}" = "1" ] && [ -d "${DAT_DIR}/IMGTHLA" ]; then
  echo "[prepare_arcas_ref] REBUILD=1 — removing existing IMGTHLA checkout"
  rm -rf "${DAT_DIR}/IMGTHLA"
fi

# Idempotent clone stage (arcasHLA reference re-runs cleanly; git-lfs pull is
# handled inside its check_ref()/build path for >=3.35.0 versions).
if [ ! -d "${DAT_DIR}/IMGTHLA" ]; then
  echo "[prepare_arcas_ref] cloning ANHIG/IMGTHLA (git-lfs enabled)"
  mkdir -p "${DAT_DIR}"
  git lfs install --local >/dev/null 2>&1 || true
  git clone https://github.com/ANHIG/IMGTHLA.git "${DAT_DIR}/IMGTHLA"
fi

echo "[prepare_arcas_ref] building dat/ref via arcasHLA reference --version ${VERSION}"
if [ "${REBUILD:-0}" = "1" ]; then
  arcasHLA reference --rebuild --version "${VERSION}"
else
  arcasHLA reference --version "${VERSION}"
fi

# Post-condition checks: the files genotype.py actually reads.
for f in "${DAT_DIR}/ref/hla.p.json" "${DAT_DIR}/ref/hla.idx" "${DAT_DIR}/ref/hla.fa"; do
  [ -s "$f" ] || { echo "[prepare_arcas_ref] MISSING: $f" >&2; exit 1; }
done

echo "[prepare_arcas_ref] OK — dat/ref built for IMGT ${VERSION}"
