#!/usr/bin/env bash
# Prepare the TRUST4 reference files (venus#31 P1 item 2, tcr_from_rna).
#
# TRUST4 needs two small reference files (both <350 KB):
#   1. -f  : V/D/J/C gene sequences WITH coordinates. For BAM input this must
#            be built from the SAME genome + GTF as the alignment (the
#            repo-tracked hg38_bcrtcr.fa is Gencode v35 + hg38; our genome is
#            GRCh38.d1.vd1 + GENCODE v47, so we rebuild for coordinate
#            consistency).
#   2. --ref  the IMGT V/D/J/C reference with annotation (human_IMGT+C.fa).
#
# This script downloads the IMGT reference from the TRUST4 GitHub repo (the
# exact file the authors ship for hg38) and builds the coordinate fasta from
# OUR genome + GTF via the package's BuildDatabaseFa.pl, using the IMGT gene
# list. Outputs land in $TRUST_REF_DIR (venus.oxoflow config key):
#   bcrtcr.fa        (the -f file)
#   human_IMGT+C.fa  (the --ref file)
#
# Usage: TRUST_REF_DIR=/home/data/reference/trust4 bash scripts/prepare_trust4_ref.sh \
#            /home/data/reference/GRCh38.d1.vd1.fa \
#            /home/data/reference/gencode_v47/gencode.v47.primary_assembly.annotation.gtf
set -euo pipefail

GENOME_FA="${1:?usage: prepare_trust4_ref.sh <genome.fa> <gencode.gtf>}"
GENCODE_GTF="${2:?usage: prepare_trust4_ref.sh <genome.fa> <gencode.gtf>}"
TRUST_REF_DIR="${TRUST_REF_DIR:?set TRUST_REF_DIR (venus.oxoflow: trust_ref_dir)}"

REPO_URL="https://api.github.com/repos/liulab-dfci/TRUST4/contents"
PIXI_TRUST4="${PIXI_TRUST4:?set PIXI_TRUST4 to the pixi env run dir (trust4 1.1.10, provides perl + BuildImgtAnnot.pl)}"

mkdir -p "$TRUST_REF_DIR"
cd "$PIXI_TRUST4"

echo "[prepare_trust4_ref] 1/3 human_vdjc.list from IMGT+C.fa gene headers"
curl -fsSL "$REPO_URL/human_IMGT%2BC.fa" -H "Accept: application/vnd.github.raw" \
    -o "$TRUST_REF_DIR/human_IMGT+C.fa"
grep ">" "$TRUST_REF_DIR/human_IMGT+C.fa" | cut -f2 -d'>' | cut -f1 -d'*' | sort | uniq \
    > "$TRUST_REF_DIR/bcr_tcr_gene_name.txt"

echo "[prepare_trust4_ref] 2/3 bcrtcr.fa from OUR genome + GTF (BuildDatabaseFa.pl)"
perl BuildDatabaseFa.pl "$GENOME_FA" "$GENCODE_GTF" \
    "$TRUST_REF_DIR/bcr_tcr_gene_name.txt" > "$TRUST_REF_DIR/bcrtcr.fa"

echo "[prepare_trust4_ref] done -> $TRUST_REF_DIR/{bcrtcr.fa,human_IMGT+C.fa}"
ls -la "$TRUST_REF_DIR"
