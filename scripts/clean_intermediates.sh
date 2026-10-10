#!/usr/bin/env bash
# One-click intermediate cleanup for the venus cohort workdir.
#
# Two modes (see docs/cleanup.md for the full rationale):
#   --mode all    clean EVERY intermediate that is not protected_output
#                 (engine `oxo-flow clean --force`). Delivery/report/ASCAT/CNV/
#                 MSI products are protected in rules/*.oxoflow and survive.
#   --mode useless  clean only "useless" scratch (engine temporary=true
#                 tombstones). NOT USABLE YET: tombstone cleanup silently
#                 skips wildcard-rule outputs on current oxo-flow engines
#                 (all venus scratch rules use {wildcard} outputs). This
#                 script refuses to rm manually — manual deletion without
#                 checkpoint tombstones breaks incremental bookkeeping.
#
# Usage:
#   scripts/clean_intermediates.sh                 # preview (dry-run) of --mode all
#   scripts/clean_intermediates.sh --mode all --apply
#   scripts/clean_intermediates.sh --mode useless  # status only, deletes nothing
#   scripts/clean_intermediates.sh --mode all --apply --no-backup
set -euo pipefail

cd "$(dirname "$0")/.."
MODE="all"
APPLY=0
BACKUP=1
while [ $# -gt 0 ]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --apply) APPLY=1; shift ;;
        --no-backup) BACKUP=0; shift ;;
        -h|--help) sed -n '2,19p' "$0"; exit 0 ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

command -v oxo-flow >/dev/null || {
    echo "oxo-flow not on PATH (server: export PATH=/home/wsx/.cargo/bin:\$PATH)" >&2
    exit 1
}

if [ "$MODE" = "useless" ]; then
    cat <<'EOF'
--mode useless relies on engine `temporary = true` tombstone cleanup, which
currently SILENTLY NO-OPS for wildcard rules (every venus scratch rule has
{wildcard} outputs). Tracked upstream (Traitome/oxo-flow#844); until it lands, use `--mode all` to
reclaim scratch (delivery/report/ASCAT/CNV products stay protected), or wait.

Nothing was deleted.
EOF
    oxo-flow clean venus.oxoflow
    exit 0
fi

STAMP="$(date +%Y%m%d-%H%M%S)"
if [ "$APPLY" = "1" ] && [ "$BACKUP" = "1" ]; then
    mkdir -p backups
    echo "== backup: deliver/ report/ ascat/ msi/ + undeclared evidence -> backups/cleanup-${STAMP}.tar.gz"
    # Belt-and-braces: deliver/*.msi.tsv and msi/reference.list are already
    # protected, but msi/ also holds the undeclared {pair_id}_all/_dis/_unstable
    # side products (the whole msi/ dir is archived, so all three land in the tar).
    # The extra paths are undeclared evidence-grade side products that `clean`
    # never sees: germline/candidate Manta VCFs, scRNA per-cell results
    # (QC_Cluster.h5ad, singlecell.csv, report HTMLs) and STAR junctions.
    tar czf "backups/cleanup-${STAMP}.tar.gz" --ignore-failed-read \
        deliver report ascat msi \
        manta/*/results/variants \
        scrna/count/*/outs/analysis \
        scrna/count/*/outs/singlecell.csv \
        scrna/count/*/outs/*_scRNA_report.html \
        rna/star/*/*.SJ.out.tab
    du -h "backups/cleanup-${STAMP}.tar.gz"
fi

if [ "$APPLY" = "1" ]; then
    echo "== cleaning: deleting every output NOT marked protected_output"
    oxo-flow clean venus.oxoflow --force
    echo "== done. Next 'oxo-flow plan' will show scratch rules as stale —"
    echo "   re-run only what you need ('oxo-flow run' regenerates lazily)."
else
    echo "== PREVIEW (nothing deleted). Re-run with --apply to execute."
    oxo-flow clean venus.oxoflow
fi
