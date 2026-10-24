#!/bin/bash
# Bootstrap materialized pixi environments for every venus module env.
#
# WHY: the repo ships only envs/<module>/pixi.toml + pixi.lock (~1MB of
# code). Running the pipeline needs each module env MATERIALIZED
# (envs/<module>/.pixi/envs/default). This script creates those
# materializations idempotently:
#
#   - fresh clone / fresh workdir  -> run once before `oxo-flow run`
#   - after any pixi.toml/pixi.lock change -> re-run; pixi re-syncs only
#     what drifted (this is also the #23 guard: the lock is the single
#     source of truth — never pip-install into a module env at rule
#     runtime, that desyncs the env from the lock and the engine's env
#     fingerprint will re-invalidate the rule on every run)
#
# Usage:
#   scripts/bootstrap_envs.sh            # all module envs
#   scripts/bootstrap_envs.sh amplicon   # one module (prefix match)
#
# Requires pixi on PATH (https://pixi.sh).
set -euo pipefail

here=$(cd "$(dirname "$0")/.." && pwd)
cd "$here"

want=${1:-}
shopt -s nullglob
manifests=(envs/*/pixi.toml)
if [ ${#manifests[@]} -eq 0 ]; then
    echo "error: no envs/*/pixi.toml found (run from the repo root)" >&2
    exit 1
fi

for manifest in "${manifests[@]}"; do
    mod=$(basename "$(dirname "$manifest")")
    if [ -n "$want" ] && [[ "$mod" != "$want"* ]]; then
        continue
    fi
    if [ -d "envs/$mod/.pixi/envs/default" ]; then
        echo "[bootstrap_envs] $mod: env exists, syncing against lock"
    else
        echo "[bootstrap_envs] $mod: materializing env"
    fi
    pixi install --manifest-path "$manifest"
done

echo "[bootstrap_envs] done: all module envs match their pixi.lock files"
