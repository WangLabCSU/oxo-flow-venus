# Environment policy — pixi.lock, installers, and fallback

## Is pixi mandatory for pipeline users?

**No.** The lock files (`envs/<name>/pixi.lock`) are the production
reproducibility reference, not an exclusivity gate:

| User | Path | Result |
|---|---|---|
| Production / engine runs | `pixi install` (respects the lock) | byte-identical to the validated environment |
| Personal exploration | any installer (conda/mamba/pip) against `pixi.toml` version ranges | approximately identical; the lock doubles as the exact-version reference table |

The engine currently drives production runs through pixi, so batch results
are always lock-backed. Nothing prevents a user from building the same
environment by hand; the lock simply tells them exactly what "the same"
means.

## What if a locked environment cannot be restored?

Policy: **strict by default, loud degradation, human arbitration**.

1. Try the lock first (`pixi install`).
2. If the lock cannot be satisfied (packages pulled from channels, platform
   mismatch, solver break), fall back to a fresh solve from `pixi.toml` —
   but:
   - write the new lock as `pixi.lock.new` (never silently overwrite the
     committed `pixi.lock`);
   - treat the environment fingerprint as CHANGED (runs before/after remain
     distinguishable);
   - record the diff (which packages drifted) in the run log and disclose
     "environment reproducibility: degraded" in the report methods section;
   - arbitration (commit the new lock as the new gold standard or revert) is
     a human decision.
3. Engine-side graceful fallback has landed upstream (Traitome/oxo-flow
   #1011, released in v0.29.0; production now runs v0.30.0): when the lock
   cannot be satisfied, the engine itself falls back to a fresh solve and
   records the degradation, instead of failing the run. The human-arbitration
   procedure above still applies to the lock files afterwards: exercise of
   the fallback MUST be recorded in the workdir notes, and committing the
   new lock as gold standard (or reverting) remains a human decision.

## Adding or changing an environment

1. Edit `envs/<name>/pixi.toml` (version ranges or exact pins).
2. Regenerate the lock **on the target platform** (server): `pixi lock
   --manifest-path envs/<name>/pixi.toml`.
3. Commit both files. The engine fingerprints toml+lock — a change invalidates
   the affected rules on the next run (intended).
