#!/usr/bin/env bash
# Reproducible launch script for the cross-method consistency run.
#
# Builds the Release binary via the mingw-release CMake preset, runs
# `--cross-check` with the diameter auto-ranging defaults and the
# `--methods 3,4,5` thesis-target filter, then post-processes the CSV
# into the agreement matrix / disagreements / perf summary artefacts.
#
# Output goes to ./<run_dir>/ (default: experiments/$(date +%F)/).
# Override via:
#   RUN_DIR=experiments/2026-04-30 scripts/run_cross_check.sh
#   METHODS=3,4,5 D_COUNT=4 D_STEP=1 scripts/run_cross_check.sh

set -euo pipefail

# Resolve repo root (this script lives in scripts/).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

RUN_DIR="${RUN_DIR:-experiments/$(date +%F)}"
METHODS="${METHODS:-3,4,5}"
D_COUNT="${D_COUNT:-4}"
D_STEP="${D_STEP:-1}"
EXTRA_ARGS=("$@")

EXE_REL="build/graph_reliability.exe"
[ -x "$EXE_REL" ] || EXE_REL="build/graph_reliability"

echo "[run_cross_check] Building Release via mingw-release preset..."
cmake --preset mingw-release >/dev/null
cmake --build --preset release

echo "[run_cross_check] ctest sanity check..."
ctest --test-dir build --output-on-failure

mkdir -p "$RUN_DIR"
CSV="$RUN_DIR/cross_check.csv"
LOG="$RUN_DIR/cross_check.log"

echo "[run_cross_check] Cross-check → $CSV (methods=$METHODS, d_count=$D_COUNT, d_step=$D_STEP)"
"$EXE_REL" --cross-check \
    --methods "$METHODS" \
    --d-count "$D_COUNT" \
    --d-step  "$D_STEP" \
    --output  "$CSV" \
    "${EXTRA_ARGS[@]}" \
    2>&1 | tee "$LOG"

echo "[run_cross_check] Post-processing $CSV..."
python scripts/cross_check_analyze.py "$CSV"

echo "[run_cross_check] Done. Artefacts in $RUN_DIR/"
ls -la "$RUN_DIR"
