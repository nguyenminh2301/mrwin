#!/usr/bin/env bash
#
# run_all.sh — End-to-end orchestrator for the mrwin (P1 v5.2) simulation engine.
#
# Usage:
#   bash run_all.sh replication       # Reproduces the audited 1-CPU reference run (~30-45 min)
#   bash run_all.sh publication       # Canonical IJE-submission scale (~2-3 hr on 32 cores)
#
# Output: 8 JSON files in the current directory:
#   results_packageA_R1.json   results_packageB_R3.json   results_packageC_item1.json
#   results_packageA_R2.json   results_packageB_R4.json   results_packageC_item2.json
#   results_aziz_benchmark.json
#   results_table4.json (master engine, replication mode runs Table 4 only)
#
#   results_packageA_R2.json
#   results_packageB_R3.json
#   results_packageB_R4.json
#   results_packageC_item1.json
#   results_packageC_item2.json
#
set -euo pipefail

MODE=${1:-replication}
if [[ "$MODE" != "replication" && "$MODE" != "publication" ]]; then
    echo "Usage: bash run_all.sh {replication|publication}"
    exit 1
fi

# Optional: override n_proc via env var. Default 1; set higher on cluster.
N_PROC=${N_PROC:-1}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$SCRIPT_DIR/python${PYTHONPATH:+:$PYTHONPATH}"

OUTDIR=$(pwd)
echo "================================================================"
echo "  P1 v5.1 Simulation Engine — running $MODE mode"
echo "  Output directory: $OUTDIR"
echo "  N_PROC for parallelisation across scenarios: $N_PROC"
echo "================================================================"
echo

t_total_start=$(date +%s)

# ----------------------------------------------------------------------------
# Package A R1 — Pleiotropy-bounded CI (analytical, instant)
# ----------------------------------------------------------------------------
echo "[1/8] Package A R1 — pleiotropy-bounded CI"
python -m p1_engine_v5.pleiotropy_bounded_ci
mv -f results_packageA_R1.json "$OUTDIR/" 2>/dev/null || true
echo

# ----------------------------------------------------------------------------
# Package A R2 — InSIDE robustness (6 scenarios)
# ----------------------------------------------------------------------------
echo "[2/8] Package A R2 — InSIDE robustness"
python -m p1_engine_v5.engine \
    --mode "$MODE" \
    --n_proc "$N_PROC" \
    --out "$OUTDIR/results_packageA_R2.json"
echo

# ----------------------------------------------------------------------------
# Package B R3 — Cox-vs-Aalen scale reconciliation
# ----------------------------------------------------------------------------
echo "[3/8] Package B R3 — Cox-vs-Aalen scale reconciliation"
python -m p1_engine_v5.engine_packageB \
    --task R3 \
    --mode "$MODE" \
    --out_dir "$OUTDIR"
echo

# ----------------------------------------------------------------------------
# Package B R4 — Bootstrap moment diagnostics (3 scenarios)
# ----------------------------------------------------------------------------
echo "[4/8] Package B R4 — bootstrap moment diagnostics"
python -m p1_engine_v5.engine_packageB \
    --task R4 \
    --mode "$MODE" \
    --out_dir "$OUTDIR"
echo

# ----------------------------------------------------------------------------
# Package C Item 1 — Q-statistic asymptotic validation
# ----------------------------------------------------------------------------
echo "[5/8] Package C Item 1 — Q-statistic parametric Gaussian validation"
python -m p1_engine_v5.q_statistic_asymptotics \
    --mode "$MODE" \
    --out "$OUTDIR/results_packageC_item1.json"
echo

# ----------------------------------------------------------------------------
# Package C Item 2 — Wall-clock benchmark + biobank-scale feasibility
# ----------------------------------------------------------------------------
echo "[6/8] Package C Item 2 — wall-clock benchmark"
python -m p1_engine_v5.wallclock_benchmark \
    --mode "$MODE" \
    --out "$OUTDIR/results_packageC_item2.json"
echo

# ----------------------------------------------------------------------------
# Package D — Aziz head-to-head benchmark vs per-component pooled MR
# ----------------------------------------------------------------------------
echo "[7/8] Package D — Aziz benchmark (cCWR vs pooled per-component MR)"
python -m p1_engine_v5.engine_aziz_benchmark \
    --mode "$MODE" \
    --out "$OUTDIR/results_aziz_benchmark.json"
echo

# ----------------------------------------------------------------------------
# Master Tables 1-4 (publication mode runs all 4; replication mode runs Table 4)
# ----------------------------------------------------------------------------
echo "[8/8] Master engine — Tables 1-4"
if [[ "$MODE" == "publication" ]]; then
    python -m p1_engine_v5.P1_simulation_engine_v5 \
        --mode publication --tables 1,2,3,4 --out_dir "$OUTDIR"
else
    # Replication mode: Table 4 only (~5 min); Tables 1-3 take longer at small N
    python -m p1_engine_v5.P1_simulation_engine_v5 \
        --mode replication --tables 4 --out_dir "$OUTDIR"
    echo "  [info] Replication mode runs Table 4 only. For full Tables 1-3,"
    echo "         use: python -m p1_engine_v5.P1_simulation_engine_v5 --mode replication --tables 1,2,3"
fi
echo

t_total_end=$(date +%s)
elapsed=$((t_total_end - t_total_start))

echo "================================================================"
echo "  All 8 phases completed in ${elapsed} seconds (mode=$MODE)"
echo "  Outputs in: $OUTDIR"
echo "================================================================"

ls -la "$OUTDIR"/results_*.json 2>/dev/null || echo "WARNING: some output files missing"

# Compute checksums for reproducibility audit
echo
echo "MD5 checksums (compare to REPLICATION_CHECKSUMS.txt for replication mode):"
md5sum "$OUTDIR"/results_*.json 2>/dev/null || true
