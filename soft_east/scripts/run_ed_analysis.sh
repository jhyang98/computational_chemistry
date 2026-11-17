#!/usr/bin/env bash
# Hypothesis analysis of the ED peak tables with east_edan (T = 1).
# Uses every table below that exists; a point present in several tables is
# taken from the first one listed.
#
# Run from the project root:
#   bash scripts/run_ed_analysis.sh
set -euo pipefail

# ---- parameters ------------------------------------------------------------
TABLES="results/peak_T1.dat results/peak_T1_fine.dat results/peak_T1_L16.dat results/peak_T1_L18.dat"
STRICT="results/strict_T1.dat results/strict_T1_L16_18.dat"
ALPHA="0.25,0.5,1,1.75" # thresholds on alpha = dln chi_KK / dln L for 'cross'
OUT=results/analysis
# ------------------------------------------------------------------------------

USE=""
for t in $TABLES; do
    if [ -f "$t" ]; then USE="$USE $t"; else echo "note: $t not found, skipped" >&2; fi
done
if [ -z "$USE" ]; then
    echo "error: no peak tables found" >&2
    exit 1
fi
mkdir -p "$OUT"
./bin/east_edan h2 $USE --out "$OUT/h2_T1.dat"
./bin/east_edan h1 $USE --out "$OUT/h1_T1.dat"
./bin/east_edan cross --alpha "$ALPHA" $USE --out "$OUT/cross_T1.dat"
SARG=""
for t in $STRICT; do
    if [ -f "$t" ]; then SARG="$SARG --strict $t"; fi
done
./bin/east_edan h3 $SARG $USE --out "$OUT/h3_T1.dat"
echo "written: $OUT/h2_T1.dat $OUT/h1_T1.dat $OUT/cross_T1.dat $OUT/h3_T1.dat"
