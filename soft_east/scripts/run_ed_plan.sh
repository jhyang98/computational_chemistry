#!/usr/bin/env bash
# Continuous-time ED runs for the interim-report plan.
# Run from the project root:  bash scripts/run_ed_plan.sh
set -euo pipefail

mkdir -p results

# 0. validation suite
./bin/east_test > results/test_log.txt
tail -n 1 results/test_log.txt

# 1. strict East coexistence line (eps = 0), crossing with the empty state
./bin/east_ed strict -L 8,10,12,14 --g -0.15:0.05:0.15 --out results/strict_T1.dat

# 2. softened East pseudo-critical points s*(g, eps) = argmax chi_KK
./bin/east_ed peak -L 8,10,12,14 --eps 0.0003,0.001,0.003,0.01,0.03 \
    --g -0.15:0.05:0.15 --s -0.5:0.005:0.8 --out results/peak_T1.dat

# 3. Elmatad-Jack comparison at g = 0, T = 4/3 (J/T = 0.75)
./bin/east_ed peak -L 8,10,12,14 --T 4/3 --eps 0.0003,0.001,0.003,0.01,0.03 \
    --g 0 --s -0.1:0.002:0.2 --out results/peak_T4over3_g0.dat

echo "done: results/strict_T1.dat results/peak_T1.dat results/peak_T4over3_g0.dat"
